"""Jet energy bookkeeping and input-file generation for the ``jet_blast`` problem.

The ``jet_blast`` problem generator drives a bipolar jet by stamping a nozzle region
every cycle (``src/pgen/jet_blast.cpp``, ``Mesh::UserWorkInLoop``).  This module answers
the two questions that keep coming up when tuning it:

1.  How much energy does a given ``<problem>`` block actually inject?
2.  What ``jet_rho`` do I need to hit a target energy for some other ``t_stop`` /
    ``t_ramp`` / ``jet_Gam``?

Two energy conventions are reported, because they do not scale the same way:

``E3d``
    The 3D-conical convention, ``2 (rho + 4p) Gam^2 v pi r_inj^2 sin^2(theta_0) t``,
    normalised by the stellar mass ``M = 1``.  This is the number usually quoted as
    ``E/Mc^2``.
``E2d``
    What the code actually injects in a 2D ``(R, phi)`` run (``nx3 = 1``): the flux
    through the two nozzle arcs, ``4 (rho + 4p) Gam^2 v theta_0 r_inj t`` per unit z,
    normalised by the in-plane stellar mass ``\\int rho 2 pi r dr``.

The two agree on ``rho_j <-> t_stop`` and ``rho_j <-> Gam^-2`` rescalings, and disagree
on ``theta_0`` and ``r_inj`` (``r_inj^2 sin^2 theta_0`` versus ``r_inj theta_0``).  Always
check both before changing the nozzle geometry.

Everything is stdlib-only so this can be imported from a notebook or a batch script.
"""

import math
import re

# Defaults matching inputs/mhd/athinput.jet_blast.
GAMMA_EOS = 4.0 / 3.0
DEFAULT_CSV = "polytrope_n3.00.csv"


def read_athinput(path):
    """Parse an Athena++ input file into a nested dict ``{section: {key: value}}``."""
    params = {}
    section = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"^<(\w+)>$", line)
            if m:
                section = m.group(1)
                if section != "comment":
                    params[section] = {}
                continue
            if section and section != "comment" and "=" in line:
                key, _, rest = line.partition("=")
                value = rest.split("#")[0].strip()
                key = key.strip()
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                params[section][key] = value
    return params


def beta_of(gam):
    """3-velocity for a Lorentz factor, clamped at rest for ``gam <= 1``."""
    if gam <= 1.0:
        return 0.0
    return math.sqrt(1.0 - 1.0 / (gam * gam))


def jet_gamma_of_time(t, gam, t_stop, t_ramp=0.0, gam_end=1.0):
    """Instantaneous injected Lorentz factor.

    Mirrors ``JetGammaOfTime`` in ``src/pgen/jet_blast.cpp``: full strength until
    ``t_stop - t_ramp``, then a raised-cosine taper down to ``gam_end`` at ``t_stop``.
    """
    if t > t_stop:
        return 0.0
    if t_ramp <= 0.0:
        return gam
    t_on = t_stop - t_ramp
    if t <= t_on:
        return gam
    s = min(max((t - t_on) / t_ramp, 0.0), 1.0)
    return gam_end + (gam - gam_end) * 0.5 * (1.0 + math.cos(math.pi * s))


def taper_factor(gam, gam_end=1.0, nsub=20000):
    """Mean of ``Gam(s)^2 v(s)`` over the ramp, normalised by its full-strength value.

    This is the fraction of a full-strength second that one second of spin-down is
    worth.  For ``gam = 31``, ``gam_end = 1`` it is 0.38311.
    """
    if gam <= 1.0:
        return 0.0
    full = gam * gam * beta_of(gam)
    acc = 0.0
    for k in range(nsub):
        s = (k + 0.5) / nsub
        g = gam_end + (gam - gam_end) * 0.5 * (1.0 + math.cos(math.pi * s))
        acc += g * g * beta_of(g)
    return (acc / nsub) / full


def effective_duration(t_stop, t_ramp=0.0, gam=31.0, gam_end=1.0):
    """Full-strength-equivalent drive duration, ``t_stop - (1 - I) * t_ramp``."""
    if t_ramp <= 0.0:
        return t_stop
    t_ramp = min(t_ramp, t_stop)
    return (t_stop - t_ramp) + t_ramp * taper_factor(gam, gam_end)


def jet_luminosity(rho_j, p_j, gam, r_inj, theta_0):
    """Return ``(L3d, L2d)``: conical and 2D-planar (per unit z) energy flux."""
    w = rho_j + (GAMMA_EOS / (GAMMA_EOS - 1.0)) * p_j   # rho * h, enthalpy density
    common = w * gam * gam * beta_of(gam)
    l3d = 2.0 * common * math.pi * r_inj ** 2 * math.sin(theta_0) ** 2
    l2d = 4.0 * common * theta_0 * r_inj
    return l3d, l2d


def m2d_from_csv(path=DEFAULT_CSV):
    """In-plane stellar mass ``\\int rho 2 pi r dr`` from the polytrope CSV."""
    r, rho = _load_profile(path)
    total = 0.0
    for i in range(1, len(r)):
        rmid = 0.5 * (r[i - 1] + r[i])
        rhomid = 0.5 * (rho[i - 1] + rho[i])
        total += rhomid * 2.0 * math.pi * rmid * (r[i] - r[i - 1])
    return total


def rho_env(radius, path=DEFAULT_CSV):
    """Linearly interpolated envelope density at ``radius`` (0 outside the star)."""
    r, rho = _load_profile(path)
    if radius <= r[0]:
        return rho[0]
    if radius >= r[-1]:
        return 0.0
    lo, hi = 0, len(r) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if r[mid] <= radius:
            lo = mid
        else:
            hi = mid
    frac = (radius - r[lo]) / (r[hi] - r[lo])
    return rho[lo] + frac * (rho[hi] - rho[lo])


_PROFILE_CACHE = {}


def _load_profile(path):
    if path not in _PROFILE_CACHE:
        r, rho = [], []
        with open(path) as f:
            f.readline()                      # header: r,rho,m,P
            for line in f:
                if not line.strip():
                    continue
                cols = line.split(",")
                r.append(float(cols[0]))
                rho.append(float(cols[1]))
        _PROFILE_CACHE[path] = (r, rho)
    return _PROFILE_CACHE[path]


def jet_energy(params, csv_path=DEFAULT_CSV, m_star=1.0):
    """Energy budget for a ``<problem>`` block (or a full parsed athinput).

    Accepts either the dict returned by :func:`read_athinput` or just its ``problem``
    sub-dict.  Returns a dict with both bookkeepings plus the pieces that go into them.
    """
    prob = params.get("problem", params)
    rho_j = float(prob["jet_rho"])
    p_j = float(prob["jet_p"])
    gam = float(prob["jet_Gam"])
    r_inj = float(prob["jet_rinj"])
    theta_0 = float(prob["theta_0"])
    t_stop = float(prob["t_stop"])
    t_ramp = float(prob.get("t_ramp", 0.0))
    gam_end = float(prob.get("jet_Gam_end", 1.0))

    t_eff = effective_duration(t_stop, t_ramp, gam, gam_end)
    l3d, l2d = jet_luminosity(rho_j, p_j, gam, r_inj, theta_0)
    m2d = m2d_from_csv(csv_path)
    rho_a = rho_env(r_inj, csv_path)
    w = rho_j + (GAMMA_EOS / (GAMMA_EOS - 1.0)) * p_j

    return {
        "h": w / rho_j,
        "beta": beta_of(gam),
        "t_eff": t_eff,
        "taper": taper_factor(gam, gam_end) if t_ramp > 0.0 else 1.0,
        "L3d": l3d,
        "E3d": l3d * t_eff,
        "E3d_over_Mc2": l3d * t_eff / m_star,
        "L2d": l2d,
        "E2d": l2d * t_eff,
        "M2d": m2d,
        "E2d_over_M2d": l2d * t_eff / m2d,
        "L_tilde": w * gam * gam * beta_of(gam) / rho_a,
        "rho_env_at_rinj": rho_a,
    }


def solve_jet_rho(e_target, t_stop, t_ramp=0.0, gam=31.0, gam_end=1.0,
                  p_j=1e-6, r_inj=0.1, theta_0=0.174533, m_star=1.0):
    """``jet_rho`` giving ``E3d/Mc^2 == e_target``.

    Exact, not iterative: the energy is linear in ``rho_j`` once the effective duration
    is known.  Raises ``ValueError`` if the requested energy is unreachable at this
    ``jet_p`` (the enthalpy floor alone already overshoots it).
    """
    t_eff = effective_duration(t_stop, t_ramp, gam, gam_end)
    geom = 2.0 * gam * gam * beta_of(gam) * math.pi * r_inj ** 2 \
        * math.sin(theta_0) ** 2 * t_eff
    w_needed = e_target * m_star / geom
    rho_j = w_needed - (GAMMA_EOS / (GAMMA_EOS - 1.0)) * p_j
    if rho_j <= 0.0:
        raise ValueError(
            "E/Mc^2 = %g is unreachable with jet_p = %g: the pressure term alone "
            "contributes more than the target." % (e_target, p_j))
    return rho_j


def write_variant(base_path, out_path, **sections):
    """Copy an athinput, overriding keys in place, preserving comments and layout.

    Called as ``write_variant(base, out, job={'problem_id': 'jb_late'},
    problem={'t_stop': 4.5, 'jet_rho': 6.0967e-4})``.  Keys absent from the source file
    are appended to the end of their section.  Returns the set of keys that were added
    rather than overwritten.
    """
    pending = {sec: dict(kv) for sec, kv in sections.items()}
    added = set()
    out_lines = []
    section = None
    last_hit = {}          # section -> index in out_lines of its last overridden key

    def flush(sec):
        """Emit any keys for ``sec`` that were not found in the file.

        New keys are placed just after the last key this call overrode in the same
        section, so they land beside related settings rather than at the very end.
        """
        leftover = pending.pop(sec, None)
        if not leftover:
            return
        at = last_hit.get(sec)
        block = ["%-15s = %s\n" % (k, _fmt(v)) for k, v in leftover.items()]
        added.update("%s/%s" % (sec, k) for k in leftover)
        if at is None:
            if out_lines and not out_lines[-1].endswith("\n"):
                out_lines[-1] += "\n"   # source file may lack a trailing newline
            out_lines.extend(block)
        else:
            out_lines[at + 1:at + 1] = block

    with open(base_path) as f:
        lines = f.readlines()

    for line in lines:
        m = re.match(r"^\s*<(\w+)>\s*$", line)
        if m:
            flush(section)
            section = m.group(1)
            out_lines.append(line)
            continue
        stripped = line.strip()
        if section in pending and stripped and not stripped.startswith("#") \
                and "=" in stripped:
            key = stripped.partition("=")[0].strip()
            if key in pending[section]:
                value = pending[section].pop(key)
                comment = ""
                if "#" in line:
                    comment = "  " + line[line.index("#"):].rstrip()
                last_hit[section] = len(out_lines)
                out_lines.append("%-15s = %s%s\n" % (key, _fmt(value), comment))
                continue
        out_lines.append(line)
    flush(section)

    with open(out_path, "w") as f:
        f.writelines(out_lines)
    return added


def _fmt(value):
    """Format a value for an athinput line, keeping small densities in sci notation."""
    if isinstance(value, float):
        if value != 0.0 and abs(value) < 1e-2:
            return "%.8e" % value
        return repr(value) if value == int(value) else "%.6g" % value
    return str(value)


def report(params, label="", csv_path=DEFAULT_CSV):
    """One-line-per-quantity summary; handy from a notebook cell."""
    e = jet_energy(params, csv_path)
    head = "jet energy budget" + (" [%s]" % label if label else "")
    print(head)
    print("  h              = %.7f" % e["h"])
    print("  beta_jet       = %.7f" % e["beta"])
    if e["taper"] < 1.0:
        print("  taper factor   = %.5f  (t_eff = %.4f)" % (e["taper"], e["t_eff"]))
    print("  L_tilde        = %.4e   (rho_env(r_inj) = %.4f)"
          % (e["L_tilde"], e["rho_env_at_rinj"]))
    print("  L_jet (3D)     = %.6e" % e["L3d"])
    print("  E/Mc^2 (3D)    = %.6e   <-- the usual quoted number" % e["E3d_over_Mc2"])
    print("  E2d/M2d c^2    = %.6e   (M2d = %.4f, what the 2D run carries)"
          % (e["E2d_over_M2d"], e["M2d"]))
    return e
