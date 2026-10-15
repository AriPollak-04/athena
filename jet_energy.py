# %%
"""Jet energy bookkeeping and input-file generation for the ``jet_blast`` problem.

The problem generator drives the jet with a luminosity

    L(t) = L_peak * a * exp(-t / t_0)

carried by the Lorentz factor.  Because L is proportional to g(Gamma) = Gamma^2 v with a
time-independent prefactor, imposing that envelope means solving g(Gamma(t)) = a
exp(-t/t_0) g(jet_Gam), which inverts in closed form (see `gamma_of_g`).  Three relations
follow, and they are the whole design:

    Gamma(0) = jet_Gam * sqrt(a)                  peak Lorentz factor
    E        = L_peak * a * t_0                   total energy   (t_jet = a * t_0)
    t_stop   = t_0 * ln( a g(jet_Gam)/g(Gam_end) ) when Gamma has decayed to Gam_end

So `a` sets the peak, `t_0` sets how fast it dies, and the energy follows.  Note that
t_jet = E / L_peak is a *definition*, not a separate input: fixing the energy fixes t_jet
and vice versa.  That leaves two free knobs, which is why `solve_exp` wants exactly two.

``jet_turnoff = enthalpy`` in the athinput (``turnoff="enthalpy"`` here) carries the same
envelope by h = 1 + 4p/rho instead, at fixed rho and Gamma = jet_Gam.  L is linear in
rho h, so h(t) = a exp(-t/t_0) h_ref gives the identical L(t) and E with no inversion.
The catch is that h cannot drop below h_ref (p below jet_p), so the decay is exhausted at
t_0 ln(a), where L has only fallen by 1/a = h_ref/h(0).  The thermal part is all it can
remove; the cold Gamma = jet_Gam beam that remains is cut off as a square pulse.

Typical use, straight from an athinput's <problem> block::

    import jet_energy as je
    je.solve_exp(a=2.0, t0=0.4875, jet_rho=1.25204135e-03)
    je.solve_exp(E=2.0e-3, Gamma0=44.0, jet_rho=1.25204135e-03)
    je.solve_exp(turnoff="enthalpy", a=2.94, t0=0.8, jet_rho=1.451e-03)

Put the `a`, `t0` and `t_stop` it returns into the athinput and the run will reproduce
the energy it reports.
"""

import math
import re

GAMMA_EOS = 4.0 / 3.0


def g_of_gamma(gam):
    """g(Gamma) = Gamma^2 v = Gamma sqrt(Gamma^2 - 1); increasing for Gamma >= 1."""
    return 0.0 if gam <= 1.0 else gam * math.sqrt(gam * gam - 1.0)


def gamma_of_g(gval):
    """Inverse of `g_of_gamma`: x = Gamma^2 solves x^2 - x - g^2 = 0, positive root."""
    if gval <= 0.0:
        return 1.0
    return math.sqrt(0.5 * (1.0 + math.sqrt(1.0 + 4.0 * gval * gval)))


def peak_luminosity(jet_rho, jet_p, jet_Gam, jet_rinj, theta_0):
    """Full-strength luminosity: both jets through the nozzle cross section.

    L_peak = 2 (rho + 4p) Gamma^2 v pi r_inj^2 sin^2(theta_0); the 4p is rho*h - rho
    for gamma_ad = 4/3.
    """
    w = jet_rho + (GAMMA_EOS / (GAMMA_EOS - 1.0)) * jet_p
    return (2.0 * w * g_of_gamma(jet_Gam) * math.pi
            * jet_rinj ** 2 * math.sin(theta_0) ** 2)


def solve_exp(jet_Gam=31.0, jet_rho=7.85e-4, jet_p=1e-6, jet_rinj=0.1, theta_0=0.174533,
              Gam_end=1.5, E=None, t_jet=None, a=None, Gamma0=None, t0=None,
              turnoff="gamma", h0=None, h_end=None):
    """Close the design from exactly two of {E or t_jet, a or Gamma0 or h0, t0}.

    `E` is E/Mc^2 for a unit stellar mass.  `Gam_end` is the Lorentz factor the jet has
    decayed to by `t_stop`; it only sets where to stop stamping, not the energy.

    `turnoff="enthalpy"` matches ``jet_turnoff = enthalpy``: Gamma stays at jet_Gam, the
    peak is given as `h0` = a h_ref instead of `Gamma0`, and `h_end` (default h_ref, the
    floor) replaces `Gam_end`.  `Gam_end` is ignored in that mode.

    Returns a dict of a, t0, Gamma0, h0, t_jet, t_stop, E, L_peak and L_stop_frac, the
    last being L(t_stop)/L(0), the step the drive ends on.
    """
    if turnoff not in ("gamma", "enthalpy"):
        raise ValueError("turnoff must be 'gamma' or 'enthalpy', not %r." % turnoff)
    by_h = turnoff == "enthalpy"
    if E is not None and t_jet is not None:
        raise ValueError("give E or t_jet, not both -- t_jet = E / L_peak.")
    if sum(x is not None for x in (a, Gamma0, h0)) > 1:
        raise ValueError("give one of a, Gamma0, h0 -- each fixes the other two.")
    if by_h and Gamma0 is not None:
        raise ValueError("Gamma is fixed at jet_Gam for turnoff='enthalpy'; give h0.")
    if not by_h and (h0 is not None or h_end is not None):
        raise ValueError("h0 and h_end only apply to turnoff='enthalpy'.")

    L_peak = peak_luminosity(jet_rho, jet_p, jet_Gam, jet_rinj, theta_0)
    h_ref = 1.0 + (GAMMA_EOS / (GAMMA_EOS - 1.0)) * jet_p / jet_rho

    if by_h:
        h_end = h_ref if h_end is None else h_end
        if h_end < h_ref:
            raise ValueError("h_end = %g is below h_ref = %g: p would drop under jet_p."
                             % (h_end, h_ref))
        r_end = h_end / h_ref
    elif Gam_end <= 1.0:
        raise ValueError("Gam_end must exceed 1: Gamma reaches 1 only asymptotically.")
    else:
        r_end = g_of_gamma(Gam_end) / g_of_gamma(jet_Gam)

    if E is not None:
        t_jet = E / L_peak
    if Gamma0 is not None:
        a = g_of_gamma(Gamma0) / g_of_gamma(jet_Gam)
    if h0 is not None:
        a = h0 / h_ref

    if sum(x is not None for x in (t_jet, a, t0)) != 2:
        raise ValueError("need exactly two of {E or t_jet, a or Gamma0 or h0, t0}.")
    if a is not None and a <= r_end:
        raise ValueError("a = %g is below r_end = %g: the jet starts under its end state."
                         % (a, r_end))

    if t_jet is None:
        t_jet = t0 * (a - r_end)        # E = L_peak a t_0, exactly, given t_stop below
    elif t0 is None:
        t0 = t_jet / (a - r_end)
    else:
        a = t_jet / t0 + r_end

    if t0 <= 0.0:
        raise ValueError("solved t_0 = %g is not positive; check the inputs." % t0)

    return {
        "a": a,
        "t0": t0,
        "t_stop": t0 * math.log(a / r_end),
        "Gamma0": jet_Gam if by_h else gamma_of_g(a * g_of_gamma(jet_Gam)),
        "h0": a * h_ref if by_h else h_ref,
        "t_jet": t_jet,
        "E": L_peak * t_jet,
        "L_peak": L_peak,
        "L_stop_frac": r_end / a,
    }


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
