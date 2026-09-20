"""Solve the exponential jet turn-off for `jet_blast`: give any two knobs, get the rest.

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

Typical use, straight from an athinput's <problem> block::

    import jet_energy as je
    je.solve_exp(a=2.0, t0=0.4875, jet_rho=1.25204135e-03)
    je.solve_exp(E=2.0e-3, Gamma0=44.0, jet_rho=1.25204135e-03)

Put the `a`, `t0` and `t_stop` it returns into the athinput and the run will reproduce
the energy it reports.
"""

import math

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
              Gam_end=1.5, E=None, t_jet=None, a=None, Gamma0=None, t0=None):
    """Close the design from exactly two of {E or t_jet, a or Gamma0, t0}.

    `E` is E/Mc^2 for a unit stellar mass.  `Gam_end` is the Lorentz factor the jet has
    decayed to by `t_stop`; it only sets where to stop stamping, not the energy.

    Returns a dict of a, t0, Gamma0, t_jet, t_stop, E and L_peak.
    """
    if E is not None and t_jet is not None:
        raise ValueError("give E or t_jet, not both -- t_jet = E / L_peak.")
    if a is not None and Gamma0 is not None:
        raise ValueError("give a or Gamma0, not both -- Gamma0 = jet_Gam sqrt(a).")
    if Gam_end <= 1.0:
        raise ValueError("Gam_end must exceed 1: Gamma reaches 1 only asymptotically.")

    L_peak = peak_luminosity(jet_rho, jet_p, jet_Gam, jet_rinj, theta_0)
    r_end = g_of_gamma(Gam_end) / g_of_gamma(jet_Gam)

    if E is not None:
        t_jet = E / L_peak
    if Gamma0 is not None:
        a = g_of_gamma(Gamma0) / g_of_gamma(jet_Gam)

    if sum(x is not None for x in (t_jet, a, t0)) != 2:
        raise ValueError("need exactly two of {E or t_jet, a or Gamma0, t0}.")
    if a is not None and a <= r_end:
        raise ValueError("a = %g is below r_end = %g: Gamma(0) is already under Gam_end."
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
        "t_stop": t0 * math.log(a * g_of_gamma(jet_Gam) / g_of_gamma(Gam_end)),
        "Gamma0": gamma_of_g(a * g_of_gamma(jet_Gam)),
        "t_jet": t_jet,
        "E": L_peak * t_jet,
        "L_peak": L_peak,
    }
