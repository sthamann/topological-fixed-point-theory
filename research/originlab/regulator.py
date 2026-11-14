"""Certified periodic heat trace control with independently bounded errors."""
from flint import arb,ctx
from fractions import Fraction


def certify(grid, cutoff, time="1/2"):
    # Central-difference eigenvalues are compared modewise to the exact
    # continuum Laplacian on the unit periodic circle.
    t = Fraction(time)
    if t <= 0 or cutoff < 0 or grid <= 2*cutoff+1:
        raise ValueError("positive time and cutoff strictly below Nyquist required")
    with ctx.workprec(160):
        tau = arb(str(t))
        pi = arb.pi()
        def weight(mode):
            return (-4*pi*pi*tau*mode**2).exp()

        continuum_low = arb(1)
        discrete_low = arb(1)
        error_low = arb(0)
        for k in range(1,cutoff+1):
            continuous = weight(k)
            discrete = (-4*tau*grid**2*(pi*k/grid).sin()**2).exp()
            continuum_low += 2*continuous
            discrete_low += 2*discrete
            error_low += 2*abs(discrete-continuous)
        # For k>=m, k²>=m²+(k-m)(2m+1).
        m = cutoff+1
        tail = 2*weight(m)/(1-weight(1)**(2*m+1))
        # Finite lattice tail, counting the Nyquist mode only once.
        discrete_tail = arb(0)
        for k in range(m,grid//2+1):
            multiplicity = 1 if 2*k==grid else 2
            discrete_tail += multiplicity*(-4*tau*grid**2*(pi*k/grid).sin()**2).exp()
        bound = error_low+tail+discrete_tail
        # Uniform lattice domination: sin(pi*k/N) >= 2*k/N on the
        # first half-period. Together with pointwise convergence this gives
        # convergence of the whole heat trace for each fixed t>0.
        lattice_uniform_tail = 2*(-16*tau*m*m).exp()/(1-(-16*tau*(2*m+1)).exp())
        return {"status":"CERTIFIED_ENCLOSURE","grid":grid,"cutoff":cutoff,"time":str(t),
                "continuum_trace_lower":str(continuum_low.lower()),
                "continuum_trace_upper":str((continuum_low+tail).upper()),
                "lattice_trace":str(discrete_low+discrete_tail),
                "comparison_error_upper":str(bound.upper()),
                "continuum_tail_upper":str(tail.upper()),
                "lattice_uniform_tail_upper":str(lattice_uniform_tail.upper()),
                "convergence_certificate":"Pointwise symbol limit, dominated by exp(-16*t*k²) for |k|<=N/2.",
                "scope":"free periodic Laplacian, positive supplied time; physical R113 regulator remains separate"}
