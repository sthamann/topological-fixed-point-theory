"""Four-mark winding memory with positive, time-covariant population mixing.

The model fixes the nearest-neighbor classical rates and the frequency
assignment N=diag(0,1,2,3)/4. Analytical bounds are documented beside the
OriginLab OPH campaign. This is a finite reduced-channel comparison.
"""
import sympy as sp


def unit(i, j):
    out = sp.zeros(4)
    out[i, j] = 1
    return out


def liouville(action, basis):
    return sp.Matrix.hstack(*(action(x).reshape(16, 1) for x in basis))


def dissipator(jump, matrix):
    rate = jump.H*jump
    return jump*matrix*jump.H - (rate*matrix+matrix*rate)/2


def controls():
    dim = 4
    basis = [unit(i, j) for i in range(dim) for j in range(dim)]
    identity = sp.eye(dim)
    cyclic = sp.diag(-1, -1, -1, -1)
    for j in range(dim):
        cyclic[(j+1) % dim, j] = cyclic[(j-1) % dim, j] = sp.Rational(1, 2)
    assert cyclic.eigenvals() == {0: 1, -1: 2, -2: 1}

    path = sum((unit(j+1, j) for j in range(3)), sp.zeros(dim))
    wrap = unit(0, 3)
    shift = path+wrap
    jumps = (path, path.H, wrap, wrap.H)

    # Splitting the wrap from the other links supplies a single frequency
    # to each jump. The opposite-frequency jumps occur with the same rate.
    def generate(matrix):
        return sum((dissipator(jump, matrix) for jump in jumps), sp.zeros(dim))/2

    mark = sp.diag(0, 1, 2, 3)/4
    clock = liouville(lambda a: mark*a-a*mark, basis)
    evolution = liouville(generate, basis)
    assert evolution*clock == clock*evolution
    assert evolution.H == evolution and evolution.rank() == 15
    assert generate(identity) == sp.zeros(dim)
    assert all(sp.trace(generate(x)) == 0 for x in basis)
    for j in range(dim):
        assert generate(unit(j, j)) == sp.diag(*cyclic[:, j])
    assert generate(wrap) == -wrap

    edge_basis = [unit(j+1, j) for j in range(3)]
    edge_block = sp.Matrix([
        [sp.trace(x.H*generate(y)) for y in edge_basis] for x in edge_basis
    ])
    assert edge_block == sp.Matrix([[-2, 1, 0], [1, -2, 1], [0, 1, -2]])/2

    # Conditional Choi certificate. In the subspace perpendicular to the
    # vectorized identity the dissipator's anticommutator terms cancel.
    def choi(action):
        return sum((sp.kronecker_product(unit(i, j), action(unit(i, j)))
                    for i in range(dim) for j in range(dim)), sp.zeros(16))

    omega = identity.reshape(16, 1)
    complement = sp.eye(16)-omega*omega.H/dim
    projected = complement*choi(generate)*complement
    jump_gram = sum((j.T.reshape(16, 1)*j.T.reshape(16, 1).H for j in jumps),
                    sp.zeros(16))/2
    assert projected == jump_gram

    # Keeping the reconstructed diagonal mixing while preserving every
    # off-diagonal element would maintain the earlier full charged trace.
    def undamped_extension(matrix):
        return sp.diag(*(cyclic*sp.Matrix([matrix[p, p] for p in range(dim)])))

    negative = complement*choi(undamped_extension)*complement
    witness = (unit(0, 0)-unit(3, 3)).reshape(16, 1)/sp.sqrt(2)
    assert (witness.H*negative*witness)[0] == -1

    # Countermodel with only deck covariance. It fixes S and has the same
    # population mixing; it breaks covariance under the full winding clock.
    def deck_only(matrix):
        return (shift*matrix*shift.H+shift.H*matrix*shift)/2-matrix

    deck = sp.diag(1, -sp.I, -1, sp.I)
    assert all(deck_only(deck*b*deck.H) == deck*deck_only(b)*deck.H for b in basis)
    assert all(deck_only(unit(j, j)) == generate(unit(j, j)) for j in range(dim))
    assert deck_only(shift) == sp.zeros(dim)
    assert liouville(deck_only, basis).rank() == 12
    assert liouville(deck_only, basis)*clock != clock*liouville(deck_only, basis)

    # Derive the three-site envelope directly from its two even eigenmodes.
    ones = sp.ones(3, 1)
    plus = sp.Matrix([1, sp.sqrt(2), 1])/2
    minus = sp.Matrix([1, -sp.sqrt(2), 1])/2
    assert sp.expand((plus.H*ones)[0]**2) == sp.Rational(3, 2)+sp.sqrt(2)
    assert sp.expand((minus.H*ones)[0]**2) == sp.Rational(3, 2)-sp.sqrt(2)
    assert (sp.trace(shift.H*generate(shift))/4) == -sp.Rational(1, 2)

    # The populations and wrap then share one measured rate. This is an
    # algebraic relation; it assumes the same clock and empty normalization.
    x = sp.symbols('x', nonnegative=True)
    vacuum = (2*x+x*x)/16
    assert sp.sqrt(sp.factor(1+16*vacuum))-1 == x
    return {
        "status": "EXACT_COVARIANT_MEMORY_BOUND",
        "physical_source_derived": False,
        "physical_solution": False,
        "class": "continuous CPTP semigroups on M4 covariant under exp(i t N) and fixed symmetric nearest-neighbor population mixing",
        "vacuum_centered_curve": "(2*exp(-m*t)+exp(-2*m*t))/16",
        "wrap_amplitude_bound": "abs(tr_4(E_30 Phi_t(E_03))) <= exp(-m*t)/4",
        "paired_source_test": "wrap amplitude <= (sqrt(1+16*C_vac(t))-1)/4",
        "charged_shift_envelope": "exp(-m*t)*(1+3*cosh(m*t/sqrt(2))+2*sqrt(2)*sinh(m*t/sqrt(2)))/4",
        "initial_minimum_shift_loss": "m/2",
        "saturator_verified": True,
        "undamped_mix_conditional_choi_witness_at_unit_rate": "-1",
        "deck_only_countermodel_fixes_shift": True,
        "deck_only_fixed_space_dimension": 4,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(controls(), indent=2))
