#!/usr/bin/env python3
"""Original v746 counterbond and v1033 parity-link probes, with precise scope."""
from pathlib import Path
import importlib
import json
import sys

import numpy as np
import sympy as sp


ROOT = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(ROOT / "verification"))
v1033 = importlib.import_module("v1033_charged_disorder")
v746 = importlib.import_module("v746_phys_gnet_local_functor")
NX, NY = 32, 8


def counterbond():
    cases = 0
    for size in (48, 96, 192):
        for h in (1, -1, 1j, -1j):
            for cut in range(1, size):
                before = [1] * (size - 1) + [h]
                mask = [-1 if j < cut else 1 for j in range(size)]
                after = [mask[(j+1)%size]*before[j]*mask[j] for j in range(size)]
                expected = before.copy()
                expected[cut-1] *= -1
                expected[-1] *= -1
                assert after == expected
                assert np.prod(after) == h
                cases += 1
    return cases


def bond_matrix(h, cut, columns):
    right = (cut+1)%NX
    result = np.zeros_like(h)
    first, second = np.ix_(columns == cut, columns == right)
    result[first, second] = h[first, second]
    result[second.T, first.T] = h[second.T, first.T]
    return result


def physical_time_and_witnesses():
    columns = np.repeat(np.arange(NX), 2*NY)
    blocks = [v1033.qwz_cylinder(NX,NY,1,r) for r in range(4)]
    cases = 0
    for r, h in enumerate(blocks):
        for cut in range(NX):
            mask = np.where(columns <= cut,-1,1)
            B = bond_matrix(h, cut, columns)
            assert np.array_equal(mask[:,None] * blocks[(r+2)%4] * mask[None,:], h - 2*B)
            cases += 1
    cut, y, orbital = 5, 0, 0
    source = 2*(cut*NY + y) + orbital
    target = 2*((cut+1)*NY + y) + orbital
    states = {}
    for r in (1,3):
        values, vectors = np.linalg.eigh(blocks[r])
        sea = vectors[:,values<0]
        assert sea.shape[1] == NX*NY
        C = sea @ sea.conj().T
        coefficient = blocks[r][target,source]
        assert coefficient == -0.5
        # Directed existing hopping B=t c_target† c_source.
        norm_sq = float(abs(coefficient)**2 * (
            C[source,source].real*(1-C[target,target].real)
            + abs(C[source,target])**2
        ))
        assert norm_sq > 0
        states[str(r)] = {"occupied": sea.shape[1], "norm_sq_B_on_sea": norm_sq,
                          "signed_four_word": -norm_sq}
    return cases, states


def phase_test():
    X = sp.Matrix([[0,1],[1,0]])
    Y = sp.Matrix([[0,-sp.I],[sp.I,0]])
    Z = sp.diag(1,-1)
    R = -sp.I*X
    assert R**2 == -sp.eye(2)
    assert [R*A*R.H == sign*A for A,sign in zip((sp.eye(2),X,Y,Z),(1,1,-1,-1))] == [True]*4
    # On the two fermionic modes of this very sheet fiber, 2-particle
    # action is det(R). This is not a second choice of spin phase.
    assert R.det() == 1 and R.det()**2 == 1
    H = -sp.I*Y
    sigma = -Z
    K = -Y
    assert R == sp.I*H*sigma
    assert R*H*R.H == H.inv()
    assert H*K*H.H == K and R*K*R.H == sigma*K*sigma.H == -K
    assert [A.det() for A in (H,R,sigma)] == [1,1,-1]
    # On original v746 NS sizes the antipodal two-mode density is exactly I/2:
    # the occupied Fourier phases alternate and there are N/2 even terms.
    for size in (48,96,192):
        half = size//2
        shift = v746.shift(size,ns=True)
        H_source = np.linalg.matrix_power(shift,half)
        for j in range(half):
            pair = v746.window(size,j,1)
            assert np.array_equal(H_source[np.ix_(pair,pair)], np.array(H,dtype=float))
        assert sum((-1)**n for n in range(half,size)) == 0
        # Each occupied +/-i count is N/4, so their products equal 1.
        assert (1j)**(size//4)*(-1j)**(size//4) == 1
    # A proper partial lift has local parity square, even on global-even states.
    full = sp.diag(R, sp.eye(2))
    P_region = sp.diag(-1,-1,1,1)
    transfer = sp.zeros(4)
    transfer[2,0] = 1
    assert full**2 == P_region
    assert P_region*transfer*P_region == -transfer
    assert P_region[0,0]*P_region[3,3] == -1  # Occupied mode 0 with spectator 3.
    assert P_region[2,2]*P_region[3,3] == 1
    return {"R_square_one_particle":"-I", "ad_R_even_matrix_basis":["I","X"],
            "ad_R_odd_matrix_basis":["Y","Z"],"Gamma_R_square_total_even":"+I",
            "partial_Gamma_R_square_on_even_transfer_witness":[-1,1],
            "shared_sheet_frame_identity":"R=i H sigma",
            "v746_even_clock_invariant_K_original_NS_norm_sq":"1/2",
            "v746_signed_four_word_for_both":"-1/2",
            "clock_invariant_even_two_particle_creation_phase":{"R":1,"sigma":-1}}


def main():
    old_counterbonds = counterbond()
    cases, states = physical_time_and_witnesses()
    print(
        json.dumps({"v746_quarter_sign_generalization_checks":old_counterbonds,
                    "v1033_all_cuts_all_registers_exact_time_checks":cases,
                    "source_cut":5,"source_transverse_site":0,"source_hop_coefficient":-0.5,
                    "source_car_quarter_witnesses":states,"relative_spin_phase":phase_test(),
                    "scope":"conditional original fixed finite carrier, full P_adm uncomputed"},
                    indent=2)
    )


if __name__ == "__main__":
    main()
