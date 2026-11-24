# H6 common-parent synthesis — 24 November 2026

Status: exploratory, scoped results. Three original standalone audit scripts and their original outputs are archived here byte-for-byte. All three scripts were rerun on 24 November with **40/40, 187/187, 201/201** grouped checks and output identical to the archived text files. The total is 428 check calls, not 428 independent proofs. No ledger, frozen prediction, paper, or benchmark promotion. The H6h31 native reconstruction supersedes earlier selection arguments in its own scope; this package complements it and does not replace it.

## Source and scope

- Source audit: NB4 stress measure, positive ladder, radial star, scale-free inverse moments, bounded modular dilation exclusion and frozen TM1 incompatibility. NB4 is a conditional radial stress readout; no four independent spacetime translations follow. Under the specified radial/KMS dictionary r=exp(-4pi), massive/zero weights are 3.48735451781e-6 and 1.39494910413e-5, inverse ratios 286751.563136 and 71688.3595345. The source audit was calculated before H6h31; its k=1 candidate agrees structurally with the native current ladder, but k=4 is not a replacement for that current channel.
- Native half-space: A_q(k)=(S†-w(k)I)(I-qS^4), w=sum_i(1-cos k_i), preserves A_q at zero momentum. The bounded invertible factor R does not change the half-space kernel condition. The unweighted free parent has one massless Weyl corner; finite sections have a counter-boundary. Nonchiral free bulk bound 1-q follows from the Wilson identity. Weighted finite-section tests preserve the kernel; a full domain/continuum construction is not established.
- Sheet reversal does not enforce equal massless and massive spectral weight. The canonical thermofield state is swap invariant for every q. This is a counterexample to that inference, not a global exclusion of q=1/sqrt(2).
- FC, inverse Coxeter compression and FJ generate M3(C); exact TM1 cannot carry both frozen TFPT angles. Free clock weights cannot establish a unique physical flavour selection.
- Common parent: the four-arm H6 star has an explicit layerwise orthogonal intertwiner BU=VA_q and weighted square identity; symbolic checks at depths 1–4, with the same pattern at arbitrary finite depth. This is linear and at zero momentum.
- Local mirror: in the stated z56 Casimir model with the symmetric singlet and particle-hole readout, the gap is 29J/2 and G(z)=z/(z²-(29J/2)²). The isolated local self-energy Delta²/z is **not** the full interacting lattice self-energy. The finite atomic embedding retains all eight corner poles globally; unwanted pole residues decay on the physical left local algebra as depth grows. No finite-J interacting mirror-gap theorem.
- Projective lift: the SU(3)-preserving SU(4) cocharacter (1/2,1/2,1/2,-3/2) compensates the Spin(10) endpoint -1 in the diagonal Z4 E8 quotient. Integral charges of the 240 roots, multiplicities 14,64,84,64,14 for charges -2,-1,0,1,2. This closes the specified clutching construction, not identification of the native physical bundle. The naive family-line Dirac index is 8 mod 16 under the stated U(5) Fock-bundle assumptions, not a full family of 16. Branching 4=1+3 alone does not eliminate the fourth component.
- Clifford area: the oriented gamma5 trace yields the z57 t=1 Plebanski potential. The trace is indefinite, so positivity alone cannot select it. Curved gauge backgrounds yield Pauli and mixed Wilson commutator terms: the square is block diagonal only in the flat commuting case.

## Additional exact inequality

For a normalized positive Stieltjes response C(x)=1/x+sum_i w_i/(x+m_i²) without contact terms, w_i>=0, m_i²>0 and finite W=sum_i w_i>0, a=sum_i w_i/m_i², b=a²+sum_i w_i/m_i^4. Cauchy–Schwarz gives a²<=W sum_i w_i/m_i^4, hence b/a²>=1+1/W if a>0. Equality holds precisely when the massive spectral measure is supported at one mass. This strengthens the existing v1037 sign test in this class, independent of scale. This is an analytic observation, not a new numerical check count.

## Replay

From this directory with numpy, scipy, sympy installed:
```sh
python tfpt_h6h30_source_audit.py
python tfpt_native_halfspace.py
python tfpt_h6_common_parent.py
```

The filenames are the original filenames. For provenance, the files were compared by Git blob identity against the original conversation artifacts. The archive is not a new OriginLab lane. All physical source identifications, interacting gauge measure, four-dimensional metric soldering, nonperturbative dynamics, flavour selection and dimensionless scale ratios remain open unless resolved by a separate source-reviewed result. H6 and H6h31 remain the baseline; no TOE declaration.
