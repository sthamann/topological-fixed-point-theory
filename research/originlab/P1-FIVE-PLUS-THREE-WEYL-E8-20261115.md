# Exact F1 algebraic bridge: P1 geometric 5+3 quarter clock is native E8

15 November 2026. Exact E8 root-lattice theorem, not a proof of the physical F1 origin. Native weighted R113 values, original fermion realization, physical time and 3+1D parent remain open.

## Original owners

The actual original P1 and P2 geometry is in documents/tfpt_1_architecture_e8/sections/013-c3-p1-p2-origin.tex, verification/v189_riemann_roch_carrier.py, v197_rr_carrier_clifford_d5.py, and v512_seam_tau_flag.py. The native E8 D8+ roots and glue grade are defined in research/originlab/affine_operators.py and marked_source.py. The independent exact execution and root census are in p1_e8_geometric_lift.py and test_p1_e8_geometric_lift.py.

## Eight complex geometric mode characters

Let D=mu4 and rotate z -> i*z.

C5=H0(P1,O(D)) has basis f_k=z^k/(z^4-1), k=0,...,4, phases 1,i,-1,-i,1.
F3=H0(P1,Omega1(log D)) has basis omega_j=z^(j-1)dz/(z^4-1), j=1,2,3, phases i,-1,-i.

The **typed direct sum** C5 direct_sum F3 is eight-dimensional with Z4 degrees (0,0,1,1,2,2,3,3). The same geometry also gives a 3-to-5 exact sequence and a natural isomorphism (with orientation twist) from F3 to a subspace of C5. Therefore physical independence of the two summands does NOT follow from their being different mathematical types. It must be established on the actual source BRST/Dirac Hilbert space, not inferred from dimensions.

Conditional on independently realizing all eight as positive chiral complex CAR fermions, realification supplies 16 Majorana fields, 4 in each Z4 phase. The eight-complex-mode even Fock halfspinor has dimension 128 with phases (32,32,32,32). The 120 so16 bilinear currents (Lambda² of the 16 Majorana vector) have phases (28,32,28,32). Thus the standard E8=so16 plus positive halfspinor construction has full 248 current phase counts:

    (28,32,28,32)+(32,32,32,32)=(60,64,60,64).

This is EXACTLY the previously audited native E8 current grade-one character, no arbitrary 256-dimensional End(16) cover and no eight-state subtraction. It is a conditional CAR/VOA construction, not an original physical source derivation.

## Exact native E8 root and fixed Lie algebra comparison

In the standard E8 D8+ chart, use the geometric torus vector

    q_geom=(0,0,1,1,2,2,3,3).

For a root alpha, its Z4 phase is i^(alpha dot q_geom). Counted on all 240 native E8 roots, this gives 52,64,60,64 roots in grades 0,1,2,3, respectively. With eight invariant Cartan currents the result is 60,64,60,64.

The 52 fixed roots split exactly into root systems D5 (40 roots, rank 5) and A3 (12 roots, rank 3), giving the native 60-dimensional fixed Lie algebra D5 plus A3.

## An explicit 31-reflection E8 Weyl equivalence

The existing repository native glue degree is sum(root_doubled[5:]) modulo 4. This is the torus character of

    q_native=(0,0,0,0,0,2,2,2).

Fix eight E8 simple roots by their *doubled* coordinates, in the order:

    ( 1,-1,-1,-1,-1,-1,-1, 1)
    ( 2, 2, 0, 0, 0, 0, 0, 0)
    (-2, 2, 0, 0, 0, 0, 0, 0)
    ( 0,-2, 2, 0, 0, 0, 0, 0)
    ( 0, 0,-2, 2, 0, 0, 0, 0)
    ( 0, 0, 0,-2, 2, 0, 0, 0)
    ( 0, 0, 0, 0,-2, 2, 0, 0)
    ( 0, 0, 0, 0, 0,-2, 2, 0).

The Gram/Cartan determinant is 1. Apply the following one-based reflection indices in order:

    4,2,3,1,4,5,4,2,3,1,4,3,5,6,5,4,
    2,3,1,7,6,5,4,2,3,1,4,5,6,8,7.

Direct exact integer arithmetic gives

    w(q_geom)=(0,0,0,0,4,-2,2,2)
              =q_native - 4*alpha_7,

where alpha_7 denotes the *actual* undoubled seventh root listed above. All E8 lattice root pairings are integral, so for each of the 240 roots alpha,

    i^(alpha dot w(q_geom)) = i^(alpha dot q_native).

This is a constructive proof of Weyl conjugacy of the TWO ENTIRE order-four root characters, not just the same character dimension. Every root phase is compared in the independent regression. The standard E8 current bracket is respected by Weyl automorphisms at the abstract Lie algebra level.

## Crucial physical limits

1. This starts with D=mu4; the v512 nonsquare V4 counterwitness remains. It does not derive the square from an unmarked spacetime-free rule.
2. The 5- and 3-dimensional geometric spaces are related by an exact sequence. Treating them as independent physical CAR fields is a new source-level claim requiring the actual original positive Gram, ghost cancellations, real structure and zero modes. Eight complex vector-space dimensions do NOT prove 16 physical Majorana fields.
3. A positive chiral free-fermion VOA with 16 Majoranas and its integral-spin spinor simple-current extension can realize E8 level one; extension locality/GSO choice and physical original operator identification remain open.
4. The native R113 Q0, Z0, ZU, ZY, ZUY, original bosonic endpoint integration, the same absolute reference, and physical time remain uncomputed.
5. This supplies no interacting 3+1D parent, complete Standard Model chirality or universal spin-two gravity.

Next decisive source experiment: identify the actual independent one-particle CAR images of the five meromorphic sections AND three logarithmic differentials before and after original BRST/Hodge reduction, and check that their full positive Gram rank is 8 rather than 5,2 or otherwise. Then bind their physical quarter action and original charged Ward currents through the explicit Weyl word on the same time, source and reference. Avoid declaring algebraic Weyl conjugacy to be physical F1.
