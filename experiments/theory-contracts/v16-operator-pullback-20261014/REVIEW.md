# Independent review

14 October 2026. Reviewed `PROOF.md`, `FOLD-AUDIT.md`, and the two local exact checkers. Both reran successfully: the packet/vector inverse reported 32,010 checks, including all 75 supplied-connection basis/lift cases and all 60 mixed directions; the fold audit reported 702 checks.

The operator inverse is correctly scoped. With the stated Clifford convention, the chiral trace Gram is −32 on matching generators and zero off-diagonal, so the formula recovers the unique `spin(16)` Lie-algebra coefficients whenever the supplied packet connection lies in that image. The explicit residual check retains components outside the image. The original marked `C_0,X_(0,1)` witness gives the stated two real vector generators in the chosen frame. This transports the supplied connection/action; it does not construct the raw Clifford multiplicity map `Xi`, choose the lift or physical connection, or identify the packet Hilbert space with the raw-field Hilbert space. The central-sign caveat is necessary and correctly stated.

The fold is a genuine finite-cover pushforward: the eight components are a unitary decomposition of the original one-particle space, the cyclic endpoint carries the original phase, and the folded covariance/state is `F P F*`, not a product of eight seas. Localizing a quotient interval means taking the union of eight separated cover arcs. Under that changed localization and the quotient conformal action, the transported short-distance stress contributions give central charge 8; the original point-local cover net remains central charge 1. This corrects the overly broad objection that a fold can never yield a `c=8` quotient net. The source still has to select this quotient localization physically.

The periodic (E_8) conjugacy is exact: the matrix is orthogonal with determinant (+1), maps all 240 roots bijectively, and the computed shift is in (E_8). The phase equality for every root establishes the stated periodic root action, not merely its centralizer type. The quarter phases are separately retained: their `g^2` root action has order 8 in this stipulated E8 lattice model, so they cannot be silently replaced by the periodic case. The folded twisted-sector weight formula is correctly presented as a conformal-sector weight, not as a physical excitation energy or field dimension without the relevant OPE.

The remaining physical gates are explicit rather than hidden: P1/P2 have not selected the Möbius quotient/locality or Pin/gauge lift; the E8 extension has not been derived from the raw source; the root automorphism does not preserve the elementary integer fermion lattice; and the fold plus packet inverse still does not supply `Xi`, its state/time transport, or the raw mixed-current Ward identity. No claim of a complete TFPT/T1–T8 solution follows.

## Addendum: source U(1) transport in the fold

The new §4a identity is correct in the fold's stated phase convention. With (b=(1,\ldots,1)), direct multiplication gives (Mb=(0,0,2,0,0,2,0,0)=2(e_3+e_6)). Since \(\nu_\alpha=\nu_0+(\alpha/4)b\), it follows exactly that

\[
M\nu_\alpha=t+\lambda+(\alpha/2)(e_3+e_6),\qquad \lambda=M\nu_0-t\in E_8.
\]

For each source root \(r\) and target root \(p=Mr\), the added phase is therefore \(\exp(\pi i\alpha(e_3+e_6)\cdot p)\). The rerun checked this identity over all 240 roots for all four \(\alpha\) values; `fold_check.py` passes 1,668 exact checks. Thus \(Q'=2(H_3+H_6)\) is the transported Cartan charge on the stipulated E8 root currents. This is a charge/action dictionary inside the already assumed folded E8 extension, not a field-level map for the original CAR fermions or proof that P1/P2 select a compensating gauge action.

One convention should remain explicit when comparing to the source parameter: source modes \(e^{i(j-a)\theta}\) have full-circle phase \(e^{-2\pi i a}\), whereas the fold boundary is written with \(e^{+2\pi i\alpha}\). Hence \(\alpha\equiv-a\pmod1\) unless the orientation is reversed. Both quarter phases are covered by the tests, but the signed identification must be fixed before transporting the physical generator \(H_a=L_0+(a-1/2)Q\).

## Addendum: packet realization of the charge-selected lift

The new §4b case is consistent with the preceding charge match. For (Q_E=2\operatorname{diag}(0,0,1,0,0)), the original packet formula with (m=(0,-1,-1,0)) yields a diagonal Hermitian charge on the 64-complex packet; after the stated realification and trace inverse, the vector connection is exactly
\[
\Omega_V=-2L_{4,5}-2L_{10,11}=iQ'
\]
in the interleaved real-vector convention. There is no residual outside \(\mathfrak{spin}(16)\). This is the actual packet formula evaluated at the matched lift, not a second fitted vector action. The fixed `e_3/e_6` and family-plane orientation are part of that dictionary; simultaneous reversal of carrier and family charges preserves the (a_i), while changing only relative family handedness changes the lift choice.

The checker reruns at 32,133 exact checks. The existing 75 connection basis/lift cases and 60 mixed actions remain, with the new source-charge packet case added. This confirms the unique lift within the assumed packet/E8 dictionary; it still does not derive that dictionary or the raw field map from P1/P2.
