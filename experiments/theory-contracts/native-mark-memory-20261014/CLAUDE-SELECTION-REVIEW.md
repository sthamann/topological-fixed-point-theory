# Audit: weaker mixed-spinor selection and the finite clock readout

## Verdict

The proposed mixed-root criterion is mathematically stronger than the synthesis section's existing *even unimodular completion* hypothesis, but its source-level premise is still an added condition: the original microscopic seam must produce a nonzero local weight-one field in a mixed D-spinor/A-weight coset. Under that premise, the arithmetic really does force (d=4), and the resulting integral even closure is already unimodular E8. This is a useful sharpening of the conditional selector, not evidence that the seam produces the field.

The proposed obstruction to an E8-current-to-(M_{16}) bridge is also valid as stated **if** the map is a grade-preserving Lie homomorphism and source/target use the same ℚ4 clock. It does not rule out the synthesis's nonmultiplicative positive expectation, a selected CP response, or an infinite-dimensional extension. The finite ℚ4 average is not a continuous torus average.

## Mixed-spinor norm test

In `tfpt_explorer/docs/session-2026-09-28/TFPT_Gesamtsynthese_2026-09-28.tex`, §globaloriginselection, Satz 1 first derives even unimodular completions only for (d=4m^2); Satz 3 then assumes such a completion and a nonzero D-spinor-sector field of weight one, using (h_D=(d+1)/8) to select (m=1). Satz 3's proof is correct within that hypothesis, but it does not cover a field whose charge is already in a mixed discriminant coset without assuming the full completion in advance.

For the proposed weaker test, take the root lattice (R_d=A_{d-1}\oplus D_{d+1}). In the mixed coset consisting of a D-spinor class and the A-class (k\omega_1) (or its conjugate), the minimum norm is

\[
\mu(d,k)=\frac{d+1}{4}+\frac{k(d-k)}{d}.
\]

This follows from the spinor-coset minimum for (D_{d+1}) and the fundamental-weight norm in (A_{d-1}). For any weight (v\in R_d^*\) and root-lattice vector (r\in R_d),

\[
\|v+r\|^2-\|v\|^2=2(v,r)+\|r\|^2\in2\mathbb Z,
\]

since (v\) pairs integrally with (R_d) and (R_d) is even. Therefore a norm-two vector in this coset requires its minimum to be 2: the minimum is positive and at most 2, and all other norms in that coset differ by an even integer. Setting \(\mu(d,k)=2\) gives

\[
4k(d-k)=d(7-d).
\]

For integer (d\ge2) and (0\le k<d), the only solutions are ((d,k)=(4,1),(4,3),(7,0)). The last is the **pure** D-spinor coset, not a mixed A/D charge; it cannot satisfy the proposed mixed-field condition. The two mixed solutions are exchanged by A-side conjugation.

At (d=4), (R_4=A_3\oplus D_5) has discriminant group \(\mathbb Z_4\oplus\mathbb Z_4\). The mixed generator ((s,\omega_1)) has order four and discriminant norm (5/4+3/4=2\equiv0\pmod{2\mathbb Z}). Its generated overlattice has index 4, hence determinant (16/4^2=1); it is positive, integral and even. Thus it is the E8 lattice (the unique positive even unimodular rank-eight lattice). No strict integral overlattice of a unimodular integral lattice exists, since its determinant would be (1/n^2\notin\mathbb Z) for index (n>1). This shows the mixed root itself completes the lattice in the (d=4) case.

At (d=7,k=0), adjoining the pure D8 spinor gives (A_6\oplus E_8), determinant 7 and rank 14. It is not mixed and is not unimodular as a whole. This is a useful counterexample to any unqualified claim that a spinor root alone always selects E8.

The exact inference is therefore:

1. **If** the local source supplies a nonzero weight-one field in a mixed D-spinor/A-coset, its charge must be a norm-two root.
2. The norm congruence and equation above force (d=4).
3. That mixed root generates the order-four isotropic glue and E8 closure.
4. Neither the existence of that local field nor its identification with the original TFPT seam follows from the lattice arithmetic.

The synthesis already states the relevant physical distinction near its current-selection discussion: the GIT lattice selector is not itself an analytic/local generator; `SEAM.SIMPLECURRENT.GENERATOR.01` remains open for the net-level identification. Atlas dossiers also retain `E8.GLU.01` at lattice scope and distinguish later formal glue classification from physical selection. So call the arithmetic a sharpened *conditional origin test*, not a newly completed origin derivation.

## Root and current counts; what the neutral average removes

For the E8 branching already printed in the synthesis around its (\operatorname{Spin}(10)\times SU(4)) decomposition,

\[
248=(45,1)\oplus(1,15)\oplus(10,6)\oplus(16,4)\oplus(\overline{16},\overline4),
\]

the root-grade counts for the stated ℚ4 mark are 52, 64, 60, 64, totaling 240 roots. The neutral grade also has the eight Cartan currents, so the weight-one Lie algebra grade dimensions are 60, 64, 60, 64, totaling 248. The distinction between root counts and full current dimensions is necessary; adding 8 only to grade zero is correct.

The finite (M_{16}) clock in the synthesis §“kleinste entscheidende Test” is explicitly

\[
u=\operatorname{diag}(1,i,-1,-i,1,i,-1,-i,\ldots),\qquad
E_0(X)=\tfrac14\sum_{r=0}^3u^rXu^{-r}.
\]

The source function `verification/v993_minimal_defect_selector.py::clock` implements `diag(i**(j % 4))`; `group_E` averages its four conjugations. This is a discrete \(\mathbb Z_4\) clock expectation, not Haar averaging over a continuous torus. It kills every nonzero clock grade (1, 2 **and** 3), not only grades 1 and 3. For the specific witness (X=|1\rangle\langle0|), (X) has grade 1 and (X^*) grade 3, while (X^*X=|0\rangle\langle0|) has grade zero; consequently (E_0(X)=E_0(X^*)=0) and (E_0(X^*X)=X^*X\ne0). The composite-response statement is correct, and demonstrates exactly why (E_0) is not multiplicative.

There are thus two different and both important facts: (E_0) erases individual charged operators; it preserves selected neutral products when the full word is formed before averaging. The counterexample does not show that all full E8-current information survives in (M_{16}), and it does not show that all source responses are lost. It is a finite-algebra test, not a physical-time or spatial interpretation.

## Conditional no-go for a grade-preserving Lie map

Assume the source Lie algebra is the full simple complex (\mathfrak e_8) of weight-one zero modes with the indicated (\mathbb Z_4) grading, the target is (M_{16}(\mathbb C)) with this same clock, and (\varphi: \mathfrak e_8\to\mathfrak{gl}_{16}(\mathbb C)) is a grade-preserving Lie homomorphism. The fixed Lie subalgebra of target clock conjugation is

\[
\mathfrak{gl}_{16}^{(0)}=\mathfrak{gl}_4\oplus\mathfrak{gl}_4\oplus\mathfrak{gl}_4\oplus\mathfrak{gl}_4,
\]

because each of the four clock eigenspaces has multiplicity four. The D5 subalgebra lies in source grade zero. Each projection of (\varphi|_{\mathfrak{so}_{10}}\) into one (\mathfrak{gl}_4) is either zero or injective, since (\mathfrak{so}_{10}\) is simple. Injectivity is impossible by dimensions (45>16), so all four projections vanish. Thus the kernel of (\varphi) contains the nonzero D5 algebra. Since (\mathfrak e_8\) is simple, that kernel is all of (\mathfrak e_8); the homomorphism is zero.

This is a genuine finite-dimensional obstruction to a **multiplicative, grading-respecting Lie-algebra representation** of the whole current algebra in this clock space. It is not a no-go for a CP map, a conditional expectation, a nonmultiplicative current response, or an infinite-dimensional extension. In particular, the synthesis's witness (X^*X) warns against replacing the complete operator word by a product of separately averaged factors. So the correct interpretation is “the finite neutral image is not a faithful graded Lie realization,” not “the finite readout cannot retain any joint response.”

## Scope and source boundary

The root/VOA arithmetic assumes the chosen family (A_{d-1}\oplus D_{d+1}), its integral norms, and that the candidate field charge belongs to the stated mixed discriminant coset. The “local” designation and the map from the original seam to that vertex operator remain additional physical/operator-identification premises. The no-go assumes a grade-preserving Lie homomorphism and exactly the finite (u) above. None of these results picks the architecture from P1/P2 or establishes physical time.
