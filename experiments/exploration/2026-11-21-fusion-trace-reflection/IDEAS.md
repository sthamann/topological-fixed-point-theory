# 21 November: PF/AF state, reflected Doob time and hard chirality gates

**Gebiete:** operator-algebraic state reconstruction, symbolic dynamical systems, Markov semigroups, modular fusion categories, anomaly theory. **Status:** stage 1 mathematical candidate, original TFPT F1 OPEN.

## 1 — Repair a false finite/infinite state mismatch

Input only the rank-two Fibonacci fusion rule. Fusion path adjacency A=[[0,1],[1,1]], PF phi=(1+sqrt5)/2, right r=(1,phi), initial vacuum. After n insertions the finite fusion AF algebra is M_{d_n(0)} ⊕ M_{d_n(tau)}, d_n=(F_{n-1},F_n). Define

    tau_n(X0⊕Xtau) = phi^(-n) [Tr(X0)+phi Tr(Xtau)].

This is a positive, normalized, **exactly compatible** trace under path refinement because A*r=phi*r. It is the natural PF trace on the corresponding stationary AF system (with chosen root). For n=5 the total-tau block is M5 and the conditional trace is exactly I5/5. Restrict to its first four fusion histories: the charges split 2+3 with conditional probabilities 2/5 and 3/5. This **equals the finite full-B5 braid-invariant trace**, whereas unconditional n=4 PF probabilities are 2/phi^4 and 3phi/phi^4. Thus the previous mismatch compared different conditioning and is **not a general no-go**. Important: imposing the total tau charge is a physical boundary choice, not derived.

This classical/AF finite-to-infinite trace bridge is generic for ALL primitive adjacency graphs: baseline **3/3** labeled primitive binary matrices, so no independent specificity for Fibonacci.

## 2 — Exactly reflected stochastic time, derived from the same graph

PF Doob matrix P=phi^-1 diag(r)^-1 A diag(r) is

    P=[[0,1],[phi^-2,phi^-1]].

It is row-Markov, reversible with pi=(1,phi²)/(1+phi²), spectrum {1,-phi^-2}. The one-step transfer has a negative eigenvalue and is not a symmetric positive Euclidean e^-H. The *paired reflected* update

    T=P^dagger P=P²=Pi_pi+phi^-4(I-Pi_pi)

is a strictly positive Markov contraction, spectrum {1,phi^-4}. Therefore the **unique positive selfadjoint** generator with a paired step of unit time is

    H = 4 log(phi) (I-Pi_pi).

The reflected kernel K_{mn}=pi0²+pi0*pi1*phi^-4(m+n) is positive semidefinite for all nonnegative integer m,n, all finite complex coefficient combinations. This is a valid *0+1 dimensional* reversible Markov semigroup / two-point reflection witness, NOT the full OS fermionic/gauge source or spacetime energy. Reflection time units are not physically calibrated.

Within **2-state 0/1 primitive adjacency** rules: 3 of 16 labeled graphs are primitive, exactly two (same unlabeled Fibonacci type) are invertible and have a bounded finite-time reflected e^-H; the full binary shift is rank-deficient. The binary alphabet and paired reflection are upstream choices.

## 3 — The two natural times still differ (the decisive bridge gate)

The gauge action alpha_t on the Cuntz-Krieger extension acts by word-length grading. Every fixed-length AF block matrix unit E_{wv} has degree zero, hence alpha_t(E_{wv})=E_{wv}. The uniquely selected gauge-KMS modular flow is **trivial on the finite braid M5 observables**. In contrast, the positive Doob H moves the two-state stochastic observables. They are not yet restrictions of a SINGLE original physical Hamiltonian. This is precisely where a source-selected local net/charged operator and positive-energy noncompact time are required.

## 4 — Restricted NO GO for bare Fibonacci -> E8/ordinary 3+1D fermions

Pure unitary Fib has R eigenphases exp(-4pi i/5), exp(3pi i/5): their squares are NOT identity, unlike permutation exchanges of ordinary 3+1D point particles. Its twist theta_tau=exp(4pi i/5), chiral charge c=14/5 mod 8. **Booker-Davydov (Journal of Algebra 355, 2012)** rigorously prove that Fib and all same-handedness tensor powers are completely anisotropic (no nontrivial condensable commutative separable ribbon algebras). Simple condensed Fib^n cannot produce trivial-anyonic E8_1 by itself. Five tau across five copies have composite twist 1 and 20 same-chiral copies have c=56=0 mod8; the no-condensation theorem remains. A mirror Fib⊗Fib^rev admits a known diagonal Lagrangian algebra and can become topologically trivial, but E8 invertible chiral sector then requires additional dynamics/choice. The separate conditional SU5_1²/Z5->E8_1 algebraic path remains mathematically valid but not physically sourced from Fib.

## 5 — Global minimum-entropy rule is NOT unique

For each n>=2, the 0/1 companion adjacency with superdiagonal ones and final-row first-two ones is primitive, unimodular, and has char polynomial x^n-x-1. Its Perron eigenvalue falls phi (n2), plastic constant (n3), then tends to 1 as n->infinity. Therefore minimizing positive entropy **across alphabet sizes** does not select Fibonacci. A minimum alphabet axiom is extra.

## 6 — Two further sharp typed obstructions

The real adjacency/PF process admits conjugate + and - Fibonacci braidings, so by itself it cannot select an oriented primitive P1 winding sign. Every finite algebraic F/R calculation yields algebraic values, whereas c3=1/(8pi) is transcendental; the coefficient needs a legitimately derived analytic/angular/continuum operation, not a finite fusion identity. A branching history tree has exponential volume growth, not the polynomial growth of ordinary spatial R3; a nontrivial localization/quotient is required.

## Research verdict

**Stage 1, mathematical candidate** for common finite/infinite state and reflected positive *one-dimensional* time, **not** Stufe 2, 3 or 4. Original P1/P2 physics, absolute source, R113, three independent SM families, 3+1 Lorentz continuum, gravitation all OPEN. Previous state mismatch is *repaired*, while bare Fibonacci condensation/identification into physical E8 is conditionally **VERWORFEN**.

## Reproduction

Run `python experiments/exploration/2026-11-21-fusion-trace-reflection/check_new.py`. The full standalone four-script test and executed run log are in the user-delivered ZIP artifact associated with this session.

## Independent sources

- Booker and Davydov, *Commutative algebras in Fibonacci categories*, J Algebra 355 (2012), https://doi.org/10.1016/j.jalgebra.2011.12.029 .
- Kong, *Anyon condensation and tensor categories*, Nuclear Physics B 886 (2014), https://doi.org/10.1016/j.nuclphysb.2014.07.003 .
- PF/KMS graph states: https://doi.org/10.1016/j.jfa.2013.09.016 .
- Doob/PF transforms: https://par.nsf.gov/servlets/purl/10411837 .


## 7 — General Takesaki gate and explicit one-qubit witness

Let N=M5(C) be a finite von Neumann subalgebra of a full von Neumann algebra M and let omega be a faithful normal state on M with omega|N=Tr/5. By Takesaki's invariant-subalgebra theorem, if sigma_t^omega(N)=N for all t, the restricted flow equals sigma_t^{omega|N}=id. Therefore a *nontrivial modular-time action on visible M5* is possible only if M5 is not globally invariant and the time flow accesses larger source/memory operators. This is a theorem specifically about modular evolution, **not** arbitrary Hamiltonian motion.

Exact minimal memory witness: on M5 tensor M2 let Z5=diag(1,-1,0,0,0), Z2=diag(1,-1), eps=1/2, and rho=(I10+eps Z5 tensor Z2)/10. It is positive and faithful, with reduced rho5=I5/5, but

    [log rho, E01 tensor I2] = log(3) E01 tensor Z2 !=0.

Thus a SINGLE hidden qubit permits nontrivial global modular flow despite tracial visible state; zero hidden qubits do not. One-parameter eps is freely chosen and this is **not** an explanation of the original TFPT 63 memory directions or a source-selection theorem. Reproduce: modular_ancilla_witness.py. The next actual physical gate is to derive eps and the entire nontrivial memory cross-sector dynamics from an original P1/OS source, not to insert them as new premises.

Takesaki source: https://doi.org/10.1016/0022-1236(72)90004-3 .
