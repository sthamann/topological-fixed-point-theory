# F1 P1: the physical seam requires a projectively consistent SOURCE, not pure GHZ at every cutoff

**16 November 2026. Exact independent mathematical continuation of the joint
boundary-matching and cyclic-constraint candidate. Physical F1 remains OPEN.**
No native R113 all-field integrals or physical Q0, P_adm, OS time, spin
structure, absolute Z0, 3+1D local QFT or spin-2 have been calculated.

## I. NEW EXACT SOURCE NO-GO: the previously selected unique pure P1 seam is NOT a projective family

Let d>=2, n>=2, and let E_n:C^d -> (C^d)^tensor n be the
basis-repetition isometry E_n|j>=|j>^tensor n. Put

    |GHZ_(d,n)> = (1/sqrt(d)) sum_j |j>^tensor n,
    rho_n = |GHZ_(d,n)><GHZ_(d,n)|,
    mu_n = (1/d) sum_j |j>^tensor n <j|^tensor n.

For d=4 the original two/three region combination
"edge matching + global Z4 cycle constraint" selects rho_2 and rho_3
uniquely up to phase, separately at each n.

However ordinary physical discarding of the last region gives, for
EVERY n>=3 and d>=2,

    Tr_last(rho_n)=mu_(n-1) != rho_(n-1),
    Tr_last(mu_n)=mu_(n-1).

The spectral difference has exactly

    spec[rho_(n-1)-mu_(n-1)] =
        {(d-1)/d (x1), -1/d (x(d-1)), 0 (remaining)}.

Therefore trace distance D(rho,mu)=(d-1)/d and pure-state overlap
Tr(rho mu)=1/d. At d=4 the exact numbers are 3/4 and 1/4. This is
an O(1) violation of ordinary quantum marginal/projective
consistency, NOT a numerical tolerance issue.

A stronger general theorem: if a bipartite marginal rho_AB of any
state rho_ABC is PURE, then rho_ABC=rho_AB tensor rho_C. A pure
pre-existing seam therefore CANNOT acquire correlations with a new
physical region while keeping the same old marginal under partial
trace. The result holds beyond GHZ and beyond four marks.

The discrepancy is visible even to an operator invariant under the
S tensor S cycle used as a constraint. Set O_AB=S tensor S, with
S|j>=|j+1 mod 4>. It commutes with the cycle generator itself, but

    Tr(O_AB rho_2)=1,
    Tr(O_AB Tr_C rho_3)=0.

This is an algebraic witness in the original cyclic-invariant
joint-observable sector, not merely an inaccessible offdiagonal
matrix component. Physical admissibility of this observable must
still be checked against native P_adm.

The result does NOT exclude an original nontracial or source-weighted
RG map that changes the observable embedding. It excludes the
specific interpretation as one compatible family under ordinary
partial traces of added independently accessible tensor factors.

## I-b. Stronger uniqueness theorem: entire infinite tower is classical under these three axioms

Fix four primitive labels and an entire tower of states rho_n on
(C4)^tensor n for every n>=2. Impose THREE independent premises:

(i) support(rho_n) subset C_n=span{|j>^tensor n : 0<=j<4},
    i.e. exactly matched sector labels across all n regions;
(ii) (S^tensor n) rho_n (S^tensor n)^dagger=rho_n, with S
     a transitive cyclic permutation of the four labels;
(iii) Tr_last(rho_(n+1))=rho_n with ordinary density-matrix
      marginalization and a fixed local observable embedding.

THEN the UNIQUE normalized projective family is

    rho_n=mu_n=(1/4)sum_j |j>^tensor n <j|^tensor n.

Proof: every supported fine state is E_(n+1) M E_(n+1)^dagger,
where E_(n+1)|j>=|j>^tensor(n+1) and M is a positive 4x4
unit-trace coefficient matrix. Partial trace eliminates all
offdiagonal entries M_jk, since <k|j>=delta_jk, yielding
E_n diag(M_jj) E_n^dagger. Because EACH rho_n has a fine
extension (condition iii), EACH rho_n is diagonal in the
matched label basis. Under condition ii cyclic covariance
forces its four weights to coincide. Normalization fixes 1/4.
Conversely this family obeys all three axioms exactly.

This excludes not only GHZ_n but ANY coherent offdiagonal
matched-label source, at ANY finite level of an indefinitely
extendable tower. The symbolic arbitrary 4x4 fine coefficient
matrix has all 12 offdiagonal entries eliminated, while the
cyclic weight constraint has rank 3 plus normalization rank 1.
This is a genuinely CONDitional source uniqueness proof,
not a physical no-go for arbitrary original TFPT weighted
projective maps or gauge-invariant nonfactorizing algebras.

Conclusion: demanding exact global matching, full cyclic symmetry
and ordinary quantum partial-trace projectivity gives exactly
the incoherent classical mixture, NOT a microscopic quantum
amplitude law. To retain physical phase information, the
actual original TFPT source MUST replace at least one premise:
physical matching can involve nontrivial edge multiplicities;
observable refinement can use a genuine quantum coding
isometry/nonlocal decoder; coarse gluing can be a conditioned
cap with absolute 1/4; or the source can store nonlocal
history memory in an enlarged algebra. All changes must be
derived from original R113 Q0 and P_adm, not guessed.

## II. THREE DISTINCT repairs, none equivalent by assumption

(A) Dephased ordinary source:

    mu_n=(1/4)sum_j |j...j><j...j|

is exactly projectively consistent under partial trace at every n.
It has four superselection/label sectors and no finite-local
n-partite GHZ coherence. On the full M_(4^n) algebra it is nonfaithful
(rank four, not 4^n). On its diagonal-code algebra it is tracial and
its Tomita modular Hamiltonian is scalar. Thus this is a valid
finite-volume consistency baseline, NOT a generator of nontrivial
TFPT physical translations.

(B) Conditioned physical cap:

    |+>=(1/2)sum_j|j>,    C=I tensor <+|.

The unnormalized output is EXACTLY

    C rho_n C^dagger = (1/4) rho_(n-1),
    p_success=1/4.

After k equivalent independent contractions the model success weight
is 4^(-k). Normalized postselection restores pure rho_(n-1), but
conditioning is not a deterministic linear CPTP coarse map.
The amplitude scaling is a potential ABSOLUTE Q0/Z0 source diagnostic,
NOT a physical predicted TFPT number. Any supposed "unit source
normalization" without the absolute scalar silently discards this
one-quarter factor.

The categorical deletion epsilon=sum_j<j| gives
(I tensor epsilon)E_2=I and recovers GHZ without a norm factor on
the code. Yet epsilon^dagger epsilon has top eigenvalue d=4, so
epsilon is NOT a legal probability-normalized measurement effect
on the unrestricted physical Hilbert space.

(C) Nonlocal isometric decoding:

    delta:C^4 -> C^4 tensor C^4, delta|j>=|jj>,
    W=I_4 tensor delta: C^16 -> C^64.

delta is an isometry, coassociative and cocommutative, i.e. the
copy map of a dagger special commutative Frobenius algebra in the
specified basis (Coecke, Pavlovic, Vicary 2008/2013;
https://arxiv.org/abs/0810.0812).
Here W rho_2 W^dagger=rho_3 and W^dagger rho_3 W=rho_2.

There even exists a deterministic CPTP decoder on the full M64:

    E(X)=W^dagger X W
      + Tr[(I64-W W^dagger)X] tau,
    tau=I16/16.

Both summands are CP; the output trace equals Tr X.
E(rho_3)=rho_2 exactly. This is not a partial trace over C:
it actively decodes the joint B,C code and CHANGES old AB
expectations. The operator W does not clone arbitrary unknown
quantum states, only a designated orthogonal pointer basis.

The original source must specify which coarse observable embedding
is physical. Merely exhibiting this mathematical CPTP decoder is
not evidence that the original R113 source uses it.

## III. Source-overlap law: unavoidable information/coherence tradeoff

Consider a pure source with an unobserved environment:

    |Psi_ABE>=(1/2) sum_j |jj> tensor |e_j>,
    ||e_j||=1.

Then

    rho_AB=(1/4)sum_(j,k)<e_k|e_j>|jj><kk|.

The offdiagonal PHYSICAL source amplitude is therefore exactly
the Gram of conditional environment states. Orthogonal records
<e_k|e_j>=delta_jk give mu_2; identical records give rho_2.

For equal real overlaps lambda in [0,1]:

    rho_AB(lambda)=lambda*rho_2+(1-lambda)*mu_2,
    fidelity with rho_2=(1+3lambda)/4,
    trace distance to rho_2=3(1-lambda)/4.

At lambda=1/2 the exact code eigenvalues are
{5/8 (x1), 1/8 (x3)} plus twelve zeros. Thus any claimed original
P1 information transfer must quantitatively determine the overlap
lambda rather than setting "pure" or "fully classical" by assertion.

The earlier six-history 6x6 Gram kernel K and the original R113
joint reflected field Gram are of this SAME SOURCE TYPE
(overlap matrices), not numerically identical matrices. No
tensor-product factorization is presumed: the actual source
must evaluate all history/field cross terms in one native Q0.

## IV. Quarter spin character: four distinct constrained states remain

The matching projector and the global cyclic projector produced
the m=0 state by demanding U|Psi>=+|Psi>.

But Z4 has FOUR characters. Replace the neutral projector by

    P_m=(1/4)sum_(k=0)^3 i^(mk) U^k,  m=0,1,2,3,
    U=S tensor S.

Every P_m has rank four and commutes with the matching rank-four
projector. Their intersections each have rank ONE:

    |Psi_m>=(1/2)sum_(j=0)^3 i^(mj)|jj>,
    U|Psi_m>=i^(-m)|Psi_m>.

All four obey equivalent mathematical projection and norm axioms.
Reflection j->-j exchanges m and -m. Thus a genuinely complex
quarter spin sector (m=1 or 3) is not selected without a physical
choice of orientation/charge/projective action. Its different
coherences are lost under ordinary tracing of an additional region.

This is a direct source-selection gate for original q0/P_adm
and the nine previous coherent permutation histories, not
a derivation of the native physical spin or E8 quarter root charges.

## V. Additional local-QFT exclusion: no clustering in raw label states

In the dephased projectively consistent mu_n let the local
four-sector label observable be
Z=diag(1,i,-1,-i). Its connected correlation at any separation
between different sites is

    <Z_a Z_b^dag>-<Z_a><Z_b^dag>=1.

For an infinite system where Z is a physical *local observable*
and a unique uniform spectral gap + short-range interactions
are demanded, this violates the standard exponential clustering
consequence (Hastings-Koma 2006,
https://arxiv.org/abs/math-ph/0507008).
Thus the simplistic mu_infinity mixture is not a standard unique
gapped clustering vacuum under those assumptions. The statement
does NOT exclude gauge/topological superselection where the raw
label Z is NOT a physical local observable, nor gapless theories.

The constraints themselves are global and require a nonlocal
all-sites clock operation to single out GHZ, so locality is not
yet recovered by this finite construction.

## VI. Exact physical F1 next gate: trace vs cap vs decoder must be source-selected

Require the native original R113 full all-field weighted operation
with unchanged P_adm, physical reflection, original q0/spin and
absolute Z0 on BOTH native blocking levels. Before promoting any
GHZ/five-point/quarter-clock source, calculate all:

1. Actual fine-to-coarse map on **the same fixed local observable
   algebra**, including whether the map is trace, postselected cap
   (and its absolute 1/4), or nonlocal CPTP decoding with modified
   observable identification. Do not substitute one for another.
2. Joint original source Gram G[(c,i),(d,j)] for the six clock
   histories and all Higgs/fermion/geometry source directions.
   Its environment overlaps and phase characters must come
   from the original field integral and physical P_adm, not an
   invented positive matrix.
3. Full Q0/Z0, ZU, ZY, ZUY, source-generated OS positive
   Dirichlet form/Higgs U7 kinetic and directed clock/stress
   responses at two genuine cutoffs with unbounded history memory
   tracked rather than erased.
4. Original physical algebraic net, faithful state and modular
   half-sided inclusion (or full OS reconstruction) to produce
   actual time and translation, then 3+1D chiral spectrum,
   native E8 current and propagating massless spin2.

A family of finite pure GHZ states is NOT enough for F1. A
categorical topological copying algebra is NOT an interacting
3+1D quantum field theory. A probabilistic cap's 1/4 is NOT
the source value until original physical Q0 is evaluated.

The underlying native source diagnostic remains
whole_source_bundle=null and spatial_tangential_circle_NOT_OS_time.
No ledger or F1 proof promotion is justified.

Reproduction:
research/originlab/f1_seam_projective_consistency.py
and test_f1_seam_projective_consistency.py.

External context for quantum multitime marginal consistency:
Milz et al., Quantum Stochastic Processes and Quantum Non-Markovian
Phenomena, PRX Quantum 2 030201 (2021),
https://doi.org/10.1103/PRXQuantum.2.030201.
