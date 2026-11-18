# Two-face collar, actual face currents, and relative endpoint transport

Research ID: `UR.COLLAR.FACE_CURRENT.01`. Date: 2026-10-21.
Verdict: **conditional exact source bridge; physical P1 intervention not selected**.
This combines existing TFPT data, not new general current-algebra mathematics.

## 1. Original data and explicit premises

Use the original selfadjoint circle operator `B_alpha=-i d_theta+alpha`,
`alpha=1/4` or `3/4`, its actual periodic/twisted domain, and the original
positive normal length `ell`. The original product equation gives the closed
transfer `T=exp(-ell B)` and its orthogonal two-face graph projection

`G=[[A,C],[C,1-A]]`, `A=(1+exp(-2ell B))^-1`, `C=1/(2 cosh(ell B))`.

For two-sided unbounded B, T is not declared bounded. G is bounded and is
the projection onto its closed graph. The circle Fourier blocks `G_r` are
rank-one 2x2 projections with `tr_faces(G_r)=1`.

As in the original collar proof, ADD the charged CAR algebra on both faces
and the pure gauge-invariant quasifree state with covariance G. This field
and statistical premise is NOT derived from P1. Its second-face occupation
is `n_r=(1+exp(beta(r+alpha)))^-1`, `beta=2ell`.
The original one-face negative sea is `P_-=1_(B<0)`.
Do not identify a geometric graph projection, a mixed occupation, a finite
native code projector or the physical admissibility projector with each other.

Original sources reused:
* `seam-full-collar-reconstruction-20261015/PROOF.txt:63-143`, especially its
  explicit CAR premise and selected labeled second face at lines96-116.
* `joint-native-algebra-selection-20261016/ORIGINAL-ORDERED-CURRENT-KERNEL-AND-SEAM-GRADING.txt:33-161`.
* `joint-native-algebra-selection-20261016/ORIGINAL-DENSITY-GEOMETRY-AND-TIME-HISTORIES.txt:36-99`.
* `joint-native-algebra-selection-20261016/ORIGINAL-P1-VARIATION-GAP-AND-EVENT-SELECTION.txt:37-134`.
All these paths are relative to `experiments/theory-contracts/`.

## 2. The same graph fixes the one-face current cocycle

Let `P0,P1` be the two labeled FACE projections, not native 3+2 marks.
For a smooth real circle function f, let `A_j(f)=Pj M_f` on the doubled
one-particle space. In the G-Fock representation, normal ordered quadratic
generators exist: `[G,A_j(f)]` is Hilbert-Schmidt. For fixed Fourier shifts
this follows because `G_r-G_(r+k)` decays at either tail and `C_r` decays;
smooth Fourier coefficients allow the convergent sum. Constant face terms
are also allowed because C is Hilbert-Schmidt on this gapped circle.

Their central commutator is obtained by the Wick contraction

`c_G(A,B)=Tr[G A(1-G)B G - G B(1-G)A G]`.

For second-face multiplication, the summand is
`n_r(1-n_s)(f_(r-s) g_(s-r)-g_(r-s) f_(s-r))`.
Combining the k and -k terms and using the absolutely convergent telescoping
sum `sum_r(n_r-n_(r+k))=k`, k>0, gives EXACTLY

`c_1(f,g)=sum_(k>0) k (conj(f_k)g_k-f_k conj(g_k))`.

It is independent of ell and alpha. The opposite face has `c_0=-c_1`;
mixed-face central commutators vanish. The latter follows from the symmetric
coefficient `C_r C_s`, whose antisymmetric Fourier sum cancels.
Consequently common multiplication on BOTH faces has `c_common=0`.
This also follows directly from `tr_faces G_r=1` in the convergent trace.
For nonzero Fourier modes G-normal-ordering and the original sea normal
ordering agree. The zero mode in section3 uses the ORIGINAL relative Q;
its finite offset from G-normal-ordering is `Tr(n_ell-P_-)` and is retained.

Thus the original second face has precisely the original current relation
`[J_m,J_-k]=m delta_(mk) I`, not a newly fitted current level. Finite windows
in the checker are controls of this proof, not finite CCR representations.
Their boundary commutators are retained; no infinite identity is inferred
from a finite trace of a commutator.

For the source profiles `f=cos(theta)`, `g=sin(theta)`, `c_1=-i/2`, while
the common two-face cocycle is zero. More stringently, choose
`f=sin(4theta)`, `g=1-cos(4theta)`. BOTH vanish at all four original marks
`theta=j pi/2` and are quarter-periodic. Their cocycle is `c_1=-2i`.
With `W(f)=exp(iJ(f))`, the current rectangle is
`W(sf)W(tg)W(-sf)W(-tg)=exp(2ist)I` on the second face,
whereas its common two-face counterpart is I. This is a conditional
coherent-history diagnostic, not an executed detector or probability law.

## 3. Complete current state and time, with all charge carry retained

Wick contraction in the SAME G-state gives, for m>0,

`<J_m J_-m>=sum_r n_r(1-n_(r+m))=m/(1-exp(-beta m))`,
`<J_-m J_m>=m/(exp(beta m)-1)`.

Each term of the second sum is exp(-beta m) times its counterpart in the
first; their difference telescopes to m. No spectral cutoff is used.
In this circle source the particle/hole product is convergent, so the
restricted thermal state is normal in the original negative-sea Fock
representation. This extra circle-specific fact is not asserted for an
arbitrary collar or arbitrary Fock representation.

Reuse the original FULL charge-sector/Hall decomposition. The unscaled
relative excitation energy is `H0=8K_alpha` with
`E_alpha(q)=q(q-1)/2+alpha q` and oscillator energies m. Thus

`Z_rel(beta)=sum_(q in Z) exp(-beta E_alpha(q)) * prod_(m>=1)(1-exp(-beta m))^-1`.

This equals the original full particle/hole product. No charge is replaced
by q mod4; all winding carry remains. The formula is relatively normal
ordered: a given absolute sea/eta contribution is retained separately,
not reset or declared zero.
For `f_0=average(f)`, the whole current characteristic function is

`omega_ell(W(f))=chi_ell(f_0) exp[-1/2 sum_(m>0) m coth(ell m)|f_m|^2]`,
`chi_ell(t)=sum_q exp(-2ell E_alpha(q)+itq)/sum_q exp(-2ell E_alpha(q))`.

The nonzero-mode expression follows from the thermal oscillator trace in
that original, complete boson-fermion/Hall decomposition. Smooth profiles
retain EVERY Fourier mode. Existing source time acts as
`U(t)J_-m U(t)^*=exp(-itm/8)J_-m`; geometric ell is not physical time.
After expressing Heisenberg profiles in this convention, every finite
ordered Weyl correlation follows by the original BCH cocycle and this
characteristic function. If a Schrodinger history retains a final U(T),
it must also be retained: in a thermal trace its relative factor is obtained
by `beta -> beta+iT/8`, including the ratio `Z_rel(beta+iT/8)/Z_rel(beta)`.
The prescribed absolute time phase remains separate. Removing it without
normalization would lose information.

Moreover `Tr|n_ell-P_-| <= [exp(-2ell alpha)+exp(-2ell(1-alpha))]/[1-exp(-2ell)]`.
Hence the long-collar limit gives the original negative sea in trace norm,
`chi_ell ->1` and `coth(ell m)->1`, recovering the existing vacuum Weyl kernel.
This selects the vacuum only within this GIVEN operator/CAR/long-collar
class. It neither selects the original source nor makes sign-only data
sufficient to reconstruct time magnitudes.

## 4. Two different full lifts with the same face marginal

For `g_phi=M_exp(i phi)`, common source conjugation gives
`G_common=(g_phi direct_sum g_phi) G (g_phi^* direct_sum g_phi^*)`.
A second-face action gives instead
`G_face=(I direct_sum g_phi) G (I direct_sum g_phi^*)`.
Their second-face occupations both equal `g_phi n g_phi^*`, but their
off-diagonal blocks are respectively `g_phi C g_phi^*` and `C g_phi^*`.
They therefore differ on the full source. Section2 also distinguishes their
coherent phase laws. Equality of a reduced covariance cannot identify them.
This does not say that the two physical contexts must both be available.
The distinction needs access to the full source or coherent control of the
lifts. If only one-face conjugation channels are available and the other face
never recouples, those channels agree; a central phase alone is not an
ordinary unmonitored-channel discriminator.

## 5. A concrete geometric relative endpoint family

The second-face lift is the graph of `T_phi=g_phi T`, not in general the
graph of `g_phi T g_phi^*`. This holds for the actual closed unbounded T:
postcomposing by bounded unitary g preserves closedness and its domain,
and `Graph(gT)=(I direct_sum g)Graph(T)`.

A geometric construction exists without choosing a new bulk interaction.
Choose any smooth interpolation chi(r), flat near the two endpoints,
with chi(0)=0, chi(ell)=1, and let `g_chi=exp(i chi(r)phi(theta))`.
For the original scalar TRANSPORT BLOCK `L=partial_r+B_alpha`,

`L_phi=g_chi L g_chi^*=partial_r+B_alpha-i chi_prime phi-chi phi_prime`.

Its Cauchy relation is exactly `T_phi=g_phi T`. The displayed source jet
is linear in phi; its mixed second derivative is zero. For L*L or the
spectral action the corresponding quadratic/contact terms are NOT zero
and must be differentiated, not discarded. In a full Dirac real/graded
module use the original lifted bundle unitary and actual normal symbol,
not an identification of this scalar block with every Dirac component.

Any two chi with the same endpoints differ by a bundle unitary equal to I
on BOTH faces. They therefore give the same boundary graph. The interpolation
is an internal frame extension, not an additional measurable coupling/profile.
Operator domains, real structure, reflection, reference and admissible
subspace must be transported by their actual bundle maps. With the ENTIRE
datum transported this can be a mere frame change, not a physical event.

This closes the MATHEMATICAL relative-endpoint construction in the declared
source class. It does not prove that P1 selects its physical execution,
that its unitary belongs to the ORIGINAL fixed admissible local net, or
that a fixed reference/domain can be left untransported without new terms.

## 6. Consequence for the actual TFPT target

The shared collar, field premise and labeled face determine state, current
level and its ordered phase consistently. Independent current normalization
or an arbitrary oscillator is unnecessary in this declared class.
The naive common-conjugation identification with a one-face current is false;
relative endpoint transport supplies a precise conditional repair.
Physical source selection, phase/reference interpretation, the original
P1-to-native/E8 field identification, event rates, pointer/record dynamics,
3+1D chirality, couplings and gravity remain unproved. No T1-T8 goal closes.
The already given native4/5/20/carry, M/O/NS and full eight-preimage field
support are not reidentified or replaced by this two-face calculation.

The primitive scalar `Delta eta` evaluated on a covariantly isospectral
family is not automatically a family connection with this curvature.
The required line-with-connection compatibility is explicitly still open in
`documents/tfpt_research_contracts/sections/017-source.tex:239-285`,
`SEAM.DETLINE.UNIFICATION.01`. No 1D calculation proves its 4D assertion.
Standard implementer-cocycle background: [Langmann, JMP35 (1994)](https://arxiv.org/abs/hep-th/9305014).
