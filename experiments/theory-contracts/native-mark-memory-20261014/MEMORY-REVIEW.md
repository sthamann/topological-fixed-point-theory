# Exact thermal-GNS memory from the native Round-20 dynamics

**Independent review, 2026-10-14.** This review is restricted to the fixed
finite spatial lattice, the neutral total-charge sector, the already-defined
Round-20 Hamiltonian, and its faithful Gibbs state. It introduces no bath,
stochastic reset, new Hamiltonian, or new source law.

## Verdict

The proposed result is mathematically valid after two wording corrections:

1. The Hilbert space is the **completion of** `B(H)` in the faithful Gibbs
   norm, not `B(H)` itself with a new claim of Hilbert completeness. It is
   canonically isometric to Hilbert--Schmidt space by
   \(A\mapsto A\rho_\beta^{1/2}\).
2. \(C''(0)=-A^2-B^*B\) is an identity for the finite compressed matrix (or
   its matrix elements). It does not assert that every \(M_a\) belongs to
   \(\operatorname{Dom}\mathcal L^2\). Only
   \(M_a\in\operatorname{Dom}\mathcal L\) is needed.

With those qualifications, the original Hamiltonian and Gibbs state supply
an exact non-Markovian memory kernel

\[
K(t)=B^*e^{it\mathcal L_Q}B,
\]

and the native hopping makes \(K(0)=B^*B\) nonzero. This proves that the
four-dimensional local-grade compression is not the autonomous unitary
motion \(e^{itA}\). It does **not** make the GNS projection a measurement,
a completely positive conditional expectation, or a physically selected
forgetting operation. The local grade itself is inherited from the earlier
chosen E8 grading; the Round-20 source does not derive why this is the
physical readout.

There is also an exact constructive dictionary on each chosen native
\(e_8\)-transfer direction: the four grades form a constant \(M_4\) algebra,
while the fourth transfer is the genuine bilateral integer carry. The
resulting concrete orbit algebra is \(M_4(C(S^1))\), not a four-state
quotient. This is an internal re-expression of the original infinite charge
orbit. It is not a dynamical tensor factor of the full Hamiltonian and is not
canonically the \(M_4\) used in any earlier CP channel.

Extending that coordinate split over the full neutral charge lattice also
gives a decisive state statement: relative to this native grade versus
carry/kernel/scalar factorization, the original Hamiltonian is not a tensor
sum and its faithful Gibbs state is not a product. The conclusion concerns
this explicit source-derived factorization only, not every factorization of
the separable Hilbert space.

## Original-source ownership

The result uses the following registered sources.

* `source:75dfcc502b5bec528074`,
  `charged-lift-round16/PROOF.md`, lines 12--18 and 34--64: the physical E8
  lattice basis has \(b_8=s=(1/2)^8\); the integer functional
  \(j(x)=4\sum_{i=6}^8x_i-2\sum_{i=1}^5x_i\) is congruent modulo four to the
  older grade \(2\sum_{i=1}^5x_i\). It has
  \(j(b_1),\ldots,j(b_8)=(-4,-4,-4,-4,-4,2,2,1)\).
* `source:5959236bb77ca96b5633`,
  `local-charge-transport-round20/PROOF.md`, lines 12--14 and 51--77: a site
  coefficient vector \(n_x\in\mathbb Z^8\) denotes the physical lattice
  vector \(B_{E8}n_x\), and the exact on-site shifts and local spectral
  projectors generate the full site algebra. Lines 94--122 define the
  neutral bond hopping with every coefficient-basis generator
  \(p=e_1,\ldots,e_8\), prove its bounded positivity, and exhibit its
  nonzero charge-transfer matrix element. Lines 183--187 prove that total
  charge is conserved while individual site labels are not.
* `source:e39426b76acca89849ef`,
  `hopping-source-domain-round20/PROOF.md`, lines 15--60: the charge-diagonal
  terms and scalar coupling are combined with the bounded hopping
  \(L_J\); the full Hamiltonian is self-adjoint on the diagonal domain, and
  only total charge sectors reduce it.
* `source:18d340aeed16115daafa`,
  `scalar-charge-energy-round20/PROOF.md`, lines 38--67: the scalar/charge
  Hamiltonian has charge-diagonal \(D+V\) and the same off-diagonal
  \(L_J\). This confirms that the local grade projectors commute with every
  unbounded term and fail to commute only with hopping.
* `source:0d7cc9da475ee4945fe8`,
  `history-sum-control-round26/PROOF.md`, lines 3--31: the later transfer
  calculation is explicitly in the neutral total-charge sector and uses the
  same normalized hopping. The already-reviewed transfer-log result supplies
  the faithful Gibbs state of the limiting Round-20 Hamiltonian.

The Atlas search found no earlier contract containing this same thermal
Liouville/GNS compression theorem. The superficially similar
`contract:charged-sector-memory-gate-20260920` is a different result: it
compresses a proposed affine Hamiltonian \(H=aL_0+c\) to D8/spinor sectors
and finds vanishing Feshbach memory there. It neither proves nor contradicts
the present nonzero local-grade memory under Round-20 hopping.

## 1. The weighted GNS Hilbert space and Liouvillian

Let \(\mathscr H_0\) be the neutral physical Hilbert space and

\[
\rho_\beta=Z_\beta^{-1}e^{-\beta H}.
\]

The previous transfer review established that \(H\) is self-adjoint,
\(e^{-\beta H}\) is trace class, and \(\rho_\beta\) is faithful. Define

\[
\langle X,Y\rangle_\rho=\operatorname{Tr}(\rho_\beta X^*Y),
\qquad
\mathscr K_\rho=\overline{B(\mathscr H_0)}^{\|\cdot\|_\rho}.
\]

Faithfulness makes this a norm on bounded operators. The map

\[
J_\rho:X\longmapsto X\rho_\beta^{1/2}
\]

is an isometry into Hilbert--Schmidt space. Its range is dense: in an
eigenbasis of \(\rho_\beta\), every matrix unit is
\(J_\rho(r_j^{-1/2}|i\rangle\langle j|)\), and finite-rank operators are
dense in Hilbert--Schmidt norm. Thus \(\mathscr K_\rho\) is exactly the
weighted GNS/Liouville Hilbert space needed here.

Because \([H,\rho_\beta]=0\),

\[
\alpha_t(X)=e^{itH}Xe^{-itH}
\]

preserves the Gibbs norm. Under \(J_\rho\), it is ordinary two-sided
conjugation on Hilbert--Schmidt operators. It therefore extends to a strongly
continuous unitary group on \(\mathscr K_\rho\). By Stone's theorem

\[
\alpha_t=e^{it\mathcal L},
\qquad
\mathcal L=H_{\rm left}-H_{\rm right}
\]

is self-adjoint. On its natural commutator core it acts as
\(\mathcal L X=[H,X]\). This is the rigorous meaning of the displayed
commutator; no unrestricted multiplication of unbounded operators is used.

## 2. The four native grade projectors really exist

At a marked site \(x\), define the grade from the **physical** E8 vector,
not directly from an arbitrary coordinate convention:

\[
k(n_x)=2\sum_{i=1}^5(B_{E8}n_x)_i\pmod4
      =j(B_{E8}n_x)\pmod4.
\]

Let \(M_a\) be multiplication by
\(\mathbf 1_{\{k(n_x)=a\}}\), restricted to the neutral total-charge
sector. All four projectors are nonzero there. For example, at two sites the
configuration \((m e_8,-m e_8)\) is neutral and has local grade
\(m\pmod4\), since \(B_{E8}e_8=b_8=s\) and \(j(s)=1\).

The projectors are mutually orthogonal and sum to the identity. Since the
state is faithful,

\[
w_a=\operatorname{Tr}(\rho_\beta M_a)>0.
\]

Hence

\[
\mathscr S=\operatorname{span}\{M_0,M_1,M_2,M_3\}
\]

has dimension four. Its GNS-orthogonal projection is

\[
PX=\sum_{a=0}^3
M_a\,\frac{\operatorname{Tr}(\rho_\beta M_aX)}{w_a}.
\]

This formula also exposes the physical limitation: because hopping generally
makes \([\rho_\beta,M_a]\ne0\), the coefficient functional need not be
positive. Thus this orthogonal \(P\) is not automatically a positive or CP
conditional expectation.

There is a stronger exact obstruction. Let
\(\mathcal A=\operatorname{span}\{M_a\}\). If a normal conditional
expectation \(E:B(\mathscr H_0)\to\mathcal A\) preserved the Gibbs state,
then bimodularity would give, for \(a\ne b\),

\[
E(M_aXM_b)=M_aE(X)M_b=0.
\]

State preservation would therefore imply

\[
0=\operatorname{Tr}(\rho_\beta M_aXM_b)
 =\operatorname{Tr}(M_b\rho_\beta M_aX)
\quad\text{for every bounded }X.
\]

Since \(M_b\rho_\beta M_a\) is trace class, duality with all bounded
\(X\) forces \(M_b\rho_\beta M_a=0\). Thus \(\rho_\beta\) must be block
diagonal in the grade projectors, equivalently
\([\rho_\beta,M_a]=0\) for every \(a\). The scalar function
\(E\mapsto Z_\beta^{-1}e^{-\beta E}\) is injective, so its spectral
projections are the spectral projections of \(H\); hence this would also
force \([H,M_a]=0\). Section 5 gives the direct contradiction when \(J>0\).

Conversely, if \(\rho_\beta\) were grade block diagonal, the formula

\[
E(X)=\sum_a
\frac{\operatorname{Tr}(\rho_\beta M_aX)}{w_a}M_a
\]

would define a normal unital CP, bimodular, idempotent and
\(\rho_\beta\)-preserving expectation onto \(\mathcal A\). Therefore the
obstruction is exact: for the grade-mixing native Hamiltonian, no
Gibbs-preserving normal conditional expectation onto the four-grade algebra
exists. A sharp nonselective grade readout cannot simply be identified with
the equilibrium-preserving GNS projection.

## 3. Domain and block-operator audit

Every unbounded piece of \(H\) is diagonal in the site charges or independent
of them. Therefore

\[
[H,M_a]=[H_{\rm hop},M_a].
\]

At fixed spatial volume \(H_{\rm hop}\) is bounded, so the commutator on the
right is bounded. It follows that each \(M_a\) belongs to
\(\operatorname{Dom}\mathcal L\), and \(\mathcal LM_a\) is a valid vector of
\(\mathscr K_\rho\).

Put \(Q=I-P\) and, as operators between the finite resolved space and its
orthogonal complement,

\[
A=P\mathcal LP,\qquad B=Q\mathcal LP.
\]

Because \(\mathscr S\) is finite dimensional and contained in
\(\operatorname{Dom}\mathcal L\), both \(A\) and \(B\) are bounded and
\(A=A^*\). Moreover \(P\mathcal LQ=B^*\) on
\(\operatorname{Dom}\mathcal L\cap Q\mathscr K_\rho\), and the right-hand
side gives its bounded extension.

The operator

\[
\mathcal D=\mathcal L-(B+B^*)
\]

on \(\operatorname{Dom}\mathcal L\) is a bounded self-adjoint perturbation
of \(\mathcal L\). Its off-diagonal blocks vanish, so it reduces both
\(P\mathscr K_\rho\) and \(Q\mathscr K_\rho\). Consequently

\[
\mathcal D=A\oplus\mathcal L_Q,
\qquad
\mathcal L_Q=Q\mathcal LQ,
\quad
\operatorname{Dom}\mathcal L_Q=\operatorname{Dom}\mathcal L\cap Q\mathscr K_\rho,
\]

and \(\mathcal L_Q\) is self-adjoint. This validates the propagator
\(e^{it\mathcal L_Q}\); a merely formal compression of an unbounded
self-adjoint operator would not have sufficed.

## 4. Exact memory equation and the second derivative

For

\[
C(t)=Pe^{it\mathcal L}P\big|_{\mathscr S},
\]

write \(X(t)=Pe^{it\mathcal L}P\) and
\(Y(t)=Qe^{it\mathcal L}P\). The block equations are

\[
X'=iAX+iB^*Y,
\qquad
Y'=iBX+i\mathcal L_QY,
\qquad Y(0)=0.
\]

Variation of constants gives

\[
C'(t)=iAC(t)-\int_0^t
B^*e^{i(t-s)\mathcal L_Q}B\,C(s)\,ds,
\qquad C(0)=I_{\mathscr S}.
\]

This is an exact finite-dimensional Volterra equation. Since \(B\) has
finite-dimensional domain, strong continuity of
\(e^{it\mathcal L_Q}\) implies norm continuity of
\(B^*e^{it\mathcal L_Q}B\). The equation can therefore be differentiated at
zero within the finite matrix space:

\[
C'(0)=iA,
\qquad
C''(0)=-A^2-B^*B.
\]

No assumption \(\mathscr S\subset\operatorname{Dom}\mathcal L^2\) is
needed. Equivalently, for \(u,v\in\mathscr S\), the second derivative is
the quadratic-form identity

\[
\frac{d^2}{dt^2}\langle u,C(t)v\rangle\big|_{t=0}
=-\langle\mathcal Lu,\mathcal Lv\rangle.
\]

The equal-time kernel is

\[
K(0)=B^*B\ge0.
\]

This kernel determines the uninserted compressed propagator \(C(t)\), hence
the corresponding two-time projection problem. It is not by itself a full
law for arbitrarily many time-resolved grade insertions. A canonical Hilbert
space sufficient for all such prescribed bounded histories is instead

\[
\mathscr K_{\rm hist}
=\overline{\operatorname{span}}\left\{
L_{M_{a_n}}e^{it_n\mathcal L}\cdots
L_{M_{a_1}}e^{it_1\mathcal L}\mathbf 1
\right\},
\]

where all finite words, real times and grades are included and
\(L_{M_a}X=M_aX\). This is a subspace canonically generated inside the same
thermal GNS representation, not an added bath. The smaller cyclic subspace
of \(Q\mathscr K_\rho\) generated by \(\operatorname{Ran}B\) under
\(e^{it\mathcal L_Q}\) is enough to minimize the displayed two-time memory
realization, but it need not contain the vectors created by intermediate
left insertions. No full event law follows without controlling that larger
history space and specifying the operational instrument.

## 5. The original neutral hopping makes the memory nonzero

The coefficient-basis hopping directions have grade increments

\[
(0,0,0,0,0,2,2,1)\pmod4.
\]

Assume the retained hopping strength satisfies \(J>0\), as in the Round-20
Hamiltonian. In particular the existing \(p=e_8\) channel transports the physical
spinor vector \(s\) and changes the marked-site grade by one. Starting from
the all-zero neutral configuration, one bond term creates
\((-e_8,+e_8)\). The marked local grade changes from zero to three while
total charge remains zero. Round 20 gives a nonzero hopping matrix element
(magnitude \(J\), or the stated doubled finite-cell multiplicity at
\(L=2\)); diagonal terms cannot cancel it.

It follows that \([H,M_0]\) has a nonzero off-diagonal charge matrix element.
Every element of \(\mathscr S\) is diagonal in that basis, so
\(\mathcal LM_0\notin\mathscr S\). Faithfulness of \(\rho_\beta\) ensures
that a nonzero bounded operator cannot become the zero GNS vector. Hence

\[
B=Q\mathcal LP\ne0,
\qquad
B^*B\ne0.
\]

Therefore the compressed dynamics cannot equal \(e^{itA}\), because its
second derivative differs by the nonzero positive matrix \(-B^*B\).
At the separately allowed Round-26 boundary value \(J=0\), this witness
vanishes and no nonzero-memory conclusion follows from the hopping argument.

The kernel has a sharper exact rank. For

\[
F=\sum_{a=0}^3 f_aM_a\in\mathscr S,
\]

the operator \(\mathcal LF=[H_{\rm hop},F]\) has only off-diagonal charge
matrix entries. If \(BF=0\), then \(\mathcal LF=P\mathcal LF\) as a GNS
vector. Faithfulness upgrades this to equality of bounded operators. The
left side is off diagonal and the right side belongs to the diagonal span
\(\mathscr S\), so both must vanish. For every grade \(a\), an \(e_8\) hop
on the selected bond has a nonzero matrix element from grade \(a\) to
\(a+1\pmod4\). Therefore \([H_{\rm hop},F]=0\) forces
\(f_{a+1}=f_a\) for all four \(a\). Conversely constants are killed by every
commutator. Thus

\[
\ker B=\mathbb C I,
\qquad
\operatorname{rank}B=3,
\qquad
\operatorname{rank}K(0)=3.
\]

Accordingly \(K(0)=B^*B\) is positive semidefinite on \(\mathscr S\), with
the identity as its one-dimensional null direction. It is positive definite
only after restricting to the three-dimensional GNS-orthogonal complement
of the constants; calling it strictly positive on all four grades would be
wrong.

## 6. Exact four-grade/carry algebra on one native transfer direction

Choose one oriented bond \(x\to y\), take the existing \(p=e_8\) channel,
and mark the **arrival** site \(y\). This fixes the sign convention:
\(W=W_{xy,e_8}\) raises \(k(n_y)\) by one. Marking the departure site would
lower the grade and would require reversing \(W\) or the grade labels; it is
the same structure with the opposite convention.

Let \(P_r=M_r\) be the arrival-grade projectors and set

\[
Z=W^4,
\qquad
R=W(I-P_3)+W^{-3}P_3.
\]

The source's cocycle phases along any bilateral \(W\)-orbit can be removed
by recursively rephasing its basis vectors. More strongly, the displayed
identities can be checked directly from the grade-shift relations, without
discarding the cocycle:

\[
R^4=I,
\qquad
[R,Z]=0,
\qquad
W=R\bigl[(I-P_3)+ZP_3\bigr].
\]

Here \(R\) cycles only the four residues. On the closing step, \(Z\) carries
the integer orbit coordinate. Define

\[
E_{ab}=R^aP_0R^{-b},
\qquad 0\le a,b\le3.
\]

They satisfy

\[
E_{ab}^*=E_{ba},
\qquad
E_{ab}E_{cd}=\delta_{bc}E_{ad},
\qquad
\sum_aE_{aa}=I.
\]

The unitary \(Z\) commutes with these matrix units. Conversely, \(W\) and
all \(P_r\) are generated by \(Z\) and the \(E_{ab}\). On each bilateral
orbit, \(Z\) is the ordinary shift of the integer carry and has full spectrum
\(S^1\). Fourier transformation of the carry therefore gives

\[
C^*(W,P_0,P_1,P_2,P_3)\cong M_4(C(S^1)),
\]

and its von Neumann closure in the bilateral orbit representation is

\[
M_4(L^\infty(S^1)).
\]

The full neutral Hilbert space contains many such \(W\)-orbits. Their
cocycle phases can be gauged separately along the nonperiodic integer line,
so they give a multiplicity of the same faithful representation rather than
extra points of the spectrum. The phrase "carry center" refers only to the
center of this generated orbit algebra. It is not central in the full native
observable algebra.

The Atlas search returned the Round-16 carry contract but no registered
statement of the concrete \(M_4(C(S^1))\) orbit algebra. The algebra above is
therefore a new consequence proved from the registered source formulas, not
a quotation of an earlier Atlas verdict.

This is precisely the infinite carry already present in Round 16. Replacing
it by a finite periodic register would impose \(W^{4q}=I\) for some \(q\),
while the source has a nontrivial translation by \(4q e_8\). It would also
conflict with the exact local charge Ward relation: a finite-dimensional
unitary satisfying \(W^*n_yW=n_y+e_8\) would make the trace of the two sides
different. Equivalently, the exponentiated Ward relation
\(e^{i\theta n_y}We^{-i\theta n_y}=e^{i\theta}W\) would give, in dimension
\(d\), the determinant identity
\(\det W=e^{i\theta d}\det W\) for every \(\theta\), which is impossible.
Thus the \(M_4\) mark cannot be separated by quotienting away the
carry while retaining the original charge action.

## 7. Why the carry center is not a conserved dynamical fiber

Let \(r_{xy}\) be the neutral configuration created by one \(e_8\) transfer:
\(-e_8\) at \(x\), \(+e_8\) at \(y\), and zero elsewhere. Along the chosen
orbit, \(n_j=j r_{xy}\). The native diagonal charge energy is a positive
quadratic form, hence

\[
D(n_{j+4})-D(n_j)
=\bigl((j+4)^2-j^2\bigr)D(r_{xy})
=(8j+16)D(r_{xy}).
\]

Since \(D(r_{xy})>0\), the difference is nonzero already at \(j=0\). It
happens to vanish at the symmetric point \(j=-2\), which is not an operator
identity and does not restore conservation. The charge-dependent scalar
coupling has the analogous homogeneous quadratic difference
\((8j+16)V(r_{xy},\phi)\); the strictly positive charge term alone already
suffices.

The charge-independent scalar kinetic term commutes with \(Z\). The complete
\(H_{\rm hop}\) also commutes with \(Z=W^4\): same-bond transfer generators
commute, disjoint ones commute, and each adjacent-bond cocycle sign is raised
to the fourth power. Its retained diagonal constant therefore cancels in the
conjugation as well. Consequently the nonzero commutator \([H,Z]\) comes
from the native diagonal charge/scalar energy, not from dropping a hopping
offset.

It follows that the Fourier parameter \(z\in S^1\) central for the orbit
algebra does not define a reducing fiber of the full Hamiltonian. Other hop
directions also move between the transverse labels of the chosen \(W\)-orbits.
Thus the algebraic decomposition supplies a native four-mark dictionary and
its mandatory carry, but no source-derived factorization of the entire
dynamics into an independent four-level system and environment.

Finally, this constant \(M_4\) subalgebra is abstractly a four-level matrix
algebra. Nothing in the construction canonically identifies it with the
earlier profile-free \(\Phi\), \(\Gamma\), or any family-carrier \(M_4\).
Such an identification would require an explicit intertwiner of the actual
operators and dynamics. The grade functional itself remains an inherited
input rather than a P1 source derivation.

## 8. Full neutral-lattice extension and unavoidable Gibbs correlation

The neutral integer charge lattice on \(N\) sites is

\[
\Lambda_0=
\left\{(n_x)\in(\mathbb Z^8)^N:\sum_xn_x=0\right\}
\cong\mathbb Z^{8(N-1)}.
\]

For the selected arrival site \(y\), lift its local grade to the integer
functional

\[
\jmath(n)=j(B_{E8}n_y).
\]

For the bond-transfer vector \(v=r_{xy}\), one has \(\jmath(v)=1\). Thus
\(\jmath\) is primitive and the map \(q\mapsto q-\jmath(q)v\) gives the
exact internal direct sum

\[
\Lambda_0=\ker\jmath\oplus\mathbb Zv,
\qquad
\ker\jmath\cong\mathbb Z^{8(N-1)-1}.
\]

Writing \(\jmath=4m+r\), \(0\le r<4\), yields the Hilbert-space coordinate
isomorphism

\[
\ell^2(\Lambda_0)\otimes L^2(\mathbb R^N)
\cong
\mathbb C^4_{\rm grade}\otimes
\ell^2\!\left(\mathbb Z_{\rm carry}\times
\mathbb Z^{8(N-1)-1}_{\ker\jmath}\right)
\otimes L^2(\mathbb R^N).
\]

Unchanged spectator species may be appended to the last factor. The cocycle
along the \(v\)-coordinate is removed by the same recursive unit phases as
in Round 16. This is an exact coordinate split of the existing Hilbert
space, not an appended environment. The choice of marked site, primitive
functional and kernel complement remains a convention inherited from the
chosen grade.

This split permits a four-vector test of dynamical factorization. Hold the
kernel coordinates at zero and choose one normalized Schwartz scalar vector
\(\psi\). The product basis vectors with

\[
(m,r)=(0,0),(0,1),(1,0),(1,1)
\]

have physical orbit coordinates \(q=0,1,4,5\), respectively. Their diagonal
Hamiltonian expectations have the form

\[
E(q)=q^2A_\psi+E_{\rm scalar}(\psi)+48JN,
\]

where

\[
A_\psi=D(v)+\frac{\lambda}{N}
\left\langle\psi,
e(e_8)(\phi_x^2+\phi_y^2)\psi\right\rangle>0.
\]

The precise nonnegative scalar-coupling term may vanish when \(\lambda=0\),
but \(D(v)>0\) always remains. The hopping shifts have zero diagonal matrix
elements in a fixed charge basis; their retained constant contributes the
same \(48JN\) to all four vectors and cancels from the mixed difference.
Therefore

\[
E(5)+E(0)-E(4)-E(1)=8A_\psi>0.
\]

If the full Hamiltonian in this coordinate split were

\[
H=H_{\rm grade}\otimes I+I\otimes H_{\rm memory},
\]

the same rectangular difference of product-vector expectations would be
zero. The displayed native-source witness therefore excludes that tensor
sum. It also shows why retaining the carry is more than bookkeeping: the
quadratic source energy couples grade and carry even before considering
other hop directions.

Now suppose the faithful Gibbs density were a product in this same split,

\[
\rho_\beta=\rho_{\rm grade}\otimes\rho_{\rm memory}.
\]

Faithfulness makes both factors injective. Since the grade factor is finite
dimensional, \(\log\rho_{\rm grade}\) is bounded; the possibly unbounded
\(\log\rho_{\rm memory}\) is defined by spectral calculus on its natural
domain. Spectral calculus for tensor products gives

\[
\log(\rho_{\rm grade}\otimes\rho_{\rm memory})
=\log\rho_{\rm grade}\otimes I
+I\otimes\log\rho_{\rm memory}.
\]

Because \(\rho_\beta=Z_\beta^{-1}e^{-\beta H}\) and the exponential is an
injective scalar function on the real spectrum, taking the spectral logarithm
would force \(H\) to be such a tensor sum, up to the scalar
\(-\beta^{-1}\log Z_\beta\,I\). This contradicts the positive rectangular
difference. Hence the original faithful Gibbs state is necessarily
grade--memory correlated in this native coordinate split.

This does not select a resetting map or a preferred state for the carry. It
shows the opposite: replacing the native carry/kernel/scalar degrees by an
independent fixed memory state is not licensed by the original equilibrium
state. It is not a no-go theorem for all alternative Hilbert-space
factorizations or all possible instruments.

## Exact scope

What is proved is an exact memory identity for a chosen four-dimensional
subspace of observables under the already-existing native Hamiltonian and
its faithful equilibrium state. The complement is internal to the same
closed quantum system; no external bath or ad hoc dissipator was added.

What is not proved:

* that \(P\) is a CP map, a Lüders operation, or a physical measurement;
* that any normal conditional expectation onto the four-grade algebra can
  preserve the thermal state (the opposite is proved for \(J>0\));
* that the unresolved complement is actually discarded by the source;
* that the local Z4 grade is selected as the physical readout rather than
  inherited as an earlier arithmetic choice;
* that the constant \(M_4\) orbit algebra is canonically the \(M_4\) of an
  earlier channel, family carrier, or instrument;
* that a fixed Fourier carry value is preserved by the full Hamiltonian;
* that grade and carry/kernel/scalar state are independent in the explicit
  native coordinate split (the Gibbs product is excluded there);
* Markovianity, exponential decay, a short-memory approximation, or a
  positive classical stochastic kernel;
* an arbitrary multi-time event law from the two-time compression alone;
* a spatial continuum, Lorentzian causality, or a microscopic TFPT origin
  for the Round-20 hopping rates.

Thus the theorem is a genuine constructive consequence of the native
Round-20--26 dynamics, but its physical interpretation remains conditional
on a source-derived choice of readout or instrument.
