# Native charge conjugation is an all-time charged-response discriminator

15 October 2026. `UR.NATIVE.CHARGED_CONJUGATION.01`. **PARTIAL** for the
physical common-source question. The following source-specific symmetry and
response statements are exact under the stated original Hamiltonian/state.
They do not select that Hamiltonian, a chemical marking, or a P1/P2 source.

## 1. Start with the original full field net

`local-charge-transport-round20/PROOF.md`, sections 1–5, constructs

    H_charge = tensor_x l2(Z8),
    T_(x,p)|{n}>=epsilon(p,n_x)|...,n_x+p,...>,
    W_(xy,p)=T_(y,p) T_(x,p)^*,
    epsilon(p,n)=(-1)^beta(p,n),
    beta(p,n)=sum_i p_i n_i G_ii/2+sum_(i>j) p_i n_j G_ij mod2.

This uses the complete original even unimodular G, not a parity truncation.
All integer charges remain. The onsite operators and diagonal projectors
generate B(l2(E8)) at each site. The source has global U(1)^8 charge
Q_tot=sum_x n_x. Only the finite frame group is gauged; fixing the later
sector C0={sum_x n_x=0} does not turn all charged onsite fields into
endomorphisms of that sector.

For the later scalar-coupled source in
`native-source-reconstruction-20261014/PROOF.md`, retain its given finite
spatial lattice, scalar space, couplings and domain:

    H = K_scalar+V(n,phi)+H_hop,
    V = D(n)+phi^T[-Delta_a+m²I+2u diag(e(n_x))]phi/2,
    e(n)=n^T G n/2,
    D = sum_x e(n_x)/N+nu sum_edges e(n_x-n_y),
    H_hop=J sum_(edges,p)(2I-W_p-W_p^*).

The original coercive diagonal/scalar Hamiltonian, bounded hopping and
their common selfadjoint finite-volume realization are reused. The source
premises include a,m>0, u,nu,J>=0 and the specified eight arithmetic hopping
directions. They are not derived afresh here.

## 2. The actual cocycle permits charge conjugation

Define the linear unitary, not antiunitary,

    C|{n}>=|{-n}>,

acting as identity on the scalar Hilbert space. C²=I. Bilinearity gives
epsilon(-p,-n)=epsilon(p,n) and epsilon(p,-n)=epsilon(-p,n), hence

    C T_p C=T_-p.

One sign must be retained. The original adjoint law is

    T_p^*=epsilon(p,-p) T_-p=(-1)^e(p) T_-p.

In the original arithmetic basis the eight e(e_i) are 2,1,1,1,1,1,1,1.
Thus T_-p and T_p^* differ by a minus sign for seven basis directions;
they cannot be identified literally. This sign cancels in the adjoint-paired
responses below. On a bond the two onsite signs cancel as well:

    C W_p C=W_-p=W_p^*.

Quadratic evenness gives D(-n)=D(n) and V(-n,phi)=V(n,phi). Scalars are
unchanged; every hopping term is individually C invariant. The exact
Hamiltonian therefore satisfies C H C=H. C preserves the diagonal/scalar
form domain and the original core, while the hopping is bounded; this is
an identity of the existing selfadjoint operator, not only a finite matrix
or formal Taylor coefficient. It implies C exp(-itH) C=exp(-itH) and the
same for bounded positive-time heat propagators.

Also C Q_tot C=-Q_tot, so C preserves C0. For any beta>0 for which the
original trace-class theorem applies, both the full Gibbs density and the
normalized sector density

    rho_beta=exp(-beta H)/Z,
    rho_beta,0=P0 exp(-beta H)P0/Z0

are C invariant. The second expression is a density on the full Hilbert
space supported in C0; the full generator is still needed for charged
intermediate propagations.

## 3. Equality of complete adjoint-paired charged responses

For any C-invariant density rho, define the bounded Euclidean return response

    G_(rho,p)(t)=Tr[rho T_p^* exp(-tH) T_p], t>=0.

The original H is nonnegative. Unitary conjugation inside the trace gives

    G_(rho,p)(t)=G_(rho,-p)(t)

for every t, not just at t=0 or to one perturbative order. The same argument
applies to the real-time ordered response
`Tr[rho T_p^* exp(itH) T_p exp(-itH)]`. C is linear, so no complex
conjugation of the amplitude or reversal of t is introduced.

More generally any specified bounded multitime word has the same expectation
as its C-conjugate word. Thus the full source symmetry, rather than a
single matching moment, constrains every sign-reversed marked history.
For Euclidean KMS correlation functions the same result holds in their
standard thermal analyticity domain; no unbounded exp(+tH) is introduced
in the definition used above.

This result requires a C-invariant starting state. A non-C-invariant
prepared state, a C-breaking boundary condition/field, or another Hamiltonian
changes the premise. It is not a theorem that every possible native vacuum
or infinite-volume phase must be C invariant.

## 4. Decisive comparison with one fixed chiral quarter twist

In the separately stipulated chiral E8 lattice-Fock representation, the
fixed primary vacuum and its charge-shift preparations yield

    C_delta(v,t)=exp[-t(|v|²/2+delta.v)],
    C_delta(v,t)/C_delta(-v,t)=exp[-2t delta.v].

Here delta_k=(k+alpha)/8-1/2, with alpha=1/4 or 3/4 in the original quarter
cases. For the *same marked* opposite-charge pair, the native response at
the invariant start is symmetric whereas this fixed-twist response is
asymmetric whenever delta.v is nonzero.

The checker uses the actual original native embedding basis B, not a root
count or a dimension match. Its eight columns span R8. Therefore a nonzero
delta cannot annihilate every B e_i. For example B e_0=2e_0 has delta pairing
-15/16 at alpha=1/4 and -13/16 at alpha=3/4. At least this charged pair
already discriminates the fixed-twist target from the symmetric native
response data.

Consequently an exact identification preserving **all** these marked
responses, their opposite-charge assignment, the specified starting states
and a common time law is excluded. This is not an exclusion of every
unmarked partial readout, state-changing preparation, different source,
or covariant family of twists. A symmetric mixture of +/-delta can remove
this particular asymmetry, but is not thereby a native source dictionary
or a physically selected state; no such mixture is inserted as a solution.

## 5. C0 and uniform chemical marking: keep the charged corners

A uniform linear marking has the form eta.Q_tot, with eta in the *native
integer-coordinate dual*. It vanishes exactly on C0. Because [H,Q_tot]=0,
the formal comparison H_eta=H+eta.Q_tot restricts to the same H_0 on C0.
It also obeys C H_eta C=H_-eta. The coercive original charge energy controls
this linear marking, so its finite-volume selfadjoint comparison is
well defined. Eta is additional input here, not selected by P1/P2.

Balanced W transports remain within C0, commute with Q_tot and have the
same sector-internal dynamics under H_eta. They still move local charges
and have the original nontrivial cocycle and dynamics. C0 never implies
that local charge information or all neutral responses are absent.

By contrast T_(x,p) maps C0 to C_p. One cannot use it as an endomorphism
inside H_0 or silently delete the exterior propagation. With the full
charge-sector generator and a normalized C-invariant start rho_0 supported
in C0, exact sector conservation gives

    G_(rho_0,p)^eta(t)
      =Tr[rho_0 T_p^* exp(-tH_eta) T_p]
      =exp[-t eta.p] G_(rho_0,p)^0(t).

The undeformed pair is equal by section 3. Since the heat operator is
strictly positive on nonzero vectors, these normalized returns are positive,
and their ratio is exactly

    G_(rho_0,p)^eta(t)/G_(rho_0,-p)^eta(t)=exp[-2t eta.p].

Thus the uniform marking is invisible to the dynamics confined to C0 but
visible in onsite insertion/removal responses with charged intermediate
sectors. For a proposed geometric label map v=Bp, the coefficients would
have to satisfy eta=B^T delta. This is a conditional coordinate relation,
not a derivation or selection of eta. It matches only this odd charged
response factor, not the full even energy, oscillator response, state,
locality or interacting field algebra.

This charged-corner distinction is the positive reusable interface. An
automated search must distinguish sector-internal W data from full onsite
T insertion data and must not claim a chiral twist has been selected merely
because a new linear chemical coefficient can be stipulated.

## 6. Exact check and scope

The checker passed **24/24 exact groups**, zero failures. It loads the
original Round20 Gram and actual beta/translate/transfer functions by AST,
without executing the old prerequisite suite or replacing the cocycle.
Finite monomial-action samples include a large-charge neutral profile;
they are implementation regressions of the all-charge bilinear identities.
The analytic C symmetry and functional-calculus argument prove the all-time
statement. No finite heat matrix, partition sum, Gibbs spectrum, oscillator
cutoff, new physics simulation or independent empirical test is claimed.

Pinned original sources are in `source_manifest.json`:

- `local-charge-transport-round20/PROOF.md`, sections 1–5: full field space,
  onsite cocycle/adjoint laws, neutral transport, global U(1)^8 and H_charge.
- Its original `checker.py`, main Gram assignments and local field functions.
- `native-source-reconstruction-20261014/PROOF.md`, sections 1–4: later C0
  restriction, scalar-coupled Hamiltonian, trace-class Gibbs family and time.

No ledger/paper/verification promotion or physical T1–T8 closure occurs.

## 7. Relation to existing charged-response work

`source-ground-response-20260920` already establishes a stronger all-ground
response comparison for a different pair of supplied constructions: the
native N=64 Hamiltonian `Delta N_b+g(b^dagger P+P^dagger b)` and a specified
quarter-shifted Vaux boundary. It admits all boundary ground mixtures and
coherences and finds incompatible charged thresholds under common positive
time calibration. That result is preserved. It is not the Round20 local
integer-charge/scalar source used in the present proof.

`source-charged-extension-selection-20260921` already proves that an added
odd-block energy can be invisible to the whole retained even dynamics while
changing charged two-point responses; full Gibbs sector weights can detect
it under their stated conditions. The chemical-corner distinction above
uses the same established logical separation. Its specific original
realization here is the actual Round20 U(1)^8 field net and its full bilinear
cocycle, with C and the original quadratic/scalar/hopping Hamiltonian.

The present contribution is therefore this concrete original-source
comparison and its exact C symmetry, not a newly discovered general
even/charged ambiguity or a conflation of the different native models.
It keeps the positive onsite T fields and their chemically sensitive
charged corners available for the next original-seam operator test.
