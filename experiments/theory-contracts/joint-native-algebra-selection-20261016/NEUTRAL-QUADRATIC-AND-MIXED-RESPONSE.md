# Neutral quadratic source fields and their first mixed response

2026-10-16. Continuation of UR.NATIVE.JOINT_OPERATOR.01, scoped/PARTIAL.
No physical gate, ledger, paper or load-bearing verification status changes.
The existing E8_1 vacuum, actual native code, and the already constructed
original neutral source fields are premises. No fitted coupling is added.
Use R_neut for the weight-two source field; R5 is the marked logical clock.

## 1. The actual source candidate, not a new Hamiltonian

Original `neutral_source_response._native_neutral_branch` starts with the
64 chiral source roots and omits one actual SU4 family weight, retaining48.
Its symmetrized normal products satisfy the original lattice identity

    B_alpha = (:E_alpha E_alpha^dagger: + :E_alpha^dagger E_alpha:)/2
            = (1/2) :(alpha.h)^2: .

The two ordered regular terms have opposite alpha(-2) contributions;
these cancel in the symmetric product. Thus

    B48 = 12 T_E8 + R_neut,
    R_neut = T_AR,  T_A=(1/2):h^T A h:,
    AR = 4 P_A2 - 8 P_U1,
    tr AR=0, tr AR^2=96, ||R_neut Omega||^2=48.

These coefficients are the ORIGINAL charge moments, not adjusted to the
logical clock. The source selection of a family omission is still an
input: the original reader accepts `omitted_family`. The existing required
map to a physically selected three-family space is not discharged here.
For each of the four choices the same trace/norm holds; their AR matrices
sum to zero. The full64 field is B64=16 T_E8. This identity is a control,
not an imposed physical averaging law. Its 5+2+1 eigenspaces live in the
eight-current source; they are not the logical five-code's 3+2 split.

## 2. Classify all local neutral weight-two zero modes at once

In the GIVEN lattice VOA, the full lattice-charge-zero weight-two space is

    h(-2)Omega + Sym^2(h(-1))Omega, dimensions 8+36.

This is lattice charge neutrality, not an assumed Standard Model singlet
condition. For an energy-preserving zero mode use o(v)=v_(wt(v)-1):

    o(h(-2)Omega)=-h(0),
    o(T_A)=sum_(n>0) h(-n)^T A h(n) + (1/2)h(0)^T A h(0).

For a Hermitian field A is real symmetric in a real orthonormal Cartan
basis. On the neutral oscillator code h(0) vanishes. All code states use
four mode-one oscillators, so only n=1 acts and o(T_A)=dGamma(A).
These identities concern the local field zero mode. The earlier individual
mode roundtrip K48=sum J_(lambda,-1)^dagger J_(lambda,-1) is NOT o(B48).

Let Pi be the actual five-code projector, f_j its normalized basis, and
J the actual native real complex structure. For any real symmetric A set

    Ac=(A-JAJ)/2,   Aa=(A+JAJ)/2.

Then the ACTUAL code has compression and total leakage

    Pi dGamma(A) Pi = (tr A)/2 I5,

    L(A)=sum_j ||(I-Pi)dGamma(A)f_j||^2
        =4||Ac-(tr A)I8/8||_HS^2 + (5/2)||Aa||_HS^2.

Derivation: split the one-particle space H+ plus H-. Write B=A|H+ and
C=P_H- A|H+. Actual one-register erasure gives Pi dGamma(B)Pi=tr B I5.
The original two-register marginal is (I+Swap)/4; expanding the squared
four-register sum gives code leakage 8||B-tr B I4/4||^2. The term C creates
one H- particle and has orthogonal leakage 5||C||^2 by the original
one-body number matrix 5I4. Reality of A gives

    ||Ac-(tr A)I8/8||^2=2||B-tr B I4/4||^2,
    ||Aa||^2=2||C||^2.

The real leakage Gram relative to the MATRIX Hilbert-Schmidt metric has
eigenvalues 0 (multiplicity1), 4 (multiplicity15), 5/2 (multiplicity20).
Hence its kernel is exactly A=aI8. The complex sesquilinear Gram also has
an imaginary antisymmetric part; it need not vanish. The displayed norm
formula uses REAL coefficients. Both parts are retained in the certificate.

The checker constructs all five ORIGINAL support polynomials in the eight
real oscillator coordinates, computes their Fock coefficients exactly,
checks all36 logical matrix actions and the full leakage Gram, and
compares it with this independently derived formula. It does not infer
the result solely from a numerical matrix rank or substitute a new code.

## 3. What the actual R_neut and B48 do on this code

For every common real orthogonal Cartan chart,

    Pi o(R_neut) Pi=0,
    240 <= L(AR) <= 384,
    Pi o(B48) Pi=48 I5,  L(A48)=L(AR),
    o(B64)|Code=64 I5,  L(A64)=0.

The inequalities are universal bounds, not a claimed sharp set of values
as source alignment varies. They use tr AR^2=96 and coefficients 5/2,4.
They therefore need no fitted Cartan identification. Concrete individual
matrix elements in different source charts still require the actual map.

In particular R_neut is NOT zero as a source field: its vacuum norm and
correlators remain as previously proved. Its first action on code leaks
out instead of giving a nontrivial INTERNAL logical Hamiltonian.
Pi here describes a mathematical restriction; repeated physical projection
or a Zeno/measurement rule is neither assumed nor derived.

There is also a global statement. Arbitrary time-dependent sequences
generated SOLELY by these Hermitian neutral weight-two zero modes preserve
each oscillator mode's occupation. On the initially mode-one-only code
their evolution is Sym^4(U8) for some one-particle unitary U8. If the final
operation preserves the code, it preserves its summed one-body matrix
5 P_H+. Thus U8 preserves H+; its restriction is a code-preserving U4.
The complete original stabilizer theorem then bounds its logical action
by phase*rho5(G31), which excludes marked R5. Intermediate exit from H+
and code is allowed in this argument. No physical time selection follows.

## 4. Do NOT combine separate exclusions: derive the mixed response

Root-current zero-mode rotations and neutral quadratic zero modes need
not remain in either separate group. In particular conjugating R_neut by
a root operation produces lattice-charged weight-two fields. The previous
inner-E8 exclusion and the new neutral-quadratic exclusion do NOT exclude
the combined source. This is the next justified connection to inspect.

Its first actual response already follows from the SAME lattice formula.
For an original norm-two lattice root beta (beta^2=2),

    E_beta,0 T_A = (v_beta.h(-1)) e^beta,
    v_beta = -A beta + (beta^T A beta) beta/2,
    beta.v_beta=0.

Indeed translate (1/2)x^T A x by -beta/z in the lattice vertex formula.
The z^-1 coefficient has -A beta.x from the first shift and
(beta^T A beta) beta.x/2 from creation times the complete shift. No
field identification, fitted coupling or new response rule enters.

For the original AR and ALL240 original doubled Chevalley roots r,
beta=r/2 (not the Construction-A r/sqrt2 convention), the exact distribution is

| squared norm of v_beta | root directions |
|---:|---:|
| 0 | 46 |
| 14 | 96 |
| 30 | 32 |
| 32 | 60 |
| 64 | 6 |

Thus194 actual first root responses are nonzero. These are complex
current-basis response directions, not194 elementary particles, independent
predictions or physically selected primitive instruments. Compact real
rotations pair the adjoint root responses.

Each nonzero response also has an exact action test on the actual code.
Let b=U^T beta, c=U^T v_beta in the original positive polarization. The
energy-preserving zero mode of (v_beta.h(-1))e^beta has, on f(x)Omega,
the PBW partition-(3) component

    [f(b) v_beta(-3) - (D_c f)(b) beta(-3)/3] e^beta.

Normal ordering gives the first term from v_beta(-3) times the completely
shifted quartic, and the second from the v_beta(1) derivative of f followed
by the completely shifted cubic and beta(-3) creation. Their signs differ
because the two shifted polynomial degrees have opposite parity.
The actual five quartics have no common nonzero projective zero; for every
real nonzero beta at least one f_j(b) is nonzero. Since v_beta is nonzero
and perpendicular to beta, these two vector contributions cannot cancel.
Every one of the194 individual nonzero directions therefore leaks from
code, and its compression is zero by its nonzero lattice charge. Existence
of leakage is invariant under common real orthogonal Cartan transport.
This is a FIRST-ORDER result, not a no-go for second responses or histories.

## 5. Actual continuation, scopes and review

The combined source now has an explicit derivative with coefficients fixed
by original AR and beta, rather than a freely chosen matrix pulse. The next
decisive object is an ORIGINAL ordered charged-pair return response on the
actual code and erasure carrier, together with its source-selected history.
An internally nontrivial return matrix alone would establish mathematical
availability, not physical permission, polar selection, a time scale or
unconditional 3+2. An inserted projection, detector, parent extension to
the leaked field sectors or fitted coupling cannot supply the missing law.

Preserve these exclusions exactly: neutral local weight-two zero modes
alone cannot give R5; individual first mixed responses leak. Combined
root/quadratic histories, individual K48 mode roundtrips, higher-weight
fields, other mode products and non-automorphism processes remain outside
the global exclusion. The prior36-direction response tomography remains
valid; it is not a claim that executing those fields yields the clock.

Original sources: neutral_source_response.py (_native_neutral_branch,
build_neutral_source_response_data), process.py (_code_basis), the prior
actual Cartan/Fock embedding and complete code stabilizer in this contract.
Primary formula: Dong and Nagatomo, Section2.2, pp4-5,
https://arxiv.org/pdf/math/9808088#page=5 .
The existing clock-jordan-response-reconstruction-20261015 provides all36
neutral quadratic response directions; it did not establish this code action.

Independent analytic review confirmed the full44-dimensional neutral
weight-two classification, zero-mode formula, leakage constants and bounds,
global endpoint-stabilizer argument, the source194-count and first mixed
PBW formula. The exact checker is replayed with Python -OO; checks do not
depend on assert statements. No physical field/coupling selection is certified.
