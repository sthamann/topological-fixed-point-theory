I have read all named sources (and located that `CHIRAL-SEWING-SELECTION.txt` and `CLAUDE-FIRST-REPORT.md` sit in the request folder, not under `R/`). Here is the corrected report.

---

# Corrected report: charged endpoint vs chiral sewing

**Scope.** Original-source analysis over the named files only. The first report (`CLAUDE-FIRST-REPORT.md`) is treated as an external proposal, not evidence. I checked each correction against the actual algebra and the cited lines.

---

## A. The six corrections, accepted or rebutted

### 1. `σ = W_ext^{1/2}` and "half-integral powers → Ramond" — **ACCEPTED; retracted.**

`SAME-SOURCE-AFFINE-DEFECT.txt:244-253` states the sector action of `W_ext` explicitly:
`0↔v`, `s↔c`, hence `0+s → v+c` labelled "**E8 vacuum representation to g-twisted representation**", with `0+v` (NS) and `s+c` (R) each preserved. So `W_ext` does **not** map into the Ramond sector; it maps the bosonic E8 module to the **g-twisted** module `v+c`.

The phrase "half-integral powers" refers to powers of the **vertex coordinate z**, not an operator root. Lines 321-325: "For D8 roots the exponent [`z^{-(Rβ).e0}`] is integral. For E8 spinor roots it is half-integral… the lattice realization of the Li Δ(−e0,z) twist," and reading it "as an ordinary local automorphism… is false." This is the half-monodromy in `z`, confirmed again at `AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:184-192` ("choose one branch of the fractional power").

`NEUTRAL-NET-AND-DEFECT.txt:214-222` fixes the module identities that the first report conflated:
- `H_{E8,g-twisted} = H_c ⊕ H_v` (the g-twisted module),
- `H_{F,R} = H_s ⊕ H_c` (the fermionic Ramond module).

These are different spaces. The first report (§2.2, §0, line 72) wrote that `W_ext`'s "half-integral powers land in the `v+c` (Ramond) representation" and set the generator `σ = W_ext^{1/2}`. **Both are wrong:** `v+c` is g-twisted, not Ramond; and no operator square root is given by the source. A bounded Borel function of the unitary `W_ext` commutes with every spectral/NS–R projection, so it cannot be a sector-changing spin field; a non-functional root is extra data and still proves nothing about locality. **I retract `σ=W_ext^{1/2}` and the identification of `v+c` with Ramond.**

### 2. A bounded `σ ∈ F(I)` in the NS CAR net — **ACCEPTED; this is the decisive retraction.**

The first report's premise (N) asked for a **bounded, `F(I)`-affiliated** `σ` in the *same* NS conformal net, claiming "only localization remains." This is self-contradictory. The first report itself requires `ψ^a(e^{2πi}z)σ(0) = −ψ^a(z)σ(0)` (its own §3, item 1): `σ` has **−1 monodromy with the vector fermions `ψ^a`**, which are the very generators of the NS net `F(I)`. An object with non-trivial relative monodromy to the generators of `F(I)` is **not** a local element of `F(I)`. Rescuing it as "local w.r.t. the full D8 currents `J^{ab}=:ψ^aψ^b:`" is a different, weaker statement: **D8-current-relative locality is not locality with the vector fermions.**

This is exactly the retained central/spinoriality obstruction (the `Spin(10)` kernel acts `−1` on the `16`), not merely a missing localization step. The correct object is an **unbounded charged intertwiner** between `H_{F,NS}` and `H_{F,R}` on the extended (spin) Hilbert space — a Ramond twist/disorder field of positive weight `h=1`, hence a smeared **unbounded** field, not a bounded element of the original `F(I)`. `AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:189-199, 276-279` says precisely this: the image "lies in the g-twisted module," and constructing a smeared intertwiner with strong convergence and a **bounded conformal-net defect** "still require energy bounds and a defect localization theorem… The Δ proposition by itself is a VOA module theorem." Energy bounds change neither sector nor monodromy.

**Correct formulation of the extra input.** The generator lives on an **extension net / field bundle**: the fermionic spin net `A_F` with its Ramond sector, equipped with (a) the defect/Arf lift (see point 3 below, `AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:201-238, 281-292`, "PHASE CONDITIONAL") and (b) an energy-bounded localized charged intertwiner `H_{F,NS} → H_{F,R}`. The actual extra physical input is therefore a **chosen spin/Ramond extension plus a defect-localization theorem**, not a bounded `σ∈F(I)`. **I retract "only localization remains"** and the claim that the residual gap is a single `F(I)`-boundedness check.

### 3. Clifford normalization, `2^8`, and Spin(16) — **ACCEPTED; corrected.**

With `{ψ_0^a, ψ_0^b} = δ_{ab}` (as in `SAME-SOURCE-AFFINE-DEFECT` conventions), the Dirac matrices are `γ_a = √2 ψ_0^a`, and the chirality is
`Γ = i^8 γ_1⋯γ_16 = i^8 (√2)^{16} ψ_0^1⋯ψ_0^16 = 2^8 ψ_0^1⋯ψ_0^16`.
The first report's `Γ = i^8 ψ^1_0⋯ψ^16_0` (§2.1) **dropped the `2^8`.** The irreducible module is `2^8 = 256`-dimensional (correct in the first report), splitting `128_+ ⊕ 128_-`. The spinor representation is carried by **Spin(16)**, not `SO(16)`; the first report's "`SO(16)` spinor" is loose. The branching `128_+ = (16,4)⊕(\overline{16},\overline{4})` under `so(10)⊕so(6)` and `h=1`, `dim=64` are correct and independently confirmed (`:209`, `NEUTRAL-NET-AND-DEFECT.txt:242`).

Separately: "**one-dimensional intertwining spaces**" (`AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:243`) are the **fusion Hom-space multiplicities** (`n_a∈{0,1}`, `SIMPLIFICATION.txt:66-79`), **not** the Ramond ground space, which is 64-dimensional. The first report's phrase "1-dim Ramond ground intertwiners" (§2.2, §7) conflated the two. **Corrected.**

### 4. "Categorical/theta route is a proven type obstruction" — **PARTLY REBUTTED; blanket refutation retracted.**

The dimension facts are right: `WEAKEST-DIRECT:349-363` gives `W_aux` of complex dimension 4 (a determinant/Heisenberg line), while `M_h` has a 64-dim ground space and infinitely many descendants. This excludes **direct identification** of `W_aux` with `M_h` — a genuine type mismatch. But that does **not** refute a *functor* or categorical *reconstruction* from finite geometric data to an infinite VOA, which is a different claim. The first report elevated the no-direct-identification into "an already-present categorical boundary completeness principle *cannot* be the origin," which overreaches. Indeed `WEAKEST-DIRECT:305-315` records a stronger auxiliary branch (theta → primitive integral E8 lattice via the 60 stabilizer rays → even unimodular rank-8 lattice), i.e. exactly a finite-data → lattice reconstruction, with an explicit arithmetic normalization premise (the `ζ8 G31` homothety ambiguity). **I retract the blanket refutation.** The accurate statement: the theta/determinant line is not *itself* the charged generator `M_h`, but a categorical reconstruction functor to the E8 VOA is not excluded and is an open alternative route (its premise being arithmetic normalization, not a dimension count).

### 5. "Chiral subnet discards half the D8 currents" — **ACCEPTED; reasoning retracted.**

The first report's Route-2 argument (§5) claimed switching to a chiral current subnet "means discarding half the D8 currents." This is wrong. `NEUTRAL-NET-AND-DEFECT.txt:220-222`: fermionization yields **16 free Majorana–Weyl (single-chirality) fermions**, whose currents `J^{ab}=:ψ^aψ^b:` generate the **full** `so(16)_1 = (D8)_1` (all 112 roots). A single chiral 16-Majorana theory already carries all of D8; switching spacetime chirality does not delete half the internal Lie algebra. So the stated reason for rejecting Route 2 is void. (Whatever is retained about the *vector-Ward* no-go and the nonchiral Dirac charges requires an explicit source map, `NEUTRAL-NET-AND-DEFECT.txt:371-392`; it cannot be transferred by name.) **I retract the "discards half the D8" reasoning.**

### 6. A3 unit spectral flow, deck `Z2`, Ramond puncture — **ACCEPTED; retracted.**

`ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:640-654`: the NS matching path gives `λ_n^{NS}(t)=n+1/2+t` with **one crossing** over `t:0→1` (`φ:0→2π`), i.e. full unit flow returns NS→NS and changes occupation by one unit. Half flow (`π`) is what exchanges NS/R. The first report (§6) called the `φ:0→2π` flow "the NS↔R interchange a spin-field insertion performs" — **wrong**; unit flow preserves the spin boundary condition. A3 and A4 both remain explicit and compatible (readings 1/2, `:657-670`).

The deck action is **free**: `tau(s,t)=(s+2π,−t)` on `M̃=([0,4π]×[−ε,ε])/∼` (`:548-554`), with the lifted core `t=0` carrying `Hol_spin = r^2 = −1` ⟹ **NS** (`:563-567`, 11.2). There is no "open fixed point of the deck `Z2`." The first report's `σ` "at an open fixed point of the deck `Z2`" (§4, §6) is unsupported — a Ramond puncture is an **additional geometric defect** unless the original prescription supplies it, which `:740-752` ("A4 still does not supply") says it does not. **I retract "open fixed point of deck Z2."** (The A4 ⟹ Pin⁻ ⟹ NS statement itself, which the first report also gave, is correct per `:563-570`.)

---

## B. Strongest still-valid alternative from actual P1

The genuine theorem is untouched: **marked generator theorem** (`SIMPLIFICATION.txt:16-108`): any connected local unitary same-conformal extension `W⊇U=V_D5⊗V_A3` with `Hom_U(M_h,W)≠0` is forced to `V_E8 = M_0+M_h+M_2h+M_3h`, each once, via `n_a∈{0,1}`, isotropy of `q(a,b)=(5a²+3b²)/8`, and `H⊆S⊆H^⊥=H`. The four-class classification `{U, D8, E8_+, E8_−}` (`:117-147`) is exact and `h=(1,1)` selects `E8_+`.

What survives as representation theory: `M_h=(16,4)` **is** the `so(10)⊕so(6)`-chirality component of the Spin(16) spinor `128_+`, weight 1, dim 64 (`NEUTRAL-NET-AND-DEFECT.txt:242`, dictionary `(12)` at `ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:204-211`). This identification is valid.

What is **not** valid is calling its realization a derived source selection. Realizing one *local* `M_h` requires **choosing a spin/Ramond extension net** (the fermionic spin net plus the Arf/defect lift `g_hat=(+,−,−,+)`, `AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:211-238`) and then proving an **energy-bounded localized charged intertwiner** (eq (32), `NEUTRAL-NET-AND-DEFECT.txt:382-392`). Per the request, a chosen extension is **not** a derived P1 selection. So the strongest honest statement is:

> **Theorem (B) + open premise (4.8):** given a local same-conformal realization containing `M_h`, `V_E8` and its vacuum follow. Supplying that `M_h` needs a chosen Ramond extension and a defect-localization theorem — an added input, not a consequence of P1.

The global operator `W_ext` (`SAME-SOURCE-AFFINE-DEFECT.txt`, EXISTENCE PROVED) is the correct *global representation intertwiner*; its boundary is exactly locality (`:327-333`), and it targets `v+c`, not the Ramond `M_h`. It is evidence that the algebra is rich, not that `M_h` is already a local field.

---

## C. Independent assessment of `CHIRAL-SEWING-SELECTION.txt`

**The finite identity is correct.** With `U=V_{D5}⊗V_{A3}`, `A=C4×C4`, four isotropic supports `H0,H2,H_+,H_−` (matching `SIMPLIFICATION.txt:117-121`), `S_{x,y}=(1/4)exp(−2πi b(x,y))`, and `v_H=Σ_{h∈H}e_h`, the x-coordinate of `Sv_H` is `(1/4)Σ_{h∈H}exp(−2πi b(x,h)) = |H|/4` for `x∈H^⊥`, else 0. Hence `Sv_H=(|H|/4)v_{H^⊥}`. Demanding a single S-stable chiral vacuum line (`Sv_H=λv_H`, `|λ|=1`, vacuum coord `λ=|H|/4>0`) forces `|H|=4`, `λ=1`, `H=H^⊥`, selecting `H_+` or `H_−` (the two E8). The table (§2 lines 63-67) is right. Scalar/paired modular invariance (`Z_pair=χ^*Iχ`, rank-16 pairing `v_Hv_H^*` vs rank-1) does **not** select it (§3, lines 84-99) — correct, and the diagonal full CFT is a valid counter to "scalar invariant ⟹ single chiral block." The caveats against `tau=i` tautology and CAR-purity (§3 lines 100-107) are sound. **No invented theorem from dimensions/purity/`tau=i` is used.**

**Does ORIGINAL P1 furnish the line and the cycle exchange? No — it is a new assumption.** The original spectral-bordism objects are `B=(A_+,H_+,D_+,J,Γ)` with morphisms "**unitary intertwiners preserving the displayed one-sided data**" (`002-section…tex:252-267`). The primitive completion (`027-subsection…tex:23-66`) outputs a *discrete* datum (`S^+=Λ^{even}(E_3⊕E_2)`, `dim 16`, `X_f≅P^1\μ_4`, square `D_4`, `N_fam=3`, `(g,b,s)=(5,3,2)`). **Neither passage defines a genus-one conformal-block functor on the 16 U-modules, a chiral vacuum branching vector `v_H`, nor a cycle-exchange operator.** This is exactly `CHIRAL-SEWING §4` lines 109-122, and the missing "algebra-fiber classification" is the same gap flagged in `ORIGINAL-SOURCE-AUDIT.txt:60-74, 162-175`.

**Overlooked candidate mechanism — honest verdict.** The nearest structure in actual P1 is the `D_4` symmetry of the square `μ_4={1,−1,i,−i}` (tex Step 2, lines 92-102), whose cyclic part is the quarter rotation `z→iz` (= the A3 quarter; `WEAKEST-DIRECT:106-116` `c`,`Q` with `QcQ^{-1}=c^{-1}`). A "two-cycle exchange" is plausibly the reflection class of this `D_4`. **But** identifying `D_4` on `P^1\μ_4` with the modular `S` cycle-exchange on `C[C4×C4]` is precisely the undefined block functor. I find **no derived mechanism** in the two original passages; the square `D_4` is suggestive, not a furnished `v_H`+exchange. So CHIRAL-SEWING's selection is a **new assumption**, internally correct, equivalent in strength to selecting holomorphy — not weaker, and not supplied by P1.

---

## D. Smallest computation that decides local-charged-field vs chiral-sewing

Both routes reduce to the *same* missing bridge: a map from the P1 source datum to the chiral `U`-theory. They differ in the object that bridge must produce:

- **Local charged-field route** needs a bounded/affiliated charged generator — concretely **eq (32)** (`NEUTRAL-NET-AND-DEFECT.txt:382-392`): `X J_h(f) X* = J_{Rh}(f)+shift`, `X*HX=H+Q+a`, with `X` a localized net defect. Its global version already holds (`W_ext`); the open content is localization/boundedness on the source net.
- **Chiral-sewing route** needs a **single chiral vacuum block `v_H`** with cycle-exchange invariance — `CHIRAL-SEWING` eq (1).

**The single smallest existing source object that discriminates** is the unresolved A3 role, `ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:647-670`: whether `D_Γ=i∂_θ`, `λ_n(φ)=n+φ/2π` is **(reading 1)** a scalar/determinant matching operator or **(reading 2)** the physical tangential spinor Dirac operator. Its resolution is "the smallest missing local dictionary" (`:668-670`) and it decides the fork:
- Reading 2 (physical spinor, with its Ramond zero mode) is the only thing that could seed a *local charged* Ramond field — route A.
- Reading 1 (scalar determinant) leaves only global determinant/sewing data on the seam — pushing toward the genus-one block/cycle route C.

**The single smallest new computation** that would settle it either way is to map the raw RP/P1 seam functional into the 16 U-module genus-one block space and compute its **branching/pairing rank**: a factorized rank-one `v_H` (with cycle-exchange invariance) closes the chiral-sewing route (selecting `E8_+`); a rank-16 pairing (as OS/Schlingemann reconstruction generically yields a *full*, non-chiral theory, `ORIGINAL-SOURCE-AUDIT.txt:124-161`) rules out a single chiral line and forces the local-charged-field route via eq (32). That one rank computation — once the block map exists — is the minimal decider; the A3 scalar-vs-spinor passage is the minimal *already-printed* source datum pointing to which map is even available.

---

### Net retractions
`σ=W_ext^{1/2}`; `v+c` = Ramond; bounded `σ∈F(I)`; "only localization remains"; `Γ` without `2^8` / `SO(16)` spinor; "1-dim Ramond ground"; blanket categorical refutation; "chiral subnet discards half D8"; A3 unit flow = NS/R interchange; "open fixed point of deck Z2." The marked-generator theorem, the representation-level identity of `M_h`, and the finite CHIRAL-SEWING identity stand; the physical realization of one local `M_h` and the cycle-exchange premise remain open added inputs, not P1-derived selections.