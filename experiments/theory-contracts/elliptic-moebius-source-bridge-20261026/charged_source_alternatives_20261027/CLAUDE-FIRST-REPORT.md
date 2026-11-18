I now have a complete, well-grounded picture. Here is the self-contained report.

---

# Alternative route to the missing local charged generator M_h = (16,4)

**Scope note.** This is original-source mathematical analysis over the files listed in the request plus the artifacts they decisively cite. I did not scan the whole repository and make no novelty claim against an exhaustive search. Atlas `program.json` is the sole current Atlas; I treat the dated gesamtkontext as orientation only.

---

## 0. The problem, restated precisely

Two arrows must not be conflated (established in `WEAKEST-DIRECT-E8-SOURCE-ROUTE.txt:48-72`):

- **Arrow B — the marked‑generator theorem (THEOREM).** With `U = V_D5 ⊗ V_A3`, `A = C4×C4`, `h=(1,1)`: any *connected local unitary same‑conformal* VOA extension `W ⊇ U` with `Hom_U(M_h,W) ≠ 0` is forced to be `V_E8`, with `W = M_0+M_h+M_2h+M_3h`, each sector once. Proof via regularity of `U`, the free‑module adjunction `n_a = dim Hom_B(M_a⊗B,B) ∈ {0,1}`, isotropy under `q(a,b)=(5a²+3b²)/8`, and `H ⊆ S ⊆ H⊥ = H` (`SIMPLIFICATION.txt:16-108`; Kirillov–Ostrik, DMNO, Huang–Kirillov–Lepowsky, Carpi–Gaudio–Giorgetti–Hillier). **This is a genuine theorem, conditional only on realizing one `M_h`.**

- **Arrow C — physical origin (OPEN, SEAM.EQUIV.01).** The raw RP/P1 seam state *is* that `V_E8` vacuum, and the seam actually realizes one local `M_h`. This is the active open parent (`WEAKEST-DIRECT:41-45, 394-432`).

So the whole question collapses to a single, sharply typed object: **where does one local charged `M_h = (16,4)` come from on the physical source side?** Everything downstream is then a theorem. The request's demand — "do not present a conditional E8 model as a physical origin theorem" — is exactly the Arrow‑B/Arrow‑C distinction, which I preserve throughout.

**What `M_h=(16,4)` is, in the lattice.** In `248 = (45,1)+(1,15)+(16,4)+(10,6)+(bar16,bar4) = 60+64+60+64` (`WEAKEST-DIRECT:277-280`), the glue degrees 0,2 are the 112 `D8` roots (`M_0 ⊕ M_2h = D8` currents, 120 at weight 1) and degrees 1,3 are the 128 `D8` spinor roots (`M_h ⊕ M_3h`). Thus
> **`M_h=(16,4)` is one chirality component of the 128‑dimensional `D8=SO(16)` spinor coset `s` of the bosonic extension `E8 = D8 + s`** — lowest weight `h_s = 1`, lowest‑weight space dim 64 (`ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:209`; `NEUTRAL-NET-AND-DEFECT.txt:242`; `WEAKEST-DIRECT:349-363`).

**Why ordinary fermionization cannot produce it (VERIFIED exclusion).** The existing NS vector route yields `(10,1)+(1,6)` only; the central kernel element of `Spin(10)→SO(10)` acts as `+1` on every local NS generator but `−1` on the internal `16`. All products/derivatives/OPE/weak closures stay centrally trivial (`PRIOR-SCOPED-RESULTS.json`, finding `unchanged-ns-vector-generators-…`; sourced to `seam-collar-multifield-20261016/P1-COVARIANT-INSERTION.txt:147-213`). **This is a Z2‑central (spinoriality) obstruction.** Its scope note explicitly reads: *"No exclusion of Ramond/twisted states, charged sector intertwiners."* That is the designated open door.

---

## 1. Ranking of the three routes

| Rank | Route | Verdict |
|---|---|---|
| **1 (best)** | **Local Ramond spin‑field endpoint** of the 16 Majorana‑Weyl fermions; `M_h` = chirality component of their joint Ramond ground state. | **Survives every scoped exclusion; already half‑built (EXISTENCE PROVED machinery). Keeps the full vector source. Develop this.** |
| 2 | Switch the complete vector source to a **chiral current subnet**, where charged spinor pairs are not Ward‑forbidden. | **Not justified by original P1** — a convenient weakening; breaks `D8`‑locality. See §4. |
| 3 | **Categorical boundary completeness / CAR dual line `K*`** as the charged generator. | **Type‑obstructed** — supplies a 3/4‑dimensional determinant line, never a 64‑dim module. See §5. |

The request's hypothesis ("a local Ramond endpoint/sector intertwiner, *or* an already present categorical boundary completeness principle") is therefore **confirmed in its first disjunct and refuted in its second by a proven type obstruction.**

---

## 2. The best route (R1): M_h is a localized Ramond spin field

### 2.1 Exact identification (VERIFIED SOURCE FACT at the representation level)

A4 read literally (Pin⁻ lift of the Möbius `Z2`) gives `r² = −I`, i.e. **NS** on the orientable‑double core (`ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:541-575`). In the marked frame, the glue‑square character
```
chi_can(E_r) = (-1)^deg(r) = G²  :  +1 on the 112 D8 roots,  -1 on the 128 spinor roots
```
is exactly the non‑anomalous `E8→D8` `Z2` of standard fermionization (`:119-168`, eqs (8)–(10)). Fermionizing produces **16 Majorana–Weyl fermions** `ψ^a (a=1..16)` with the four‑module dictionary (`NEUTRAL-NET-AND-DEFECT.txt:214-222`):
```
H_F,NS = H_0 ⊕ H_v ,     H_F,R = H_s ⊕ H_c ,     E8_bosonic = H_0 ⊕ H_s .
```
The currents `J^{ab} = :ψ^aψ^b:` generate `so(16)_1 = (D8)_1`; the marked `D5+A3 = so(10)⊕so(6) ⊂ so(16)` splits `a` into `1..10` and `11..16`.

Now the decisive numerics, all independently fixed by the source and all matching:

- **Ramond ground state = spin field σ.** `σ` carries the zero‑mode Clifford algebra `{ψ^a_0, ψ^b_0} = δ^{ab}` on `2^8 = 256` states; the chirality `Γ = i^8 ψ^1_0…ψ^16_0` (= Ramond fermion parity) splits `256 = 128_+ ⊕ 128_-`.
- **Branching** `SO(16) → SO(10)×SO(6)`: `Spin(SO(16)) = Spin(SO(10))⊗Spin(SO(6)) = (16⊕\bar16)⊗(4⊕\bar4)`, so
 ```
 128_+ = (16,4) ⊕ (bar16, bar4) = M_h ⊕ M_3h .
 ```
 This reproduces the degree‑1/3 branching `64+64` of (4.7) exactly.
- **Conformal weight** of the 16‑fermion spin field: `16 × (1/16) = 1 = h_s = h_c` — matching `ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:209` and `NEUTRAL-NET-AND-DEFECT.txt:242` on the nose.
- **Lowest‑weight dimension** of `(16,4)` = `16·4 = 64` — matching `WEAKEST-DIRECT:349-351`.
- **Central obstruction resolved automatically:** `σ` is a genuine `Spin(16)` (double‑cover) object; the `Spin(10)` central element is `−1` on the `16` *because* `σ` has `−1` monodromy with each `ψ^a`. This is precisely why no single‑valued NS vector product can reach it (finding 5) — the escape is sectorial, not algebraic.

> **Conclusion (representation level):** `M_h = (16,4)` is, uniquely and without fitting, the `so(10)⊕so(6)`‑chirality component of the joint **Ramond ground state (spin/order field) of the 16 Majorana–Weyl fermions** produced by the A4 fermionization.

### 2.2 How much is already built (EXISTENCE PROVED)

This is not a fresh posit. The source stack already contains every algebraic ingredient except localization:

- **Global sector‑changing operator `W_ext`** (`PROOF.txt:157-176`, eq (8)): `W_ext|n+c0 s⟩ = η_R(n) κ^{c0} |Rn+e0+c0 s⟩`, `κ=±i`, which maps `0↔v, s↔c`, restricts to the original NS winding `W1`, and sends the `E8` representation to its twisted module; its half‑integral powers land in the `v+c` (Ramond) representation (`WEAKEST-DIRECT:419`).
- **Li Δ‑twist functor** realizing the spinor vertex as a `g`‑twisted field (`AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:141-199`, eqs (8)–(13)), with `V_{E8+e0}=H_v⊕H_c` and the correct fractional‑power monodromy.
- **Four‑sector abelian intertwining algebra** with **one‑dimensional Ramond ground intertwiners**; associativity + braiding fix all products up to a 2‑cochain (`AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:240-265, 283-288`). Status recorded there: **EXISTENCE PROVED**; only the Arf/defect lift phase `(16)` vs `(14)` is a convention choice.

What is *not* yet built, verbatim from the sources (`AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:184-199, 290-294`; `NEUTRAL-NET-AND-DEFECT.txt:371-392`; `WEAKEST-DIRECT:420-432`): a **bounded, point‑local, `F(I)`‑affiliated** realization of this intertwiner on the physical source net — i.e. energy bounds + a defect‑localization theorem, plus selection from raw P1.

---

## 3. (b) First newly assumed input, and (c) smallest decisive check

### (b) The single new input
Everything in §2.1–2.2 — dimensions, weights, charges, monodromy signs, the global intertwiner, the twisted functor — already exists. The **one genuinely new premise** is:

> **(N)** The primitive P1 seam supplies a **bounded, seam‑interval‑localized Ramond spin field** `σ ∈ F(I)` (equivalently: `W_ext^{1/2}` promoted from a global representation intertwiner to a point‑local smeared field obeying energy bounds), with the `(16,4)` covariance, inside the *same* conformal net.

This is strictly smaller and better‑typed than the previously demanded "full P1→C3 bulk Dirac equivalence" (`WEAKEST-DIRECT:394-432`) and is exactly the active content of SEAM.EQUIV.01.

### (c) Smallest decisive check
Because weights/dimensions/charges/global intertwiner are already verified, the decisive residual check is **localization + the mixed NS/R zero‑mode cocycle**, in three sub‑items of which only the last is new:

1. **Ramond defining relation** (fixes the sector): `ψ^a(e^{2πi}z)\,σ(0) = -\,ψ^a(z)\,σ(0)` for all 16 `a` — half‑monodromy with the fermions, integer monodromy with the currents `J^{ab}` (so `σ` is **local w.r.t. the full `D8`**, since `s·(±e_i±e_j) ∈ Z`).
2. **Clifford/branching**: the 16 zero modes act on `σ`'s 256‑space as `Cl(16)`, restricting under `so(10)⊕so(6)` to `(16,4)⊕(bar16,bar4)`; `h(σ)=1`. *(Already fixed.)*
3. **NEW decisive item — boundedness/localization**: the twisted module `H_s` admits a **conformal‑net energy bound making `σ` a localized bounded defect endpoint affiliated to `F(I)`** (`AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:189-199, 294` names exactly this as "the smallest decisive remaining test"). Concretely, verify the mixed NS/R cocycle matrix of `NEUTRAL-NET-AND-DEFECT.txt:384-392`, eq (32):
 ```
 X J_h(f) X* = J_{Rh}(f) + (source zero-mode shift),   X* H X = H + Q + a,   with X ∈ F(I) bounded.
 ```
The global version of (32) already holds (that is `W_ext`); **the one‑line new content is `X ∈ F(I)` + energy bound.**

---

## 4. (d) Why R1 survives every exact older exclusion

R1 is categorically different from every excluded mechanism, and each exclusion is scoped to *exclude a different object*:

- **vs. NS vector products** (finding 5): scoped to the *unchanged NS generator algebra*; explicitly leaves "Ramond/twisted states, charged sector intertwiners" open. `σ` is precisely such an object. ✔
- **vs. uniform scalar `Z4` twist** (finding `uniform-z4-d8-locality`; `FULL-FINITE-ENERGY:142-186`): `σ` is not a uniform lattice shift. The `r2` half‑class there is noted to be "the only... D8‑local" uniform candidate, but "the actual local endpoints, chirality selection/Gauss law and common time have not been produced" — `σ` *is* exactly that missing local endpoint, not a competing uniform shift. ✔
- **vs. cover half‑flux** (`excluded_scoped`): a connection holonomy in the fixed charge dictionary, `C(-h)^8=-hI ≠ hI`. `σ` is an internal‑sector twist of the *same* currents, different domain entirely. ✔
- **vs. complete Dirac vector phases** (`excluded_scoped`, zero spinor Ward charge): those are `V_vec(f)=exp(iJ_vec(f))` histories; `σ` is an order/disorder field, *not* of exponential‑current form, so the vanishing‑Ward computation (`NEUTRAL-DEFECT-PAIRS:76-120`) does not apply. ✔
- **vs. continuous‑time nonunitary limits** (`excluded_scoped`): those stay in the old folium of the *unchanged* source map; `σ` generates a genuinely disjoint sector, and that finding explicitly states "disjoint‑state limits do not automatically violate B" and claims "no classification of all local scaling nets." ✔
- **vs. the smooth‑unitary energy divergence** (`NEUTRAL-DEFECT-PAIRS:207-268`, eq (12), `ΔK ≥ 2π/ε → ∞`): that excludes *smooth unitary phase separation at fixed global energy*; `σ` is a charged intertwiner in a new positive‑energy sector, which §9–10 there explicitly leave open ("does not exclude the finite positive charged generators of the standard sector representations"). The Ramond ground energy is finite and fixed (`h_s=1`). ✔

**No contradiction with A4.** A4 fixes **NS on the closed double‑core cycle**; `σ` is a **Ramond twist at the seam puncture** (an open fixed point of the deck `Z2`). NS vacuum + a localized spin‑field insertion is the standard, mutually consistent CFT picture. The sources confirm A4 does *not* fix the global spin‑surface data or forbid the twisted sector (`ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:740-752`), and the clock implementer already satisfies `U² = (-1)^F` (its twisted sector *is* the `σ`‑sector) (`:344-353`).

---

## 5. Why routes 2 and 3 are inferior — and the direct answer to the vector→chiral/categorical test

**R3 (categorical/dual line) is type‑obstructed.** The CAR dual line `K* = det(V)^{-1}` already exists by adjunction with `β|K* = +i` (`P1-ETA-SEWING-DUAL-LINE-TEST.txt:52-101`), and the theta/Heisenberg functor gives a **4‑dimensional** `W_aux` with `γ³=I, α²=iI`. But `M_h` has a **64‑dimensional** lowest‑weight space and infinitely many oscillator descendants; the theta functor supplies "neither the Heisenberg Fock space, nor the twisted group algebra and cocycle, nor vertex operators, OPE locality, positive‑energy `L0`, nor a vacuum functional" (`WEAKEST-DIRECT:343-367`). **This is a proven type obstruction, not a missing numerical check.** Hence an "already present categorical boundary completeness principle" *cannot* be the origin of `M_h`. The dual line is the right home for the 3‑dimensional charged *defect* `V`, not for the generator.

**R2 (chiral subnet) is a convenient replacement, not justified by original P1.** The only escape the vector‑Ward exclusion leaves is "charged spinor fields relative to a **chiral** current subnet" (`NEUTRAL-DEFECT-PAIRS:118-120`). But switching to a chiral subnet means *discarding half the `D8` currents*; the quarter class `r1` is then precisely the one that is **not** `D8`‑local (`FULL-FINITE-ENERGY:146-162`). Original P1 supplies the *full* vector/`D8` structure (the marked `D5+A3 ⊂ D8` glue, the complete Dirac/CAR source), so truncating to a chiral subnet is a weakening introduced *for convenience*, not something the seam data forces.

> **Direct answer to the request's test.** The vector→chiral/categorical switch is **not** justified by original P1. The correct move is not a change of *source type* but a change of *sector*: go to the **Ramond sector of the same full vector source**. The spin field `σ` is local with respect to **all** `D8` currents (`s·D8‑root ∈ Z`), so the complete vector source is retained intact. R1 is the only one of the three that keeps the source and still resolves the central‑kernel obstruction.

---

## 6. An original‑source anchor, and its honest status (CONJECTURE)

There is a real anchor for R1 inside the original seam data, which I flag carefully because the sources did **not** adopt it. A3 prints `D_Γ = i∂_θ` with `λ_n(φ) = n + φ/2π`; **at the reference `φ=0` this is the periodic Ramond spectrum with a zero mode, `dim ker = 1, ξ = 1/2`** (`ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:648-680`). The spectral flow `φ: 0→2π` has exactly one crossing — the toggling of a single zero mode, i.e. the NS↔R interchange a spin‑field insertion performs.

The sources resolve the A3/A4 tension by **reading‑1** (A3 is a scalar matching/determinant operator, so NS is preserved). **Reading‑2** (A3 *is* the physical tangential spinor Dirac operator, and its Ramond zero mode is physical) is the natural home of `σ`: the seam endpoint carries the Ramond twist whose zero mode, in the fermionized `SO(16)` frame, is the `256`‑state Clifford module whose `128_+` chirality is `(16,4)⊕(bar16,bar4)`. Reading‑2 is **conjecture**: the sources lean to reading‑1, and both only agree that the spectral‑flow magnitude is 1. I do not claim reading‑2; I record it as the sharpest available original‑source motivation that (N) is *not* an arbitrary posit but the physical content of the already‑printed A3 zero mode.

---

## 7. Status taxonomy (explicit)

- **THEOREM.** Arrow B: one same‑conformal local `M_h` ⟹ `V_E8` with all structure (`SIMPLIFICATION.txt`).
- **VERIFIED SOURCE FACTS.** A4⟹NS on the double core; `chi_can=G²` = non‑anomalous `E8→D8` `Z2`; four‑module dictionary and `H_s⊕H_c` = Ramond; `h_s=h_c=1`; the `SO(16)` spinor `128_+ = (16,4)⊕(bar16,bar4)` under `so(10)⊕so(6)`; central‑kernel exclusion of NS vector products; **EXISTENCE PROVED** for `W_ext`, the Li Δ‑twist functor, and the four‑sector abelian intertwining algebra with 1‑dim Ramond ground intertwiners.
- **CONDITIONAL CONSTRUCTION (R1).** `M_h` realized as a localized Ramond spin field `σ`, assuming (N): `σ` is a bounded `F(I)`‑affiliated field with energy bounds. All non‑localization data are already verified; only boundedness/localization and the Arf/defect lift phase remain.
- **CONJECTURE.** SEAM.EQUIV.01 (raw RP seam state = `V_E8` vacuum); and reading‑2 of A3 as the physical spinor operator supplying the Ramond zero mode.
- **REFUTED as origin of `M_h`.** Categorical/theta/dual‑line route (type obstruction); chiral‑subnet switch (unjustified weakening).

---

## 8. Sources and line ranges

- Generator theorem: `e-series-source-reconstruction-20261024/SIMPLIFICATION.txt:16-169`.
- Three‑arrow separation, weak field theorem, type obstruction, `W_ext` boundary: `elliptic-moebius-source-bridge-20261026/theta_lattice_closure_20261027/direct/WEAKEST-DIRECT-E8-SOURCE-ROUTE.txt:48-72, 236-297, 343-367, 409-432`.
- `chi_can=G²`, A4⟹NS, A3 Ramond zero mode, operator dictionary: `canonical-seam-fermionization-20261025/ORIGINAL-SPIN-AND-SEAM-DICTIONARY.txt:119-168, 541-575, 648-680, 740-752`.
- Four‑module dictionary, neutral net, mixed NS/R cocycle test (eq 32): `canonical-seam-fermionization-20261025/NEUTRAL-NET-AND-DEFECT.txt:214-222, 242, 371-392`.
- Affine defect: existence proved, Li Δ‑twist, Ramond ground intertwiners, localization as the decisive open test: `canonical-seam-fermionization-20261025/AFFINE-DEFECT-INDEPENDENT-REVIEW.txt:141-199, 240-294`; `.../PROOF.txt:139-176`.
- Sector menu, `D8`‑locality of the half/spinor class, chiral‑subnet caveat: `mixed-current-origin-20261014/FULL-FINITE-ENERGY-AND-SECTOR-MENU.txt:142-186, 208-230`; `.../NEUTRAL-DEFECT-PAIRS-AND-VECTOR-WARDS.txt:76-120, 207-268`.
- Dual line `K*` type/role, categorical refutation: `elliptic-moebius-source-bridge-20261026/theta_lattice_closure_20261027/eta/P1-ETA-SEWING-DUAL-LINE-TEST.txt:52-101, 269-345`.
- Direct seam defect, relative determinant `1−z`: `original-admissible-source-closure-20261025/DIRECT-DEFECT-CLOSURE.txt:20-132`.
- Nearest scoped exclusions and their explicit openings: `PRIOR-SCOPED-RESULTS.json` (findings 1–7, especially `unchanged-ns-vector-…` scope).

**Bottom line.** The strongest alternative route is not a change of source type but a change of sector: `M_h=(16,4)` is the `so(10)⊕so(6)`‑chirality component of the **Ramond spin field of the 16 Majorana–Weyl fermions**, weight 1, dimension 64 — the exact object left open by the central‑kernel exclusion, local with respect to the full vector `D8` currents, and already built globally as `W_ext`/the Δ‑twisted intertwiner. The remaining physical gap reduces to one well‑typed check — **bounded `F(I)`‑localization with energy bounds** — which is strictly smaller than the full P1→C3 Dirac equivalence and is the honest residue of SEAM.EQUIV.01. The categorical completeness principle is refuted as the origin by a proven type obstruction; the chiral‑subnet switch is an unjustified weakening of the original vector source.