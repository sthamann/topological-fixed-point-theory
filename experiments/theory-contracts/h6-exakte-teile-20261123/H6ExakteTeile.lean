import Mathlib

/-!
# H6: exact, finite parts of the exploration of 23 November 2026

Source: `experiments/exploration/2026-11-23-phasenraum-f5-holocode/VERDICT.md`
(H6h3 to H6h8) and the scripts `z1_kritik_pruefung.py`, `z5_gabel.py`,
`z8_sechzehn_und_saetze.py`, `z10_kritik_zwei.py`, `z3_hintergrundfreiheit.py`.

This file proves only finite algebra and arithmetic:

* A. the lifted quarter turn on a momentum orbit (z1, point 2): `g⁴ = -1`,
  `tr g = tr g² = tr g³ = 0`, `[g, D] = 0`, characteristic polynomial
  `(X⁴ + 1)²`, every one of the four eighth-root eigenspaces has dimension
  exactly `2 = 8 / 4`, and each of them meets each eigenspace of `D` in exactly
  one dimension (a quarter of it), on every orbit; no invariant sector; orbits
  of momenta in `(ℤ + 1/2)²` have length exactly 4;
* B. the matrices of `P(z) ↦ i P(i z)` on `{1, z, z²}` and of `Q(z) ↦ Q(i z)`
  on `{1, …, z⁴}` (z10, point 3), the half-turn split `3 + 2`, and the
  dimension count `1 + 10 + 5 = 16`;
* C. arithmetic of the central charge, the twist weight and the GSO phases
  (z1 point 4, z5, z8, z10 point 4);
* D. the charge `X = 5(B-L) - 4Y` of one Standard Model generation (z8);
* E. every level-one weight-6 modular form (hence `E₆`) vanishes at `τ = i`,
  and the arithmetic of the valence formula.  The valence formula itself is
  **not** proved here; it enters only as an explicit hypothesis.

Nothing here proves a physical identification (seam, families, E₈ phase,
coupling of physical fermions, state, dynamics).  No axiom is declared.
-/

set_option autoImplicit false

namespace TFPT.H6Exact

open Matrix Polynomial
open scoped Kronecker

/-! ## A. Lifted quarter turn on a momentum orbit -/

section QuarterTurn

variable {R : Type*} [CommRing R]

/-- Weighted cyclic shift of the four momentum lines: line `j` goes to line `j + 1`
with phase `p (j + 1)`.  Rows are targets, columns sources (as in `z1`, `block`). -/
def shift4 (p : Fin 4 → R) : Matrix (Fin 4) (Fin 4) R :=
  !![0, 0, 0, p 0;
     p 1, 0, 0, 0;
     0, p 2, 0, 0;
     0, 0, p 3, 0]

/-- Spinor rotation by `π/2`, `diag(ζ⁻¹, ζ)`, written with `ζ⁻¹ = -ζ³` (valid when `ζ⁴ = -1`). -/
def spinRot (ζ : R) : Matrix (Fin 2) (Fin 2) R := diagonal ![-ζ ^ 3, ζ]

/-- Lifted quarter turn on one orbit with translation phases `p`:
8×8, index `(momentum line, spinor component)`. -/
def gOrb (ζ : R) (p : Fin 4 → R) : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) R :=
  shift4 p ⊗ₖ spinRot ζ

/-- Translation phases `exp(-2πi k_j·a)` for `k₀ = (1/2, 3/2)`, `a = (1/2, 0)`:
`(-i, -i, i, i)` with `i = ζ²` (see `phases_zeta8` for the identification in `ℂ`). -/
def phases (ζ : R) : Fin 4 → R := ![-ζ ^ 2, -ζ ^ 2, ζ ^ 2, ζ ^ 2]

/-- The operator `g̃` of `z1` (orbit of `k₀ = (1/2, 3/2)`, rotation centre `(1+i)/4`). -/
def gQ (ζ : R) : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) R := gOrb ζ (phases ζ)

theorem kronecker_pow {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (A : Matrix m m R) (B : Matrix n n R) (k : ℕ) :
    (A ⊗ₖ B) ^ k = (A ^ k) ⊗ₖ (B ^ k) := by
  induction k with
  | zero => simp
  | succ k ih => rw [pow_succ, ih, ← mul_kronecker_mul, pow_succ, pow_succ]

theorem shift4_pow_four (p : Fin 4 → R) :
    shift4 p ^ 4 = (p 0 * p 1 * p 2 * p 3) • (1 : Matrix (Fin 4) (Fin 4) R) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [shift4, pow_succ, Matrix.mul_apply, Fin.sum_univ_four] <;> ring

theorem trace_shift4_pow_one (p : Fin 4 → R) : trace (shift4 p ^ 1) = 0 := by
  simp [shift4, trace, Fin.sum_univ_four]

theorem trace_shift4_pow_two (p : Fin 4 → R) : trace (shift4 p ^ 2) = 0 := by
  simp [shift4, trace, Fin.sum_univ_four, pow_succ]

theorem trace_shift4_pow_three (p : Fin 4 → R) : trace (shift4 p ^ 3) = 0 := by
  simp [shift4, trace, Fin.sum_univ_four, pow_succ]

variable (ζ : R)

/-- **`g̃⁴ = -1`** on every orbit whose phases multiply to one. -/
theorem gOrb_pow_four (hζ : ζ ^ 4 = -1) (p : Fin 4 → R) (hp : p 0 * p 1 * p 2 * p 3 = 1) :
    gOrb ζ p ^ 4 = -1 := by
  have hS : spinRot ζ ^ 4 = -1 := by
    rw [spinRot, diagonal_pow]
    ext i j
    fin_cases i <;> fin_cases j <;> simp
    · linear_combination (ζ ^ 8 - ζ ^ 4 + 1) * hζ
    · exact hζ
  rw [gOrb, kronecker_pow, shift4_pow_four, hp, one_smul, hS]
  ext ⟨a, b⟩ ⟨c, d⟩
  simp [one_apply]
  split_ifs <;> simp_all

theorem phases_prod (hζ : ζ ^ 4 = -1) :
    phases ζ 0 * phases ζ 1 * phases ζ 2 * phases ζ 3 = 1 := by
  simp [phases]
  linear_combination (ζ ^ 4 - 1) * hζ

/-- **`g̃⁴ = -1`** (four quarter turns give `-1` on spinors). -/
theorem gQ_pow_four (hζ : ζ ^ 4 = -1) : gQ ζ ^ 4 = -1 :=
  gOrb_pow_four ζ hζ (phases ζ) (phases_prod ζ hζ)

/-- `tr g̃ᵏ = 0` for `k = 1, 2, 3`, on every orbit. -/
theorem trace_gOrb_pow (p : Fin 4 → R) (k : ℕ) (hk1 : 1 ≤ k) (hk3 : k ≤ 3) :
    trace (gOrb ζ p ^ k) = 0 := by
  rw [gOrb, kronecker_pow, trace_kronecker]
  interval_cases k
  · rw [trace_shift4_pow_one, zero_mul]
  · rw [trace_shift4_pow_two, zero_mul]
  · rw [trace_shift4_pow_three, zero_mul]

theorem trace_gQ_pow (k : ℕ) (hk1 : 1 ≤ k) (hk3 : k ≤ 3) : trace (gQ ζ ^ k) = 0 :=
  trace_gOrb_pow ζ (phases ζ) k hk1 hk3

/-- Dirac symbol `kx σx + ky σy` with `i = ζ²`. -/
def diracSym (kx ky : R) : Matrix (Fin 2) (Fin 2) R :=
  !![0, kx - ζ ^ 2 * ky; kx + ζ ^ 2 * ky, 0]

/-- The spinor rotation turns the Dirac symbol of `k` into that of `rot k = (-ky, kx)`. -/
theorem spinRot_intertwines (hζ : ζ ^ 4 = -1) (kx ky : R) :
    spinRot ζ * diracSym ζ kx ky = diracSym ζ (-ky) kx * spinRot ζ := by
  have h5 : ζ ^ 5 = -ζ := by linear_combination ζ * hζ
  ext i j
  fin_cases i <;> fin_cases j <;> simp [spinRot, diracSym, Matrix.mul_apply, Fin.sum_univ_two] <;>
    ring_nf <;> simp only [h5] <;> ring

omit ζ in
/-- x-components of the orbit of `(u, v)` under `rot (x, y) = (-y, x)`. -/
def orbX (u v : R) : Fin 4 → R := ![u, -v, -u, v]

omit ζ in
/-- y-components of the orbit of `(u, v)`. -/
def orbY (u v : R) : Fin 4 → R := ![v, u, -v, -u]

omit ζ in
/-- The four lines form one orbit of the quarter turn `rot (x, y) = (-y, x)`. -/
theorem orbit_is_rotation (u v : R) (j : Fin 4) :
    orbX u v (j + 1) = -orbY u v j ∧ orbY u v (j + 1) = orbX u v j := by
  fin_cases j <;> simp [orbX, orbY]

/-- Doubled momenta `2 k_j` of the orbit of `k₀ = (1/2, 3/2)`: `(1,3), (-3,1), (-1,-3), (3,-1)`. -/
def kx2 : Fin 4 → ℤ := orbX 1 3

/-- See `kx2`. -/
def ky2 : Fin 4 → ℤ := orbY 1 3

/-- Quarter turn on doubled momenta. -/
def rot2 (k : ℤ × ℤ) : ℤ × ℤ := (-k.2, k.1)

/-- No momentum in `(ℤ + 1/2)²` (doubled: both coordinates odd) is fixed by the quarter
or the half turn, while `rot⁴ = id`: every orbit has exactly four elements. -/
theorem rot2_orbit_length_four (k : ℤ × ℤ) (hx : Odd k.1) :
    rot2 (rot2 (rot2 (rot2 k))) = k ∧ rot2 k ≠ k ∧ rot2 (rot2 k) ≠ k := by
  obtain ⟨a, b⟩ := k
  obtain ⟨t, ht⟩ := hx
  simp only at ht
  refine ⟨by simp [rot2], ?_, ?_⟩ <;> simp [rot2, Prod.ext_iff] <;> omega

/-- Dirac operator on the orbit of `(u, v)`, block diagonal in the momentum line. -/
def diracOrb (u v : R) : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) R :=
  Matrix.of fun a b => if a.1 = b.1 then diracSym ζ (orbX u v a.1) (orbY u v a.1) a.2 b.2 else 0

/-- Doubled Dirac operator `2 (k·σ)` on the orbit of `k₀ = (1/2, 3/2)`; `z1` uses `π` times it. -/
def dirac2 : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) R := diracOrb ζ 1 3

set_option maxHeartbeats 2000000 in
/-- **`[g̃, D] = 0`** on the orbit of any `(u, v)`, for any phases. -/
theorem gOrb_comm_diracOrb (hζ : ζ ^ 4 = -1) (p : Fin 4 → R) (u v : R) :
    gOrb ζ p * diracOrb ζ u v = diracOrb ζ u v * gOrb ζ p := by
  have h5 : ζ ^ 5 = -ζ := by linear_combination ζ * hζ
  ext ⟨a, σ⟩ ⟨b, τ⟩
  fin_cases a <;> fin_cases b <;> fin_cases σ <;> fin_cases τ <;>
    simp [gOrb, diracOrb, shift4, spinRot, diracSym, orbX, orbY, Matrix.mul_apply,
      Fintype.sum_prod_type, Fin.sum_univ_two] <;>
    ring_nf <;> (try simp only [h5]) <;> ring_nf

theorem gQ_comm_dirac2 (hζ : ζ ^ 4 = -1) : gQ ζ * dirac2 ζ = dirac2 ζ * gQ ζ :=
  gOrb_comm_diracOrb ζ hζ (phases ζ) 1 3

/-- `D² = |2k|² = u² + v²` on the orbit. -/
theorem diracOrb_sq (hζ : ζ ^ 4 = -1) (u v : R) :
    diracOrb ζ u v ^ 2 = (u ^ 2 + v ^ 2) • (1 : Matrix (Fin 4 × Fin 2) (Fin 4 × Fin 2) R) := by
  ext ⟨a, σ⟩ ⟨b, τ⟩
  fin_cases a <;> fin_cases b <;> fin_cases σ <;> fin_cases τ <;>
    simp [diracOrb, diracSym, orbX, orbY, pow_two, Matrix.mul_apply, Fintype.sum_prod_type,
      Fin.sum_univ_two] <;>
    ring_nf <;> simp only [hζ] <;> ring

/-- `tr (g̃ᵏ D) = 0` for `k = 0, …, 3`. -/
theorem trace_gOrb_pow_mul_diracOrb (p : Fin 4 → R) (u v : R) (k : ℕ) (hk : k ≤ 3) :
    trace (gOrb ζ p ^ k * diracOrb ζ u v) = 0 := by
  rw [gOrb, kronecker_pow]
  interval_cases k <;>
    simp [trace, diracOrb, diracSym, shift4, spinRot, orbX, orbY, pow_succ,
      Matrix.mul_apply, Fintype.sum_prod_type, Fin.sum_univ_two]

/-! ### Characteristic polynomial -/

omit ζ in
theorem charpoly_shift4 (q : Fin 4 → R) :
    (shift4 q).charpoly = X ^ 4 - C (q 0 * q 1 * q 2 * q 3) := by
  rw [charpoly, det_succ_row_zero]
  simp [Fin.sum_univ_succ, det_fin_three, charmatrix_apply, shift4, Fin.succAbove]
  ring

omit ζ in
theorem charpoly_blockDiagonal {m o : Type*} [Fintype m] [DecidableEq m] [Fintype o]
    [DecidableEq o] (M : o → Matrix m m R) :
    (blockDiagonal M).charpoly = ∏ k, (M k).charpoly := by
  have h : charmatrix (blockDiagonal M) = blockDiagonal fun k => charmatrix (M k) := by
    ext ⟨i, k⟩ ⟨j, k'⟩
    simp only [charmatrix_apply, blockDiagonal_apply, diagonal_apply, Prod.mk.injEq]
    by_cases hk : k = k' <;> by_cases hij : i = j <;> simp [hk, hij]
  simp only [charpoly, h, det_blockDiagonal]

omit ζ in
theorem shift4_op_smul (p : Fin 4 → R) (s : R) :
    MulOpposite.op s • shift4 p = shift4 (fun j => p j * s) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [shift4]

/-- **Characteristic polynomial `(X⁴ + 1)²`** on every orbit whose phases multiply to one. -/
theorem charpoly_gOrb (hζ : ζ ^ 4 = -1) (p : Fin 4 → R) (hp : p 0 * p 1 * p 2 * p 3 = 1) :
    (gOrb ζ p).charpoly = (X ^ 4 + 1) ^ 2 := by
  rw [gOrb, spinRot, kronecker_diagonal, charpoly_blockDiagonal, Fin.prod_univ_two]
  simp only [shift4_op_smul, charpoly_shift4]
  have h1 : p 0 * (-ζ ^ 3) * (p 1 * (-ζ ^ 3)) * (p 2 * (-ζ ^ 3)) * (p 3 * (-ζ ^ 3)) = -1 := by
    linear_combination ζ ^ 12 * hp + (ζ ^ 8 - ζ ^ 4 + 1) * hζ
  have h2 : p 0 * ζ * (p 1 * ζ) * (p 2 * ζ) * (p 3 * ζ) = -1 := by
    linear_combination ζ ^ 4 * hp + hζ
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one]
  rw [h1, h2]
  simp [C_neg, C_1]
  ring

theorem charpoly_gQ (hζ : ζ ^ 4 = -1) : (gQ ζ).charpoly = (X ^ 4 + 1) ^ 2 :=
  charpoly_gOrb ζ hζ (phases ζ) (phases_prod ζ hζ)

/-- `X⁴ + 1` splits into the four sector eigenvalues `ζ, ζ³, ζ⁵, ζ⁷`. -/
theorem X4_add_one_eq_prod (hζ : ζ ^ 4 = -1) :
    (X ^ 4 + 1 : R[X]) =
      (X - C ζ) * (X - C (ζ ^ 3)) * (X - C (ζ ^ 5)) * (X - C (ζ ^ 7)) := by
  have hz : (C ζ : R[X]) ^ 4 = -1 := by rw [← C_pow, hζ, C_neg, C_1]
  simp only [C_pow]
  linear_combination
    (-(-X ^ 3 * C ζ ^ 3 - X ^ 3 * C ζ + X ^ 2 * C ζ ^ 6 - X * C ζ ^ 11 - X * C ζ ^ 9 +
      C ζ ^ 12 + C ζ ^ 8 * (X ^ 2 - 1) + C ζ ^ 4 * (X ^ 2 + 1) - 1)) * hz

/-- Hence `charpoly g̃ = ∏_{m=0}^{3} (X - ζ^{2m+1})²`: each sector has algebraic multiplicity 2. -/
theorem charpoly_gQ_sectors (hζ : ζ ^ 4 = -1) :
    (gQ ζ).charpoly =
      ((X - C ζ) * (X - C (ζ ^ 3)) * (X - C (ζ ^ 5)) * (X - C (ζ ^ 7))) ^ 2 := by
  rw [charpoly_gQ ζ hζ, X4_add_one_eq_prod ζ hζ]

end QuarterTurn

/-! ### Spectral sectors: each eigenspace is exactly one quarter -/

section Sectors

variable {R : Type*} [CommRing R] {A : Type*} [Ring A] [Algebra R A]

/-- `4 •` the spectral projector of `g` (with `g⁴ = -1`) onto the eigenvalue `l` (`l⁴ = -1`):
`Σ_{k<4} (l⁻¹ g)ᵏ` with `l⁻¹ = -l³`. -/
def sectorQ (g : A) (l : R) : A :=
  (1 : A) + (-l ^ 3) • g + (-l ^ 3) ^ 2 • g ^ 2 + (-l ^ 3) ^ 3 • g ^ 3

theorem mul_sectorQ (g : A) (hg : g ^ 4 = -1) (l : R) (hl : l ^ 4 = -1) :
    g * sectorQ g l = l • sectorQ g l := by
  have e2 : g * g = g ^ 2 := (sq g).symm
  have e3 : g * g ^ 2 = g ^ 3 := (pow_succ' g 2).symm
  have e4 : g * g ^ 3 = -1 := by rw [← pow_succ', hg]
  simp only [sectorQ, mul_add, mul_one, mul_smul_comm, e2, e3, e4, smul_add, smul_smul, smul_neg]
  have c1 : l * -l ^ 3 = 1 := by linear_combination -hl
  have c2 : l * (-l ^ 3) ^ 2 = -l ^ 3 := by linear_combination (l ^ 3) * hl
  have c3 : l * (-l ^ 3) ^ 3 = (-l ^ 3) ^ 2 := by linear_combination (-l ^ 6) * hl
  have c4 : (-l ^ 3) ^ 3 = -l := by linear_combination (-l ^ 5 + l) * hl
  rw [c1, c2, c3, c4, one_smul]
  module

theorem pow_mul_sectorQ (g : A) (hg : g ^ 4 = -1) (l : R) (hl : l ^ 4 = -1) (k : ℕ) :
    g ^ k * sectorQ g l = l ^ k • sectorQ g l := by
  induction k with
  | zero => simp
  | succ k ih => rw [pow_succ', mul_assoc, ih, mul_smul_comm, mul_sectorQ g hg l hl, smul_smul,
      ← pow_succ]

theorem sectorQ_mul_sectorQ (g : A) (hg : g ^ 4 = -1) (l : R) (hl : l ^ 4 = -1) :
    sectorQ g l * sectorQ g l = (4 : R) • sectorQ g l := by
  have h2 := pow_mul_sectorQ g hg l hl 2
  have h3 := pow_mul_sectorQ g hg l hl 3
  rw [show sectorQ g l * sectorQ g l = ((1 : A) + (-l ^ 3) • g + (-l ^ 3) ^ 2 • g ^ 2 +
      (-l ^ 3) ^ 3 • g ^ 3) * sectorQ g l from rfl]
  simp only [add_mul, one_mul, smul_mul_assoc, mul_sectorQ g hg l hl, h2, h3, smul_smul]
  have c1 : -l ^ 3 * l = 1 := by linear_combination -hl
  have c2 : (-l ^ 3) ^ 2 * l ^ 2 = 1 := by linear_combination (l ^ 4 - 1) * hl
  have c3 : (-l ^ 3) ^ 3 * l ^ 3 = 1 := by linear_combination (-l ^ 8 + l ^ 4 - 1) * hl
  rw [c1, c2, c3]
  module

end Sectors

section SectorDimension

/-- For an operator with `G⁴ = -1` and `tr G = tr G² = tr G³ = 0` (char 0), every eigenspace
for an eigenvalue `l` with `l⁴ = -1` has exactly a quarter of the dimension. -/
theorem four_mul_finrank_eigenspace {K V : Type*} [Field K] [CharZero K] [AddCommGroup V]
    [Module K V] [FiniteDimensional K V] (G : Module.End K V) (hG : G ^ 4 = -1)
    (h1 : LinearMap.trace K V G = 0) (h2 : LinearMap.trace K V (G ^ 2) = 0)
    (h3 : LinearMap.trace K V (G ^ 3) = 0) (l : K) (hl : l ^ 4 = -1) :
    4 * Module.finrank K (Module.End.eigenspace G l) = Module.finrank K V := by
  set P : Module.End K V := (4 : K)⁻¹ • sectorQ G l with hP
  have h4 : (4 : K) ≠ 0 := by norm_num
  have hPP : IsIdempotentElem P := by
    unfold IsIdempotentElem
    rw [hP, smul_mul_smul_comm, sectorQ_mul_sectorQ G hG l hl, smul_smul]
    congr 1
    field_simp
  have hGP : G * P = l • P := by
    rw [hP, mul_smul_comm, mul_sectorQ G hG l hl, smul_comm]
  have hrange : Module.End.eigenspace G l = LinearMap.range P := by
    ext v
    rw [Module.End.mem_eigenspace_iff]
    constructor
    · intro hv
      have hk : ∀ k : ℕ, (G ^ k) v = l ^ k • v := by
        intro k
        induction k with
        | zero => simp
        | succ k ih => rw [pow_succ', Module.End.mul_apply, ih, map_smul, hv, smul_smul, pow_succ]
      refine ⟨v, ?_⟩
      simp only [hP, sectorQ, LinearMap.smul_apply, LinearMap.add_apply, Module.End.one_apply,
        hk 2, hk 3, hv, smul_smul]
      have c1 : -l ^ 3 * l = 1 := by linear_combination -hl
      have c2 : (-l ^ 3) ^ 2 * l ^ 2 = 1 := by linear_combination (l ^ 4 - 1) * hl
      have c3 : (-l ^ 3) ^ 3 * l ^ 3 = 1 := by linear_combination (-l ^ 8 + l ^ 4 - 1) * hl
      rw [c1, c2, c3]
      rw [show v + (1 : K) • v + (1 : K) • v + (1 : K) • v = (4 : K) • v by module, smul_smul,
        inv_mul_cancel₀ h4, one_smul]
    · rintro ⟨w, rfl⟩
      have := congrArg (fun F : Module.End K V => F w) hGP
      simpa using this
  have htr : LinearMap.trace K V P = (4 : K)⁻¹ * Module.finrank K V := by
    simp [hP, sectorQ, map_add, map_smul, h1, h2, h3, LinearMap.trace_one]
  have hproj := (LinearMap.isProj_range_iff_isIdempotentElem P).mpr hPP
  have htr' := hproj.trace
  rw [hrange]
  have h' : ((4 * Module.finrank K (LinearMap.range P) : ℕ) : K) = (Module.finrank K V : K) := by
    push_cast
    rw [← htr', htr]
    field_simp
  exact_mod_cast h'

/-- With `G⁴ = -1` (char 0) there is no invariant sector: the eigenvalue `1` does not occur. -/
theorem eigenspace_one_trivial {K V : Type*} [Field K] [CharZero K] [AddCommGroup V]
    [Module K V] (G : Module.End K V) (hG : G ^ 4 = -1) :
    Module.End.eigenspace G 1 = ⊥ := by
  rw [eq_bot_iff]
  intro v hv
  rw [Module.End.mem_eigenspace_iff, one_smul] at hv
  have h4 : (G ^ 4) v = v := by
    simp only [pow_succ, Module.End.mul_apply, hv, pow_zero, Module.End.one_apply]
  rw [hG] at h4
  simp only [LinearMap.neg_apply, Module.End.one_apply] at h4
  have : (2 : K) • v = 0 := by rw [two_smul]; nth_rewrite 1 [← h4]; simp
  rw [Submodule.mem_bot]
  exact (smul_eq_zero.mp this).resolve_left two_ne_zero

theorem commute_sectorQ {K V : Type*} [Field K] [AddCommGroup V] [Module K V]
    (G D : Module.End K V) (hGD : Commute G D) (l : K) : Commute (sectorQ G l) D := by
  unfold sectorQ
  refine (((Commute.one_left D).add_left ?_).add_left ?_).add_left ?_
  · exact hGD.smul_left _
  · exact (hGD.pow_left 2).smul_left _
  · exact (hGD.pow_left 3).smul_left _

/-- Joint sectors: with `[G, D] = 0`, `G⁴ = -1`, `D² = s²`, and vanishing traces of
`Gᵏ` (`k = 1,2,3`) and `Gᵏ D` (`k = 0,…,3`), every joint eigenspace has `1/8` of the
dimension; in particular each sector carries exactly a quarter of each eigenspace of `D`. -/
theorem eight_mul_finrank_joint {K V : Type*} [Field K] [CharZero K] [AddCommGroup V]
    [Module K V] [FiniteDimensional K V] (G D : Module.End K V) (hGD : Commute G D)
    (hG : G ^ 4 = -1)
    (h1 : LinearMap.trace K V G = 0) (h2 : LinearMap.trace K V (G ^ 2) = 0)
    (h3 : LinearMap.trace K V (G ^ 3) = 0)
    (d0 : LinearMap.trace K V D = 0) (d1 : LinearMap.trace K V (G * D) = 0)
    (d2 : LinearMap.trace K V (G ^ 2 * D) = 0) (d3 : LinearMap.trace K V (G ^ 3 * D) = 0)
    (s : K) (hs : s ≠ 0) (hD : D ^ 2 = s ^ 2 • (1 : Module.End K V)) (l : K) (hl : l ^ 4 = -1) :
    8 * Module.finrank K ↥(Module.End.eigenspace G l ⊓ Module.End.eigenspace D s) =
      Module.finrank K V := by
  set P : Module.End K V := (4 : K)⁻¹ • sectorQ G l with hP
  set Pd : Module.End K V := (2 : K)⁻¹ • (1 + s⁻¹ • D) with hPdDef
  have h4 : (4 : K) ≠ 0 := by norm_num
  have h2' : (2 : K) ≠ 0 := by norm_num
  have hPP : P * P = P := by
    rw [hP, smul_mul_smul_comm, sectorQ_mul_sectorQ G hG l hl, smul_smul]
    congr 1
    field_simp
  have hGP : G * P = l • P := by
    rw [hP, mul_smul_comm, mul_sectorQ G hG l hl, smul_comm]
  have hDPd : D * Pd = s • Pd := by
    rw [hPdDef, mul_smul_comm, mul_add, mul_one, mul_smul_comm, ← sq, hD, smul_smul,
      show s⁻¹ * s ^ 2 = s by field_simp]
    match_scalars <;> field_simp
  have hPdPd : Pd * Pd = Pd := by
    have : Pd * Pd = (2 : K)⁻¹ • (Pd + s⁻¹ • (D * Pd)) := by
      rw [hPdDef]
      rw [smul_mul_assoc, add_mul, one_mul, smul_mul_assoc]
    rw [this, hDPd, smul_smul, inv_mul_cancel₀ hs, one_smul, ← two_smul K Pd, smul_smul,
      inv_mul_cancel₀ h2', one_smul]
  have hPD : Commute P D := (commute_sectorQ G D hGD l).smul_left _
  have hPPd : Commute P Pd := by
    rw [hPdDef]
    exact (((Commute.one_right P).add_right (hPD.smul_right _))).smul_right _
  set E : Module.End K V := P * Pd with hE
  have hEE : IsIdempotentElem E := by
    unfold IsIdempotentElem
    rw [hE, mul_assoc, ← mul_assoc Pd P Pd, ← hPPd.eq, mul_assoc, hPdPd, ← mul_assoc, hPP]
  have hrange : Module.End.eigenspace G l ⊓ Module.End.eigenspace D s = LinearMap.range E := by
    ext v
    rw [Submodule.mem_inf, Module.End.mem_eigenspace_iff, Module.End.mem_eigenspace_iff]
    constructor
    · rintro ⟨hv, hdv⟩
      have hk : ∀ k : ℕ, (G ^ k) v = l ^ k • v := by
        intro k
        induction k with
        | zero => simp
        | succ k ih => rw [pow_succ', Module.End.mul_apply, ih, map_smul, hv, smul_smul, pow_succ]
      have hPdv : Pd v = v := by
        rw [hPdDef]
        simp only [LinearMap.smul_apply, LinearMap.add_apply, Module.End.one_apply, hdv,
          smul_smul, inv_mul_cancel₀ hs, one_smul]
        rw [← two_smul K v, smul_smul, inv_mul_cancel₀ h2', one_smul]
      have hPv : P v = v := by
        simp only [hP, sectorQ, LinearMap.smul_apply, LinearMap.add_apply, Module.End.one_apply,
          hk 2, hk 3, hv, smul_smul]
        have c1 : -l ^ 3 * l = 1 := by linear_combination -hl
        have c2 : (-l ^ 3) ^ 2 * l ^ 2 = 1 := by linear_combination (l ^ 4 - 1) * hl
        have c3 : (-l ^ 3) ^ 3 * l ^ 3 = 1 := by linear_combination (-l ^ 8 + l ^ 4 - 1) * hl
        rw [c1, c2, c3]
        rw [show v + (1 : K) • v + (1 : K) • v + (1 : K) • v = (4 : K) • v by module, smul_smul,
          inv_mul_cancel₀ h4, one_smul]
      exact ⟨v, by rw [hE, Module.End.mul_apply, hPdv, hPv]⟩
    · rintro ⟨w, rfl⟩
      constructor
      · have := congrArg (fun F : Module.End K V => F (Pd w)) hGP
        simpa [hE, Module.End.mul_apply] using this
      · have h' : D * E = s • E := by
          rw [hE, ← mul_assoc, ← hPD.eq, mul_assoc, hDPd, mul_smul_comm]
        have := congrArg (fun F : Module.End K V => F w) h'
        simpa using this
  have htr : LinearMap.trace K V E = (8 : K)⁻¹ * Module.finrank K V := by
    rw [hE, hP, hPdDef]
    simp only [mul_add, mul_one, mul_smul_comm, sectorQ, add_mul, one_mul,
      smul_mul_assoc, map_smul, map_add, LinearMap.trace_one, h1, h2, h3, d0, d1, d2, d3]
    simp
    field_simp
    ring
  have hproj := (LinearMap.isProj_range_iff_isIdempotentElem E).mpr hEE
  have htr' := hproj.trace
  rw [hrange]
  have h' : ((8 * Module.finrank K (LinearMap.range E) : ℕ) : K) = (Module.finrank K V : K) := by
    push_cast
    rw [← htr', htr]
    field_simp
  exact_mod_cast h'

theorem sector_eigenvalue_pow_four {K : Type*} [Field K] (ζ : K) (hζ : ζ ^ 4 = -1) (m : Fin 4) :
    (ζ ^ (2 * m.val + 1)) ^ 4 = -1 := by
  rw [← pow_mul, show (2 * m.val + 1) * 4 = 4 * (2 * m.val + 1) by ring, pow_mul, hζ,
    Odd.neg_one_pow ⟨m.val, rfl⟩]

/-- **Each of the four spinor sectors is exactly one quarter** on every orbit (phases with
product one): for every `m`, the eigenspace of `g̃` for `ζ^{2m+1}` has dimension `2` (of `8`). -/
theorem finrank_sector_gOrb {K : Type*} [Field K] [CharZero K] (ζ : K) (hζ : ζ ^ 4 = -1)
    (p : Fin 4 → K) (hp : p 0 * p 1 * p 2 * p 3 = 1) (m : Fin 4) :
    Module.finrank K
      (Module.End.eigenspace (Matrix.toLin' (gOrb ζ p)) (ζ ^ (2 * m.val + 1))) = 2 := by
  have hG : (Matrix.toLin' (gOrb ζ p)) ^ 4 = -1 := by
    rw [← Matrix.toLin'_pow, gOrb_pow_four ζ hζ p hp, map_neg, Matrix.toLin'_one]
    rfl
  have ht : ∀ k : ℕ, 1 ≤ k → k ≤ 3 →
      LinearMap.trace K _ ((Matrix.toLin' (gOrb ζ p)) ^ k) = 0 := by
    intro k hk1 hk3
    rw [← Matrix.toLin'_pow, Matrix.trace_toLin'_eq, trace_gOrb_pow ζ p k hk1 hk3]
  have h := four_mul_finrank_eigenspace (Matrix.toLin' (gOrb ζ p)) hG
    (by simpa using ht 1 le_rfl (by norm_num)) (ht 2 (by norm_num) (by norm_num))
    (ht 3 (by norm_num) le_rfl) _ (sector_eigenvalue_pow_four ζ hζ m)
  rw [Module.finrank_fintype_fun_eq_card] at h
  simp at h
  omega

theorem finrank_sector_gQ {K : Type*} [Field K] [CharZero K] (ζ : K) (hζ : ζ ^ 4 = -1)
    (m : Fin 4) :
    Module.finrank K (Module.End.eigenspace (Matrix.toLin' (gQ ζ)) (ζ ^ (2 * m.val + 1))) = 2 :=
  finrank_sector_gOrb ζ hζ (phases ζ) (phases_prod ζ hζ) m

/-- **Each sector carries exactly a quarter of each eigenspace of `D`** (z1, point 2): on the
orbit of any `(u, v)` with `s² = u² + v² ≠ 0`, the joint eigenspace of `g̃` (eigenvalue
`ζ^{2m+1}`) and `D` (eigenvalue `s`, either sign) has dimension `1`; the `D`-eigenspace for `s`
has dimension `4`. -/
theorem finrank_joint_gOrb {K : Type*} [Field K] [CharZero K] (ζ : K) (hζ : ζ ^ 4 = -1)
    (p : Fin 4 → K) (hp : p 0 * p 1 * p 2 * p 3 = 1) (u v s : K) (hs : s ≠ 0)
    (hsuv : s ^ 2 = u ^ 2 + v ^ 2) (m : Fin 4) :
    Module.finrank K ↥(Module.End.eigenspace (Matrix.toLin' (gOrb ζ p)) (ζ ^ (2 * m.val + 1)) ⊓
      Module.End.eigenspace (Matrix.toLin' (diracOrb ζ u v)) s) = 1 := by
  have hG : (Matrix.toLin' (gOrb ζ p)) ^ 4 = -1 := by
    rw [← Matrix.toLin'_pow, gOrb_pow_four ζ hζ p hp, map_neg, Matrix.toLin'_one]
    rfl
  have hGD : Commute (Matrix.toLin' (gOrb ζ p)) (Matrix.toLin' (diracOrb ζ u v)) := by
    show _ * _ = _ * _
    rw [Module.End.mul_eq_comp, Module.End.mul_eq_comp, ← Matrix.toLin'_mul,
      ← Matrix.toLin'_mul, gOrb_comm_diracOrb ζ hζ p u v]
  have hD : (Matrix.toLin' (diracOrb ζ u v)) ^ 2 = s ^ 2 • (1 : Module.End K _) := by
    rw [← Matrix.toLin'_pow, diracOrb_sq ζ hζ u v, ← hsuv, map_smul, Matrix.toLin'_one]
    rfl
  have ht : ∀ k : ℕ, 1 ≤ k → k ≤ 3 →
      LinearMap.trace K _ ((Matrix.toLin' (gOrb ζ p)) ^ k) = 0 := by
    intro k hk1 hk3
    rw [← Matrix.toLin'_pow, Matrix.trace_toLin'_eq, trace_gOrb_pow ζ p k hk1 hk3]
  have hd : ∀ k : ℕ, k ≤ 3 → LinearMap.trace K _
      ((Matrix.toLin' (gOrb ζ p)) ^ k * Matrix.toLin' (diracOrb ζ u v)) = 0 := by
    intro k hk
    rw [← Matrix.toLin'_pow, Module.End.mul_eq_comp, ← Matrix.toLin'_mul,
      Matrix.trace_toLin'_eq, trace_gOrb_pow_mul_diracOrb ζ p u v k hk]
  have h := eight_mul_finrank_joint (Matrix.toLin' (gOrb ζ p)) (Matrix.toLin' (diracOrb ζ u v))
    hGD hG (by simpa using ht 1 le_rfl (by norm_num)) (ht 2 (by norm_num) (by norm_num))
    (ht 3 (by norm_num) le_rfl) (by simpa using hd 0 (by norm_num))
    (by simpa using hd 1 (by norm_num)) (hd 2 (by norm_num)) (hd 3 le_rfl) s hs hD _
    (sector_eigenvalue_pow_four ζ hζ m)
  rw [Module.finrank_fintype_fun_eq_card] at h
  simp at h
  omega

end SectorDimension

/-! ### Identification in `ℂ`: `ζ = e^{iπ/4}` -/

section Complex

open Complex Real

/-- `ζ₈ = e^{iπ/4}`. -/
noncomputable def zeta8 : ℂ := Complex.exp (π / 4 * I)

theorem zeta8_pow_four : zeta8 ^ 4 = -1 := by
  rw [zeta8, ← Complex.exp_nat_mul]
  rw [show ((4 : ℕ) : ℂ) * (π / 4 * I) = π * I by push_cast; ring]
  exact Complex.exp_pi_mul_I

theorem zeta8_sq : zeta8 ^ 2 = I := by
  rw [zeta8, ← Complex.exp_nat_mul]
  rw [show ((2 : ℕ) : ℂ) * (π / 4 * I) = π / 2 * I by push_cast; ring]
  exact Complex.exp_pi_div_two_mul_I

/-- The spinor rotation is `diag(e^{-iπ/4}, e^{iπ/4})`, as in `z1`. -/
theorem spinRot_zeta8 :
    spinRot zeta8 = diagonal ![Complex.exp (-(π / 4 * I)), Complex.exp (π / 4 * I)] := by
  have h : -zeta8 ^ 3 = Complex.exp (-(π / 4 * I)) := by
    rw [Complex.exp_neg, ← zeta8]
    have hz : zeta8 ≠ 0 := Complex.exp_ne_zero _
    field_simp
    linear_combination (-1 : ℂ) * zeta8_pow_four
  rw [spinRot, h, zeta8]

/-- The phases are `(-i, -i, i, i) = exp(-2πi k_j·a)` for `a = (1/2, 0)` and the orbit
`kx = (1/2, -3/2, -1/2, 3/2)` (that is `kx2 / 2`), as in `z1`. -/
theorem phases_zeta8 :
    phases zeta8 = ![Complex.exp (-(2 * π * I) * ((1 / 2 : ℂ) * (1 / 2))),
                     Complex.exp (-(2 * π * I) * ((-3 / 2 : ℂ) * (1 / 2))),
                     Complex.exp (-(2 * π * I) * ((-1 / 2 : ℂ) * (1 / 2))),
                     Complex.exp (-(2 * π * I) * ((3 / 2 : ℂ) * (1 / 2)))] := by
  have hI : Complex.exp (π / 2 * I) = I := Complex.exp_pi_div_two_mul_I
  have hpi : Complex.exp (π * I) = -1 := Complex.exp_pi_mul_I
  have e0 : Complex.exp (-(2 * π * I) * ((1 / 2 : ℂ) * (1 / 2))) = -I := by
    rw [show (-(2 * π * I) * ((1 / 2 : ℂ) * (1 / 2)) : ℂ) = -(π / 2 * I) by ring,
      Complex.exp_neg, hI, Complex.inv_I]
  have e1 : Complex.exp (-(2 * π * I) * ((-3 / 2 : ℂ) * (1 / 2))) = -I := by
    rw [show (-(2 * π * I) * ((-3 / 2 : ℂ) * (1 / 2)) : ℂ) = π * I + π / 2 * I by ring,
      Complex.exp_add, hI, hpi]
    ring
  have e2 : Complex.exp (-(2 * π * I) * ((-1 / 2 : ℂ) * (1 / 2))) = I := by
    rw [show (-(2 * π * I) * ((-1 / 2 : ℂ) * (1 / 2)) : ℂ) = π / 2 * I by ring, hI]
  have e3 : Complex.exp (-(2 * π * I) * ((3 / 2 : ℂ) * (1 / 2))) = I := by
    rw [show (-(2 * π * I) * ((3 / 2 : ℂ) * (1 / 2)) : ℂ) = -(π * I + π / 2 * I) by ring,
      Complex.exp_neg, Complex.exp_add, hI, hpi]
    rw [show (-1 * I : ℂ) = -I by ring, inv_neg, Complex.inv_I, neg_neg]
  rw [e0, e1, e2, e3, phases, zeta8_sq]

theorem kx2_half : ((kx2 0 : ℂ) / 2, (kx2 1 : ℂ) / 2, (kx2 2 : ℂ) / 2, (kx2 3 : ℂ) / 2) =
    (1 / 2, -3 / 2, -1 / 2, 3 / 2) := by
  simp [kx2, orbX]

/-- **z1, point 2, in `ℂ`:** each of the four sectors `e^{iπ(2m+1)/4}` of the lifted quarter
turn has dimension exactly `2` of `8`. -/
theorem complex_sector_dim (m : Fin 4) :
    Module.finrank ℂ
      (Module.End.eigenspace (Matrix.toLin' (gQ zeta8)) (zeta8 ^ (2 * m.val + 1))) = 2 :=
  finrank_sector_gQ zeta8 zeta8_pow_four m

theorem complex_no_invariant_sector :
    Module.End.eigenspace (Matrix.toLin' (gQ zeta8)) 1 = ⊥ := by
  apply eigenspace_one_trivial
  rw [← Matrix.toLin'_pow, gQ_pow_four zeta8 zeta8_pow_four, map_neg, Matrix.toLin'_one]
  rfl

/-- The translation phases of any orbit multiply to one (`Σ_j k_j = 0`). -/
theorem orbit_phase_product (u v ax ay : ℂ) :
    (∏ j : Fin 4, Complex.exp (-(2 * π * I) * (orbX u v j * ax + orbY u v j * ay))) = 1 := by
  rw [← Complex.exp_sum]
  convert Complex.exp_zero using 2
  simp [Fin.sum_univ_four, orbX, orbY]
  ring

/-- **z1, point 2, for every orbit in `ℂ`:** on the orbit of any `(u, v)` with translation
`a = (ax, ay)` and `s² = u² + v² ≠ 0`, each spinor sector `e^{iπ(2m+1)/4}` meets the
eigenspace of `D` for `s` in exactly one dimension.  This is the exact finite content of
"jeder Sektor trägt exakt ein Viertel jedes Eigenraums". -/
theorem complex_every_orbit_quarter (u v ax ay s : ℂ) (hs : s ≠ 0) (hsuv : s ^ 2 = u ^ 2 + v ^ 2)
    (m : Fin 4) :
    Module.finrank ℂ ↥(Module.End.eigenspace (Matrix.toLin' (gOrb zeta8
        (fun j => Complex.exp (-(2 * π * I) * (orbX u v j * ax + orbY u v j * ay)))))
        (zeta8 ^ (2 * m.val + 1)) ⊓
      Module.End.eigenspace (Matrix.toLin' (diracOrb zeta8 u v)) s) = 1 := by
  apply finrank_joint_gOrb zeta8 zeta8_pow_four _ _ u v s hs hsuv m
  have h := orbit_phase_product u v ax ay
  rw [Fin.prod_univ_four] at h
  exact h

/-- The orbit of `z1` (`k₀ = (1/2, 3/2)`, doubled `(1, 3)`, `|2k|² = 10`): each sector meets
each `D`-eigenspace `±√10` in exactly one dimension. -/
theorem complex_z1_orbit_quarter (m : Fin 4) (sgn : ℂ) (hsgn : sgn = 1 ∨ sgn = -1) :
    Module.finrank ℂ ↥(Module.End.eigenspace (Matrix.toLin' (gQ zeta8))
        (zeta8 ^ (2 * m.val + 1)) ⊓
      Module.End.eigenspace (Matrix.toLin' (dirac2 zeta8)) (sgn * (Real.sqrt 10 : ℂ))) = 1 := by
  have h10 : ((Real.sqrt 10 : ℝ) : ℂ) ^ 2 = 10 := by
    rw [← Complex.ofReal_pow, Real.sq_sqrt (by norm_num)]
    norm_num
  have hs : sgn * (Real.sqrt 10 : ℂ) ≠ 0 := by
    apply mul_ne_zero
    · rcases hsgn with rfl | rfl <;> norm_num
    · exact_mod_cast (Real.sqrt_pos.mpr (by norm_num : (0 : ℝ) < 10)).ne'
  have hsgn2 : sgn ^ 2 = 1 := by rcases hsgn with rfl | rfl <;> norm_num
  apply finrank_joint_gOrb zeta8 zeta8_pow_four (phases zeta8) (phases_prod zeta8 zeta8_pow_four)
    1 3 _ hs _ m
  rw [mul_pow, hsgn2, h10]
  norm_num

/-- The four sector eigenvalues are pairwise different. -/
theorem zeta8_sectors_injective :
    Function.Injective fun m : Fin 4 => zeta8 ^ (2 * m.val + 1) := by
  have h4 := zeta8_pow_four
  have h2 := zeta8_sq
  have hz : zeta8 ≠ 0 := Complex.exp_ne_zero _
  intro a b hab
  simp only at hab
  have key : ∀ x y : ℕ, zeta8 ^ x = zeta8 ^ (x + y) → zeta8 ^ y = 1 := by
    intro x y hxy
    rw [pow_add] at hxy
    exact (mul_eq_left₀ (pow_ne_zero x hz)).mp hxy.symm
  have n2 : zeta8 ^ 2 ≠ 1 := by
    rw [h2]
    intro h
    exact absurd (congrArg Complex.im h) (by simp)
  have n4 : zeta8 ^ 4 ≠ 1 := by rw [h4]; norm_num
  have n6 : zeta8 ^ 6 ≠ 1 := by
    rw [show (6 : ℕ) = 4 + 2 from rfl, pow_add, h4, h2]
    intro h
    have : I = -1 := by linear_combination -h
    exact absurd (congrArg Complex.im this) (by simp)
  fin_cases a <;> fin_cases b <;> simp at hab ⊢
  · exact n2 (key 1 2 (by simpa using hab))
  · exact n4 (key 1 4 (by simpa using hab))
  · exact n6 (key 1 6 (by simpa using hab))
  · exact n2 (key 1 2 (by simpa using hab.symm))
  · exact n2 (key 3 2 (by simpa using hab))
  · exact n4 (key 3 4 (by simpa using hab))
  · exact n4 (key 1 4 (by simpa using hab.symm))
  · exact n2 (key 3 2 (by simpa using hab.symm))
  · exact n2 (key 5 2 (by simpa using hab))
  · exact n6 (key 1 6 (by simpa using hab.symm))
  · exact n4 (key 3 4 (by simpa using hab.symm))
  · exact n2 (key 5 2 (by simpa using hab.symm))

end Complex

/-! ## B. Family characters from the mark index (z10, point 3) -/

section Families

open Complex

/-- Pullback of `ω = P(z) dz/(z⁴-1)` under `z ↦ i z`: the coefficients become `i^{k+1} P_k`. -/
theorem quarter_on_differential_coeff (P : ℂ[X]) (k : ℕ) :
    (Polynomial.C I * P.comp (Polynomial.C I * X)).coeff k = I ^ (k + 1) * P.coeff k := by
  rw [coeff_C_mul, comp_C_mul_X_coeff]; ring

/-- Pullback of `f = Q(z)/(z⁴-1)`: the coefficients become `i^k Q_k`. -/
theorem quarter_on_function_coeff (Q : ℂ[X]) (k : ℕ) :
    (Q.comp (Polynomial.C I * X)).coeff k = I ^ k * Q.coeff k := by
  rw [comp_C_mul_X_coeff]; ring

/-- The mark polynomial `z⁴ - 1` (marks `D = μ₄`) is invariant under `z ↦ i z`. -/
theorem marks_polynomial_invariant :
    (X ^ 4 - 1 : ℂ[X]).comp (Polynomial.C I * X) = X ^ 4 - 1 := by
  simp [mul_pow, ← Polynomial.C_pow, Complex.I_pow_four]

/-- The space `deg P ≤ 2` (model of `H⁰(K(D))`) is invariant. -/
theorem quarter_on_differential_degree (P : ℂ[X]) (hP : P.natDegree ≤ 2) :
    (Polynomial.C I * P.comp (Polynomial.C I * X)).natDegree ≤ 2 := by
  rw [natDegree_le_iff_coeff_eq_zero]
  intro k hk
  rw [quarter_on_differential_coeff, coeff_eq_zero_of_natDegree_lt (lt_of_le_of_lt hP hk),
    mul_zero]

/-- Matrix of `P(z) ↦ i P(i z)` on the basis `{1, z, z²}`. -/
noncomputable def famMatrix : Matrix (Fin 3) (Fin 3) ℂ := diagonal fun k => I ^ (k.val + 1)

theorem famMatrix_eq : famMatrix = diagonal ![I, -1, -I] := by
  unfold famMatrix
  congr 1
  funext k
  fin_cases k <;> simp [pow_succ]

theorem famMatrix_represents (P : ℂ[X]) :
    (fun k : Fin 3 => (Polynomial.C I * P.comp (Polynomial.C I * X)).coeff k) =
      famMatrix *ᵥ (fun k : Fin 3 => P.coeff k) := by
  funext k
  rw [famMatrix, mulVec_diagonal, quarter_on_differential_coeff]

/-- **Family characters `{i, -1, -i}`.** -/
theorem famMatrix_charpoly :
    famMatrix.charpoly = (X - Polynomial.C I) * (X + 1) * (X + Polynomial.C I) := by
  rw [famMatrix_eq, charpoly_diagonal, Fin.prod_univ_three]
  simp [sub_eq_add_neg]

/-- Matrix of `Q(z) ↦ Q(i z)` on the basis `{1, z, …, z⁴}` (model of `H⁰(O(D))`). -/
noncomputable def carrierMatrix : Matrix (Fin 5) (Fin 5) ℂ := diagonal fun k => I ^ k.val

theorem carrierMatrix_represents (Q : ℂ[X]) :
    (fun k : Fin 5 => (Q.comp (Polynomial.C I * X)).coeff k) =
      carrierMatrix *ᵥ (fun k : Fin 5 => Q.coeff k) := by
  funext k
  rw [carrierMatrix, mulVec_diagonal, quarter_on_function_coeff]

/-- Half turn on the carrier: `diag(1, -1, 1, -1, 1)`. -/
theorem carrierMatrix_sq : carrierMatrix ^ 2 = diagonal ![1, -1, 1, -1, 1] := by
  rw [carrierMatrix, diagonal_pow]
  congr 1
  funext k
  fin_cases k <;> simp [pow_succ]

theorem carrierMatrix_pow_four : carrierMatrix ^ 4 = 1 := by
  rw [carrierMatrix, diagonal_pow, ← diagonal_one]
  congr 1
  funext k
  fin_cases k <;> simp [pow_succ]

/-- **Split `3 + 2`** of the carrier under the half turn. -/
theorem carrier_half_turn_charpoly :
    (carrierMatrix ^ 2).charpoly = (X - 1) ^ 3 * (X + 1) ^ 2 := by
  rw [carrierMatrix_sq, charpoly_diagonal, Fin.prod_univ_five]
  simp [sub_eq_add_neg]
  ring

theorem carrier_split_three_two :
    (Finset.univ.filter fun k : Fin 5 => Even k.val).card = 3 ∧
    (Finset.univ.filter fun k : Fin 5 => Odd k.val).card = 2 := by
  decide

/-- `Λ^even ℂ⁵ = 1 + 10 + 5 = 16`. -/
theorem lambda_even_C5 : Nat.choose 5 0 + Nat.choose 5 2 + Nat.choose 5 4 = 16 := by decide

theorem lambda_even_C5_sum : ∑ k ∈ (Finset.range 6).filter Even, Nat.choose 5 k = 16 := by
  decide

/-- `3 · 16 = 48` components of `H⁰(K(D)) ⊗ Λ^even H⁰(O(D))`. -/
theorem matter_components :
    3 * (Nat.choose 5 0 + Nat.choose 5 2 + Nat.choose 5 4) = 48 := by decide

end Families

/-! ## C. Arithmetic of central charge, twist weight and GSO phases -/

section Arithmetic

open Complex Real

/-- `c = ν/2 ∈ 8ℤ ↔ 16 ∣ ν`. -/
theorem central_charge_in_8Z_iff (ν : ℤ) : (∃ k : ℤ, (ν : ℚ) / 2 = 8 * k) ↔ 16 ∣ ν := by
  constructor
  · rintro ⟨k, hk⟩
    refine ⟨k, ?_⟩
    have h : (ν : ℚ) = 16 * k := by linarith
    exact_mod_cast h
  · rintro ⟨k, rfl⟩
    exact ⟨k, by push_cast; ring⟩

theorem exp_two_pi_I_mul_eq_one_iff (x : ℂ) :
    Complex.exp (2 * π * I * x) = 1 ↔ ∃ n : ℤ, x = n := by
  have h2 : (2 * π * I : ℂ) ≠ 0 := by simp [Real.pi_ne_zero, I_ne_zero]
  rw [Complex.exp_eq_one_iff]
  constructor
  · rintro ⟨n, hn⟩
    exact ⟨n, mul_left_cancel₀ h2 (by rw [hn]; ring)⟩
  · rintro ⟨n, rfl⟩
    exact ⟨n, by ring⟩

/-- The wall phase `exp(2πi c σ/8)` with `c = ν/2` is trivial for every signature `σ`
exactly when `16 ∣ ν` (z5, z8). -/
theorem cut_phase_trivial_iff (ν : ℤ) :
    (∀ σ : ℤ, Complex.exp (2 * π * I * (((ν : ℂ) / 2) * σ / 8)) = 1) ↔ 16 ∣ ν := by
  constructor
  · intro h
    obtain ⟨n, hn⟩ := (exp_two_pi_I_mul_eq_one_iff _).mp (h 1)
    refine ⟨n, ?_⟩
    have h' : (ν : ℂ) = 16 * n := by
      push_cast at hn
      linear_combination 16 * hn
    exact_mod_cast h'
  · rintro ⟨k, rfl⟩ σ
    rw [exp_two_pi_I_mul_eq_one_iff]
    exact ⟨k * σ, by push_cast; ring⟩

/-- Ramond twist weight `h_R = N/16` equals `1` only for `N = 16` (z10, point 4). -/
theorem twist_weight_one_iff (N : ℕ) : (N : ℚ) / 16 = 1 ↔ N = 16 := by
  constructor
  · intro h
    have h' : (N : ℚ) = 16 := by linarith [h]
    exact_mod_cast h'
  · rintro rfl; norm_num

/-- T-covariance of the GSO sum: the phase `e^{iπN/8}` of `θ₂^{N/2}` is `1` iff `16 ∣ N`. -/
theorem T_covariance_iff (N : ℤ) : Complex.exp (π * I * N / 8) = 1 ↔ 16 ∣ N := by
  rw [show (π * I * N / 8 : ℂ) = 2 * π * I * ((N : ℂ) / 16) by ring,
    exp_two_pi_I_mul_eq_one_iff]
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨n, ?_⟩
    have h' : (N : ℂ) = 16 * n := by linear_combination 16 * hn
    exact_mod_cast h'
  · rintro ⟨k, rfl⟩
    exact ⟨k, by push_cast; ring⟩

/-- Strict T-invariance: `e^{-2πi c/24} = 1` with `c = N/2` iff `48 ∣ N`. -/
theorem T_invariance_iff (N : ℤ) :
    Complex.exp (-(2 * π * I) * (((N : ℂ) / 2) / 24)) = 1 ↔ 48 ∣ N := by
  rw [show (-(2 * π * I) * (((N : ℂ) / 2) / 24) : ℂ) = 2 * π * I * (-(N : ℂ) / 48) by ring,
    exp_two_pi_I_mul_eq_one_iff]
  constructor
  · rintro ⟨n, hn⟩
    refine ⟨-n, ?_⟩
    have h' : (N : ℂ) = 48 * (-n : ℤ) := by push_cast; linear_combination -48 * hn
    exact_mod_cast h'
  · rintro ⟨k, rfl⟩
    exact ⟨-k, by push_cast; ring⟩

/-- `c = N/2 ∈ 24ℤ ↔ 48 ∣ N`. -/
theorem central_charge_in_24Z_iff (N : ℤ) : (∃ k : ℤ, (N : ℚ) / 2 = 24 * k) ↔ 48 ∣ N := by
  constructor
  · rintro ⟨k, hk⟩
    refine ⟨k, ?_⟩
    have h : (N : ℚ) = 48 * k := by linarith
    exact_mod_cast h
  · rintro ⟨k, rfl⟩
    exact ⟨k, by push_cast; ring⟩

/-- For `E₈` (`N = 16`) the T phase `e^{-2πi/3}` is not `1` … -/
theorem E8_T_phase_nontrivial :
    Complex.exp (-(2 * π * I) * (((16 : ℤ) : ℂ) / 2 / 24)) ≠ 1 := by
  rw [Ne, T_invariance_iff]
  decide

/-- … and equals `e^{4πi/3}` (argument 240°). -/
theorem E8_T_phase_eq : Complex.exp (-(2 * π * I) * (((16 : ℤ) : ℂ) / 2 / 24)) =
    Complex.exp (4 * π * I / 3) := by
  rw [show (4 * π * I / 3 : ℂ) = -(2 * π * I) * (((16 : ℤ) : ℂ) / 2 / 24) + 2 * π * I by
      push_cast; ring, Complex.exp_add, Complex.exp_two_pi_mul_I, mul_one]

/-- `n` sets of 16 (`c = 8n`): strictly T-invariant alone iff `3 ∣ n`. -/
theorem sets_strict_T_iff (n : ℕ) : 24 ∣ 8 * n ↔ 3 ∣ n := by omega

theorem sets_minimal_strict : ∀ n : ℕ, 1 ≤ n → n < 3 → ¬ 24 ∣ 8 * n := by decide

/-- Counts: `120 + 128 = 248`, `120 = 45 + 15 + 60`, `128 = 2·(16·4)`,
`248 = 78 + 8 + 2·27·3`, `16 = 10 + 6`. -/
theorem E8_counts :
    Nat.choose 16 2 + 2 ^ 7 = 248 ∧
    Nat.choose 10 2 + Nat.choose 6 2 + 10 * 6 = Nat.choose 16 2 ∧
    2 * (16 * 4) = 2 ^ 7 ∧
    78 + 8 + 2 * 27 * 3 = 248 ∧
    10 + 6 = 16 := by decide

/-- `⟨L₀⟩ = c/24 = 1/3` for `c = 8`. -/
theorem kms_energy_c8 : ((8 : ℚ) / 24) = 1 / 3 := by norm_num

end Arithmetic

/-! ## D. Charges of one Standard Model generation (z8) -/

namespace SM

/-- Left-handed Weyl multiplets of one generation (conjugates for right-handed fields). -/
inductive SMField
  | Q | uc | dc | L | ec | nuc
  deriving DecidableEq, Fintype, Repr

open SMField

/-- Number of Weyl components (colour × isospin). -/
def mult : SMField → ℕ
  | Q => 6 | uc => 3 | dc => 3 | L => 2 | ec => 1 | nuc => 1

/-- `B - L`. -/
def BL : SMField → ℚ
  | Q => 1/3 | uc => -1/3 | dc => -1/3 | L => -1 | ec => 1 | nuc => 1

/-- Hypercharge `Y`. -/
def Y : SMField → ℚ
  | Q => 1/6 | uc => -2/3 | dc => 1/3 | L => -1/2 | ec => 1 | nuc => 0

/-- `X = 5(B - L) - 4Y`. -/
def chargeX (f : SMField) : ℚ := 5 * BL f - 4 * Y f

/-- Integer values of `X`. -/
def Xint : SMField → ℤ
  | Q => 1 | uc => 1 | dc => -3 | L => -3 | ec => 1 | nuc => 5

theorem chargeX_eq_Xint (f : SMField) : chargeX f = Xint f := by
  cases f <;> norm_num [chargeX, BL, Y, Xint]

theorem Xint_mod_four (f : SMField) : Xint f % 4 = 1 := by
  cases f <;> decide

/-- **`X ≡ 1 mod 4` for every field of the generation.** -/
theorem chargeX_one_mod_four (f : SMField) : ∃ n : ℤ, chargeX f = 4 * n + 1 := by
  refine ⟨Xint f / 4, ?_⟩
  rw [chargeX_eq_Xint]
  have h4 := Xint_mod_four f
  have h : Xint f = 4 * (Xint f / 4) + 1 := by omega
  exact_mod_cast h

/-- `6 + 3 + 3 + 2 + 1 + 1 = 16`. -/
theorem generation_count : ∑ f, mult f = 16 := by decide

theorem generation_count_without_nuc :
    ∑ f ∈ Finset.univ.filter (· ≠ nuc), mult f = 15 := by decide

/-- The ℤ₁₆ count holds with `ν^c` and fails without it. -/
theorem z16_condition :
    16 ∣ ∑ f, mult f ∧ ¬ 16 ∣ ∑ f ∈ Finset.univ.filter (· ≠ nuc), mult f := by
  rw [generation_count, generation_count_without_nuc]; decide

/-- Consistency of the charge table: `Σ Y = Σ Y³ = Σ (B-L) = Σ (B-L)³ = 0`. -/
theorem anomaly_sums :
    ∑ f, (mult f : ℚ) * Y f = 0 ∧
    ∑ f, (mult f : ℚ) * Y f ^ 3 = 0 ∧
    ∑ f, (mult f : ℚ) * BL f = 0 ∧
    ∑ f, (mult f : ℚ) * BL f ^ 3 = 0 := by
  have huniv : (Finset.univ : Finset SMField) = {Q, uc, dc, L, ec, nuc} := by decide
  rw [huniv]
  simp [mult, Y, BL]
  norm_num

/-- Without `ν^c`, `Σ (B - L) = -1`. -/
theorem BL_sum_without_nuc :
    ∑ f ∈ Finset.univ.filter (· ≠ nuc), (mult f : ℚ) * BL f = -1 := by
  have huniv : Finset.univ.filter (· ≠ nuc) = ({Q, uc, dc, L, ec} : Finset SMField) := by
    decide
  rw [huniv]
  simp [mult, BL]
  norm_num

end SM

/-! ## E. `E₆` at `τ = i` (z3) -/

section Eisenstein

theorem S_smul_I : ModularGroup.S • UpperHalfPlane.I = UpperHalfPlane.I := by
  rw [UpperHalfPlane.modular_S_smul]
  ext
  simp

/-- Every level-one modular form of weight 6 vanishes at `τ = i`
(`f(-1/τ) = τ⁶ f(τ)` and `i⁶ = -1`). -/
theorem weight6_vanishes_at_I (f : ModularForm (CongruenceSubgroup.Gamma 1) 6) :
    f UpperHalfPlane.I = 0 := by
  have h := SlashInvariantForm.slash_action_eqn_SL'' f
    (CongruenceSubgroup.mem_Gamma_one ModularGroup.S) UpperHalfPlane.I
  rw [S_smul_I] at h
  simp only [UpperHalfPlane.denom, ModularGroup.S] at h
  have hI6 : Complex.I ^ (6 : ℤ) = -1 := by
    rw [show (6 : ℤ) = ((4 + 2 : ℕ) : ℤ) from rfl, zpow_natCast, pow_add, Complex.I_pow_four,
      one_mul, Complex.I_sq]
  simp at h
  rw [hI6] at h
  linear_combination h / 2

/-- In particular Mathlib's level-one Eisenstein series `E₆` vanishes at `i`. -/
theorem E6_vanishes_at_I : ModularForm.E (k := 6) (by norm_num) UpperHalfPlane.I = 0 :=
  weight6_vanishes_at_I _

/-- Arithmetic of the valence formula `k/12 = v_∞ + v_i/2 + v_ρ/3 + Σ v_p` for `k = 6`
(multiplied by 12).  The formula itself is a hypothesis here, not a Mathlib theorem:
if it holds, the only zero of a nonzero weight-6 form in the fundamental domain is a
simple zero at `i`. -/
theorem valence_weight6 (vInf vI vRho vRest : ℕ)
    (hval : 12 * vInf + 6 * vI + 4 * vRho + 12 * vRest = 6) :
    vInf = 0 ∧ vI = 1 ∧ vRho = 0 ∧ vRest = 0 := by
  omega

/-- Same for `k = 4`: the only zero of a nonzero weight-4 form (`E₄`, hence `χ_{E₈}`) is at `ρ`. -/
theorem valence_weight4 (vInf vI vRho vRest : ℕ)
    (hval : 12 * vInf + 6 * vI + 4 * vRho + 12 * vRest = 4) :
    vInf = 0 ∧ vI = 0 ∧ vRho = 1 ∧ vRest = 0 := by
  omega

end Eisenstein

end TFPT.H6Exact

/-! ## Axiom audit -/

#print axioms TFPT.H6Exact.gOrb_pow_four
#print axioms TFPT.H6Exact.gQ_pow_four
#print axioms TFPT.H6Exact.trace_gQ_pow
#print axioms TFPT.H6Exact.gOrb_comm_diracOrb
#print axioms TFPT.H6Exact.gQ_comm_dirac2
#print axioms TFPT.H6Exact.diracOrb_sq
#print axioms TFPT.H6Exact.trace_gOrb_pow_mul_diracOrb
#print axioms TFPT.H6Exact.rot2_orbit_length_four
#print axioms TFPT.H6Exact.charpoly_gOrb
#print axioms TFPT.H6Exact.charpoly_gQ
#print axioms TFPT.H6Exact.charpoly_gQ_sectors
#print axioms TFPT.H6Exact.four_mul_finrank_eigenspace
#print axioms TFPT.H6Exact.finrank_sector_gOrb
#print axioms TFPT.H6Exact.finrank_sector_gQ
#print axioms TFPT.H6Exact.spinRot_zeta8
#print axioms TFPT.H6Exact.phases_zeta8
#print axioms TFPT.H6Exact.complex_sector_dim
#print axioms TFPT.H6Exact.complex_no_invariant_sector
#print axioms TFPT.H6Exact.eight_mul_finrank_joint
#print axioms TFPT.H6Exact.finrank_joint_gOrb
#print axioms TFPT.H6Exact.complex_every_orbit_quarter
#print axioms TFPT.H6Exact.complex_z1_orbit_quarter
#print axioms TFPT.H6Exact.orbit_phase_product
#print axioms TFPT.H6Exact.zeta8_sectors_injective
#print axioms TFPT.H6Exact.famMatrix_represents
#print axioms TFPT.H6Exact.famMatrix_charpoly
#print axioms TFPT.H6Exact.marks_polynomial_invariant
#print axioms TFPT.H6Exact.carrierMatrix_represents
#print axioms TFPT.H6Exact.carrier_half_turn_charpoly
#print axioms TFPT.H6Exact.carrier_split_three_two
#print axioms TFPT.H6Exact.lambda_even_C5
#print axioms TFPT.H6Exact.matter_components
#print axioms TFPT.H6Exact.central_charge_in_8Z_iff
#print axioms TFPT.H6Exact.cut_phase_trivial_iff
#print axioms TFPT.H6Exact.twist_weight_one_iff
#print axioms TFPT.H6Exact.T_covariance_iff
#print axioms TFPT.H6Exact.T_invariance_iff
#print axioms TFPT.H6Exact.E8_T_phase_eq
#print axioms TFPT.H6Exact.sets_strict_T_iff
#print axioms TFPT.H6Exact.E8_counts
#print axioms TFPT.H6Exact.SM.chargeX_one_mod_four
#print axioms TFPT.H6Exact.SM.z16_condition
#print axioms TFPT.H6Exact.SM.anomaly_sums
#print axioms TFPT.H6Exact.weight6_vanishes_at_I
#print axioms TFPT.H6Exact.E6_vanishes_at_I
#print axioms TFPT.H6Exact.valence_weight6
#print axioms TFPT.H6Exact.valence_weight4
