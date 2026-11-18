/-
  Algebraic C11 continuation for the crossed-word construction.

  This file is intentionally independent of the production certificate.  It
  proves the degree-zero conditional-expectation identity for the concrete
  normal form a0 + a1 X.  It does not choose a physical state or construct a
  C*-completion.
-/

import Mathlib

namespace SeamStateExpectation

/-! The actual ordered eight-channel cocycle at the monomial level. -/

abbrev Phase := ZMod 4
abbrev LatticeCharge := Fin 8 → ℤ

def epsExp (n m : LatticeCharge) : Phase :=
  2 * (n 1 : Phase) * (m 0 : Phase) +
  2 * (n 2 : Phase) * ((m 0 : Phase) + m 1) +
  2 * (n 3 : Phase) * ((m 0 : Phase) + m 1 + m 2) +
  2 * (n 4 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3) +
  2 * (n 5 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3 + m 4) +
  2 * (n 6 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3 + m 4 + m 5) +
  2 * (n 7 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3 + m 4 + m 5 + m 6)

theorem eps_add_left (n m l : LatticeCharge) :
    epsExp (n + m) l = epsExp n l + epsExp m l := by
  simp [epsExp, Pi.add_apply, Int.cast_add]
  ring

theorem eps_add_right (n m l : LatticeCharge) :
    epsExp n (m + l) = epsExp n m + epsExp n l := by
  simp [epsExp, Pi.add_apply, Int.cast_add]
  ring

theorem eps_cocycle (n m l : LatticeCharge) :
    epsExp n m + epsExp (n + m) l =
      epsExp m l + epsExp n (m + l) := by
  rw [eps_add_left, eps_add_right]
  ring

theorem eps_neg_swap (n : LatticeCharge) :
    epsExp (-n) n = epsExp n (-n) := by
  simp [epsExp]
  ring

structure LatticeBasisWord where
  phase : Phase
  charge : LatticeCharge

@[ext] theorem latticeBasisWord_ext {x y : LatticeBasisWord}
    (hphase : x.phase = y.phase) (hcharge : x.charge = y.charge) : x = y := by
  cases x
  cases y
  simp_all

def basisMul (x y : LatticeBasisWord) : LatticeBasisWord :=
  ⟨x.phase + y.phase + epsExp x.charge y.charge, x.charge + y.charge⟩

def basisOne : LatticeBasisWord := ⟨0, 0⟩

/-- The inverse/star of i^phase E_n for the actual ordered cocycle. -/
def basisStar (x : LatticeBasisWord) : LatticeBasisWord :=
  ⟨-x.phase - epsExp (-x.charge) x.charge, -x.charge⟩

theorem basisMul_assoc (x y z : LatticeBasisWord) :
    basisMul (basisMul x y) z = basisMul x (basisMul y z) := by
  rcases x with ⟨px, nx⟩
  rcases y with ⟨py, ny⟩
  rcases z with ⟨pz, nz⟩
  apply latticeBasisWord_ext
  · simp only [basisMul]
    rw [eps_add_left, eps_add_right]
    ring
  · simp [basisMul]
    abel

theorem basis_unitary_left (x : LatticeBasisWord) :
    basisMul (basisStar x) x = basisOne := by
  rcases x with ⟨p, n⟩
  apply latticeBasisWord_ext
  · simp [basisMul, basisStar, basisOne]
    ring
  · simp [basisMul, basisStar, basisOne]

theorem basis_unitary_right (x : LatticeBasisWord) :
    basisMul x (basisStar x) = basisOne := by
  rcases x with ⟨p, n⟩
  apply latticeBasisWord_ext
  · simp [basisMul, basisStar, basisOne, eps_neg_swap]
    ring
  · simp [basisMul, basisStar, basisOne]

/-- Fixed algebraic data for a two-graded crossed word algebra.  The relation
    represented by this data is X a X^-1 = phi(a), X^2 = u. -/
structure CrossedData (A : Type*) [Ring A] [StarRing A] where
  phi : A ≃+* A
  u : A
  uInv : A
  uInv_mul_u : uInv * u = 1
  u_mul_uInv : u * uInv = 1
  phi_sq : ∀ a : A, phi (phi a) = u * a * uInv
  phi_u : phi u = u
  phi_star : ∀ a : A, phi (star a) = star (phi a)
  u_unitary_left : star u * u = 1
  u_unitary_right : u * star u = 1

variable {A : Type*} [Ring A] [StarRing A]

/-- The unique two-component normal form a0 + a1 X. -/
abbrev CrossedWord (A : Type*) := A × A

/-- Multiplication obtained by moving every X to the right. -/
def crossedMul (d : CrossedData A) (x y : CrossedWord A) : CrossedWord A :=
  (x.1 * y.1 + x.2 * d.phi y.2 * d.u,
   x.1 * y.2 + x.2 * d.phi y.1)

/-- Star obtained from (a X)^* = X^-1 a^* and X^-1 = u^-1 X. -/
def crossedStar (d : CrossedData A) (x : CrossedWord A) : CrossedWord A :=
  (star x.1, d.uInv * d.phi (star x.2))

/-- Projection onto the degree-zero coefficient. -/
def E0 (x : CrossedWord A) : A := x.1

def crossedOne : CrossedWord A := (1, 0)

def crossedX : CrossedWord A := (0, 1)

theorem uInv_eq_star_u (d : CrossedData A) : d.uInv = star d.u := by
  calc
    d.uInv = d.uInv * 1 := by simp
    _ = d.uInv * (d.u * star d.u) := by rw [d.u_unitary_right]
    _ = (d.uInv * d.u) * star d.u := by noncomm_ring
    _ = star d.u := by rw [d.uInv_mul_u]; simp

theorem star_uInv (d : CrossedData A) : star d.uInv = d.u := by
  rw [uInv_eq_star_u, star_star]

theorem uInv_mul_u_assoc (d : CrossedData A) (a : A) :
    d.uInv * (d.u * a) = a := by
  rw [← mul_assoc, d.uInv_mul_u, one_mul]

theorem u_mul_uInv_assoc (d : CrossedData A) (a : A) :
    d.u * (d.uInv * a) = a := by
  rw [← mul_assoc, d.u_mul_uInv, one_mul]

theorem uInv_mul_u_assoc₂ (d : CrossedData A) (a b : A) :
    d.uInv * (d.u * a * b) = a * b := by
  calc
    d.uInv * (d.u * a * b) = (d.uInv * d.u) * a * b := by noncomm_ring
    _ = a * b := by rw [d.uInv_mul_u]; simp

theorem phi_uInv (d : CrossedData A) : d.phi d.uInv = d.uInv := by
  calc
    d.phi d.uInv = d.phi d.uInv * 1 := by simp
    _ = d.phi d.uInv * (d.u * d.uInv) := by rw [d.u_mul_uInv]
    _ = (d.phi d.uInv * d.phi d.u) * d.uInv := by
      rw [d.phi_u]
      noncomm_ring
    _ = d.phi (d.uInv * d.u) * d.uInv := by rw [d.phi.map_mul]
    _ = d.uInv := by rw [d.uInv_mul_u, map_one]; simp

theorem phi_sq_mul_u (d : CrossedData A) (a : A) :
    d.phi (d.phi a) * d.u = d.u * a := by
  rw [d.phi_sq]
  calc
    (d.u * a * d.uInv) * d.u = d.u * a * (d.uInv * d.u) := by noncomm_ring
    _ = d.u * a := by rw [d.uInv_mul_u]; simp

/-- The adjoined defect generator X is unitary in the displayed star algebra. -/
theorem crossedX_unitary_left (d : CrossedData A) :
    crossedMul d (crossedStar d crossedX) crossedX = crossedOne := by
  apply Prod.ext
  · simp [crossedMul, crossedStar, crossedX, crossedOne, d.uInv_mul_u]
  · simp [crossedMul, crossedStar, crossedX, crossedOne]

theorem crossedX_unitary_right (d : CrossedData A) :
    crossedMul d crossedX (crossedStar d crossedX) = crossedOne := by
  apply Prod.ext
  · simp [crossedMul, crossedStar, crossedX, crossedOne, phi_uInv, d.uInv_mul_u]
  · simp [crossedMul, crossedStar, crossedX, crossedOne]

/-- The concrete normal-form multiplication is associative. -/
theorem crossedMul_assoc (d : CrossedData A) (x y z : CrossedWord A) :
    crossedMul d (crossedMul d x y) z = crossedMul d x (crossedMul d y z) := by
  rcases x with ⟨a, b⟩
  rcases y with ⟨c, e⟩
  rcases z with ⟨f, g⟩
  apply Prod.ext
  · simp only [crossedMul, d.phi.map_add, d.phi.map_mul, d.phi_u]
    simp only [d.phi_sq]
    noncomm_ring [d.uInv_mul_u, d.u_mul_uInv]
  · simp only [crossedMul, d.phi.map_add, d.phi.map_mul, d.phi_u]
    simp only [d.phi_sq]
    noncomm_ring [d.uInv_mul_u, d.u_mul_uInv]

/-- The displayed star really is involutive on the normal forms. -/
theorem crossedStar_involutive (d : CrossedData A) (x : CrossedWord A) :
    crossedStar d (crossedStar d x) = x := by
  rcases x with ⟨a, b⟩
  apply Prod.ext
  · simp [crossedStar]
  · simp only [crossedStar, star_mul, d.phi_star, star_star, star_uInv,
      d.phi.map_mul, d.phi_sq, d.phi_u]
    simp [mul_assoc, d.uInv_mul_u, d.u_mul_uInv, uInv_mul_u_assoc]

/-- The displayed star reverses the concrete crossed multiplication. -/
theorem crossedStar_mul (d : CrossedData A) (x y : CrossedWord A) :
    crossedStar d (crossedMul d x y) =
      crossedMul d (crossedStar d y) (crossedStar d x) := by
  rcases x with ⟨a, b⟩
  rcases y with ⟨c, e⟩
  apply Prod.ext
  · simp only [crossedStar, crossedMul, star_add, star_mul, d.phi_star,
      star_uInv, d.phi.map_mul, phi_uInv, d.phi_sq]
    rw [← uInv_eq_star_u d]
    simp [mul_assoc, d.uInv_mul_u, d.u_mul_uInv, uInv_mul_u_assoc]
  · simp only [crossedStar, crossedMul, star_add, star_mul, d.phi_star,
      star_uInv, d.phi.map_add, d.phi.map_mul, phi_uInv, d.phi_sq, d.phi_u]
    rw [← uInv_eq_star_u d]
    rw [mul_add, uInv_mul_u_assoc₂]
    noncomm_ring

theorem phi_inv_apply (d : CrossedData A) (a : A) :
    d.phi.symm a = d.uInv * d.phi a * d.u := by
  have h := d.phi_sq (d.phi.symm a)
  rw [d.phi.apply_symm_apply] at h
  calc
    d.phi.symm a = 1 * d.phi.symm a * 1 := by simp
    _ = (d.uInv * d.u) * d.phi.symm a * (d.uInv * d.u) := by
      rw [d.uInv_mul_u]
    _ = d.uInv * (d.u * d.phi.symm a * d.uInv) * d.u := by noncomm_ring
    _ = d.uInv * d.phi a * d.u := by rw [← h]

/-- Exact finite-word SOS identity.  No positivity or state is assumed. -/
theorem expectation_star_mul (d : CrossedData A) (y : CrossedWord A) :
    E0 (crossedMul d (crossedStar d y) y) =
      star y.1 * y.1 + d.phi.symm (star y.2 * y.2) := by
  rw [phi_inv_apply]
  simp only [E0, crossedMul, crossedStar]
  rw [d.phi.map_mul, d.phi_star]
  noncomm_ring

/-- The minimal algebraic interface needed from a positive base functional. -/
structure PositiveFunctional (A : Type*) [Ring A] [StarRing A] where
  toFun : A → ℝ
  map_add : ∀ a b : A, toFun (a + b) = toFun a + toFun b
  nonneg_star_mul_self : ∀ a : A, 0 ≤ toFun (star a * a)

/-- Positivity of omega0 o E0 follows from the base positivity alone.  No
    dual-Z2 invariance assumption is used. -/
theorem induced_positive (d : CrossedData A) (omega : PositiveFunctional A)
    (y : CrossedWord A) :
    0 ≤ omega.toFun (E0 (crossedMul d (crossedStar d y) y)) := by
  rw [expectation_star_mul, omega.map_add]
  have h0 := omega.nonneg_star_mul_self y.1
  have h1 := omega.nonneg_star_mul_self (d.phi.symm y.2)
  have hstar : star (d.phi.symm y.2) = d.phi.symm (star y.2) := by
    apply d.phi.injective
    rw [d.phi.apply_symm_apply, d.phi_star, d.phi.apply_symm_apply]
  rw [hstar, ← d.phi.symm.map_mul] at h1
  exact add_nonneg h0 h1

/-! Dual-Z2 symmetry and the exact scope of uniqueness. -/

def dualAction (x : CrossedWord A) : CrossedWord A := (x.1, -x.2)

theorem dualAction_involutive (x : CrossedWord A) :
    dualAction (dualAction x) = x := by
  rcases x with ⟨a, b⟩
  simp [dualAction]

theorem dualAction_mul (d : CrossedData A) (x y : CrossedWord A) :
    dualAction (crossedMul d x y) =
      crossedMul d (dualAction x) (dualAction y) := by
  rcases x with ⟨a, b⟩
  rcases y with ⟨c, e⟩
  apply Prod.ext
  · simp [dualAction, crossedMul]
  · simp [dualAction, crossedMul]
    abel

theorem dualAction_star (d : CrossedData A) (x : CrossedWord A) :
    dualAction (crossedStar d x) = crossedStar d (dualAction x) := by
  rcases x with ⟨a, b⟩
  apply Prod.ext <;> simp [dualAction, crossedStar, map_neg] <;> noncomm_ring

theorem E0_dual_invariant (x : CrossedWord A) : E0 (dualAction x) = E0 x := rfl

/-- Among additive real-valued extensions of omega0, dual-Z2 invariance kills
    the odd coefficient and therefore singles out omega0 o E0. -/
theorem dual_invariant_extension_unique
    (omega0 : A → ℝ) (Omega : CrossedWord A → ℝ)
    (hadd : ∀ x y, Omega (x + y) = Omega x + Omega y)
    (hext : ∀ a, Omega (a, 0) = omega0 a)
    (hinv : ∀ x, Omega (dualAction x) = Omega x)
    (x : CrossedWord A) : Omega x = omega0 (E0 x) := by
  rcases x with ⟨a, b⟩
  have hzero : Omega (0, 0) = 0 := by
    have h := hadd (0, 0) (0, 0)
    have hp : ((0, 0) + (0, 0) : CrossedWord A) = (0, 0) := by
      apply Prod.ext <;> simp
    rw [hp] at h
    linarith
  have hoddneg : Omega (0, -b) = Omega (0, b) := by
    simpa [dualAction] using hinv (0, b)
  have hodd : Omega (0, b) = 0 := by
    have hsum := hadd (0, b) (0, -b)
    have hp : ((0, b) + (0, -b) : CrossedWord A) = (0, 0) := by
      apply Prod.ext <;> simp
    rw [hp] at hsum
    rw [hzero, hoddneg] at hsum
    linarith
  have hsplit := hadd (a, 0) (0, b)
  have hp : ((a, 0) + (0, b) : CrossedWord A) = (a, b) := by
    apply Prod.ext <;> simp
  rw [hp] at hsplit
  rw [hext, hodd, add_zero] at hsplit
  exact hsplit

/-! Time invariance is conditional and separate from C11 positivity. -/

/-- Any grade-preserving evolution whose even component is `time0` commutes
    with the degree-zero expectation. -/
theorem E0_intertwines_grade_preserving_time
    (time0 : A → A) (time : CrossedWord A → CrossedWord A)
    (hgrade : ∀ x, (time x).1 = time0 x.1) (x : CrossedWord A) :
    E0 (time x) = time0 (E0 x) := hgrade x

/-- Consequently the induced state is time invariant only when the base state
    is time invariant and the evolution preserves the grading. -/
theorem induced_time_invariant
    (omega0 : A → ℝ) (time0 : A → A) (time : CrossedWord A → CrossedWord A)
    (hgrade : ∀ x, (time x).1 = time0 x.1)
    (hbase : ∀ a, omega0 (time0 a) = omega0 a)
    (x : CrossedWord A) :
    omega0 (E0 (time x)) = omega0 (E0 x) := by
  rw [E0_intertwines_grade_preserving_time time0 time hgrade, hbase]

#check crossedMul_assoc
#check crossedStar_involutive
#check crossedStar_mul
#check crossedX_unitary_left
#check crossedX_unitary_right
#check basisMul_assoc
#check basis_unitary_left
#check basis_unitary_right
#check expectation_star_mul
#check induced_positive
#check dual_invariant_extension_unique
#check induced_time_invariant
#print axioms crossedMul_assoc
#print axioms crossedStar_involutive
#print axioms crossedStar_mul
#print axioms crossedX_unitary_left
#print axioms crossedX_unitary_right
#print axioms basisMul_assoc
#print axioms basis_unitary_left
#print axioms basis_unitary_right
#print axioms expectation_star_mul
#print axioms induced_positive
#print axioms dual_invariant_extension_unique
#print axioms induced_time_invariant

end SeamStateExpectation
