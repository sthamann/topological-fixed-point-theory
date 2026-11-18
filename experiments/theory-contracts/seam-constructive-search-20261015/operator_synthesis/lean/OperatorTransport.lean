import Mathlib

/-!
# Exact algebraic transport for synthesized operators

This file proves the fixed, candidate-independent part of the operator-search
contract.  A solver must still check its concrete matrix equations.  The
theorems below say what follows from those equations: simultaneous generator
intertwining extends to the noncommutative ring they generate, unitary
conjugation preserves adjoints and positive squares, and an antilinear branch
must apply its stated conjugation automorphism to every source generator.

No source representation, physical state, analytic logarithm, or P1/P2
selection is constructed here.
-/

namespace TFPT.OperatorTransport

/-- `U` transports `a` to `b` when multiplication by `U` intertwines them. -/
def Intertwines {R : Type*} [Mul R] (U a b : R) : Prop := U * a = b * U

/-- A finite syntax for noncommutative integer polynomials in generators. -/
inductive RingExpr (ι : Type*) where
  | gen : ι → RingExpr ι
  | zero : RingExpr ι
  | one : RingExpr ι
  | add : RingExpr ι → RingExpr ι → RingExpr ι
  | neg : RingExpr ι → RingExpr ι
  | mul : RingExpr ι → RingExpr ι → RingExpr ι
  deriving Repr

namespace RingExpr

/-- Evaluate a noncommutative polynomial expression in a ring. -/
def eval {ι R : Type*} [Ring R] (ρ : ι → R) : RingExpr ι → R
  | gen i => ρ i
  | zero => 0
  | one => 1
  | add p q => eval ρ p + eval ρ q
  | neg p => -eval ρ p
  | mul p q => eval ρ p * eval ρ q

/-- Ring homomorphisms commute with evaluation of the fixed expression syntax. -/
theorem map_eval {ι R S : Type*} [Ring R] [Ring S]
    (σ : R →+* S) (ρ : ι → R) (p : RingExpr ι) :
    σ (eval ρ p) = eval (fun i => σ (ρ i)) p := by
  induction p with
  | gen i => rfl
  | zero => exact map_zero σ
  | one => exact map_one σ
  | add p q hp hq => simp only [eval, map_add, hp, hq]
  | neg p hp => simp only [eval, map_neg, hp]
  | mul p q hp hq => simp only [eval, map_mul, hp, hq]

end RingExpr

section RingTransport

variable {ι R : Type*} [Ring R]

/-- Intertwining is closed under sums. -/
theorem intertwines_add {U a b c d : R}
    (hab : Intertwines U a b) (hcd : Intertwines U c d) :
    Intertwines U (a + c) (b + d) := by
  calc
    U * (a + c) = U * a + U * c := mul_add _ _ _
    _ = b * U + d * U := by rw [hab, hcd]
    _ = (b + d) * U := (add_mul _ _ _).symm

/-- Intertwining is closed under additive inverses. -/
theorem intertwines_neg {U a b : R} (h : Intertwines U a b) :
    Intertwines U (-a) (-b) := by
  calc
    U * (-a) = -(U * a) := mul_neg _ _
    _ = -(b * U) := by rw [h]
    _ = (-b) * U := (neg_mul _ _).symm

/-- Intertwining is closed under products, with order retained. -/
theorem intertwines_mul {U a b c d : R}
    (hab : Intertwines U a b) (hcd : Intertwines U c d) :
    Intertwines U (a * c) (b * d) := by
  calc
    U * (a * c) = (U * a) * c := (mul_assoc _ _ _).symm
    _ = (b * U) * c := by rw [hab]
    _ = b * (U * c) := mul_assoc _ _ _
    _ = b * (d * U) := by rw [hcd]
    _ = (b * d) * U := (mul_assoc _ _ _).symm

/-- One checked generator relation transports every natural power. -/
theorem intertwines_pow {U a b : R} (h : Intertwines U a b) :
    ∀ n : ℕ, Intertwines U (a ^ n) (b ^ n) := by
  intro n
  induction n with
  | zero => simp [Intertwines]
  | succ n ih =>
      simpa [pow_succ] using intertwines_mul ih h

/-- Simultaneous generator relations transport every generated ring expression. -/
theorem intertwines_eval (U : R) (src dst : ι → R)
    (hgen : ∀ i, Intertwines U (src i) (dst i)) :
    ∀ p : RingExpr ι, Intertwines U (p.eval src) (p.eval dst) := by
  intro p
  induction p with
  | gen i => exact hgen i
  | zero => simp [Intertwines, RingExpr.eval]
  | one => simp [Intertwines, RingExpr.eval]
  | add p q hp hq => exact intertwines_add hp hq
  | neg p hp => exact intertwines_neg hp
  | mul p q hp hq => exact intertwines_mul hp hq

end RingTransport

section StarTransport

variable {R : Type*} [Ring R] [StarRing R]

/-- The two exact equations used from a unitary candidate. -/
structure UnitaryLike (U : R) : Prop where
  star_mul_self : star U * U = 1
  mul_star_self : U * star U = 1

/-- Conjugation by a unitary-like element. -/
def conjugateBy (U a : R) : R := U * a * star U

/-- An intertwining equation and right unitarity identify the target explicitly. -/
theorem intertwines_target_eq_conjugate {U a b : R}
    (hU : UnitaryLike U) (h : Intertwines U a b) :
    b = conjugateBy U a := by
  calc
    b = b * 1 := (mul_one b).symm
    _ = b * (U * star U) := by rw [hU.mul_star_self]
    _ = (b * U) * star U := (mul_assoc _ _ _).symm
    _ = (U * a) * star U := by rw [h]
    _ = conjugateBy U a := rfl

/-- Unitary intertwining transports adjoints as well as the checked generators. -/
theorem intertwines_star_of_unitary {U a b : R}
    (hU : UnitaryLike U) (h : Intertwines U a b) :
    Intertwines U (star a) (star b) := by
  have hs : star a * star U = star U * star b := by
    simpa only [star_mul] using congrArg star h
  calc
    U * star a = U * star a * 1 := (mul_one _).symm
    _ = U * star a * (star U * U) := by rw [hU.star_mul_self]
    _ = U * (star a * star U) * U := by simp only [mul_assoc]
    _ = U * (star U * star b) * U := by rw [hs]
    _ = (U * star U) * star b * U := by simp only [mul_assoc]
    _ = star b * U := by rw [hU.mul_star_self, one_mul]

/-- Unitary conjugation commutes with the adjoint operation. -/
theorem conjugateBy_star (U a : R) :
    conjugateBy U (star a) = star (conjugateBy U a) := by
  simp only [conjugateBy, star_mul, star_star, mul_assoc]

/-- Unitary conjugation transports an algebraic positive square exactly. -/
theorem conjugateBy_positive_square (U a : R) (hU : UnitaryLike U) :
    conjugateBy U (star a * a) =
      star (conjugateBy U a) * conjugateBy U a := by
  calc
    conjugateBy U (star a * a) =
        U * star a * (star U * U) * a * star U := by
          simp [conjugateBy, hU.star_mul_self, mul_assoc]
    _ = star (conjugateBy U a) * conjugateBy U a := by
          simp only [conjugateBy, star_mul, star_star, mul_assoc]

end StarTransport

section ConjugateLinearTransport

variable {ι R : Type*} [Ring R]

/-- Algebraic model of an antilinear branch: conjugate the source first. -/
def ConjugateIntertwines (σ : R ≃+* R) (U a b : R) : Prop :=
  U * σ a = b * U

/-- Generator checks with an explicit conjugation automorphism transport all
noncommutative polynomial expressions. -/
theorem conjugate_intertwines_eval (σ : R ≃+* R) (U : R)
    (src dst : ι → R)
    (hgen : ∀ i, ConjugateIntertwines σ U (src i) (dst i))
    (p : RingExpr ι) :
    U * σ (p.eval src) = p.eval dst * U := by
  change U * σ.toRingHom (p.eval src) = p.eval dst * U
  rw [RingExpr.map_eval σ.toRingHom src p]
  exact intertwines_eval U (fun i => σ (src i)) dst hgen p

end ConjugateLinearTransport

section ImaginarySign

variable {R : Type*} [Ring R] [StarRing R]

/-- Conjugation followed by the candidate unitary. -/
def conjugateLinearTransport (σ : R ≃+* R) (U a : R) : R :=
  conjugateBy U (σ a)

/-- A generator reversed by conjugation is transported with the opposite sign.
This is the exact guard against treating an antilinear candidate as linear on
an imaginary generator. -/
theorem imaginary_generator_changes_sign (σ : R ≃+* R) (U j : R)
    (hj : σ j = -j) :
    conjugateLinearTransport σ U j = -conjugateBy U j := by
  simp [conjugateLinearTransport, conjugateBy, hj]

end ImaginarySign

#print axioms RingExpr.map_eval
#print axioms intertwines_add
#print axioms intertwines_neg
#print axioms intertwines_mul
#print axioms intertwines_pow
#print axioms intertwines_eval
#print axioms intertwines_target_eq_conjugate
#print axioms intertwines_star_of_unitary
#print axioms conjugateBy_star
#print axioms conjugateBy_positive_square
#print axioms conjugate_intertwines_eval
#print axioms imaginary_generator_changes_sign

end TFPT.OperatorTransport
