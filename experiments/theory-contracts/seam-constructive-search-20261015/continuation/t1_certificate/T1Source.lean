import Mathlib

namespace OriginalWinding

abbrev Charge := Fin 8 → ℤ
abbrev Phase := ZMod 4

def cycle (n : Charge) : Charge := ![n 7, n 0, n 1, n 2, n 3, n 4, n 5, n 6]
def halfSwap (n : Charge) : Charge := ![n 4, n 5, n 6, n 7, n 0, n 1, n 2, n 3]
def q : Charge := fun _ => 1
def eps (n m : Charge) : Phase :=
  2 * (n 1 : Phase) * m 0 +
  2 * (n 2 : Phase) * ((m 0 : Phase) + m 1) +
  2 * (n 3 : Phase) * ((m 0 : Phase) + m 1 + m 2) +
  2 * (n 4 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3) +
  2 * (n 5 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3 + m 4) +
  2 * (n 6 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3 + m 4 + m 5) +
  2 * (n 7 : Phase) * ((m 0 : Phase) + m 1 + m 2 + m 3 + m 4 + m 5 + m 6)
def gamma (n : Charge) : Phase :=
  2 * (n 7 : Phase) * ((n 0 : Phase) + n 1 + n 2 + n 3 + n 4 + n 5 + n 6)
def eta (n : Charge) : Phase :=
  2 * ((n 0 : Phase) + n 1 + n 2 + n 3) * ((n 4 : Phase) + n 5 + n 6 + n 7)

/-- The ordered CAR cocycle is carried by the actual one-channel cycle. -/
theorem cycle_cocycle (n m : Charge) :
    gamma n + gamma m + eps (cycle n) (cycle m) = eps n m + gamma (n + m) := by
  simp [gamma, eps, cycle, Pi.add_apply, Int.cast_add]
  ring_nf
  simp only [show (4 : Phase) = 0 by decide, mul_zero, add_zero]

theorem cycle_four (n : Charge) : cycle (cycle (cycle (cycle n))) = halfSwap n := by
  funext j
  fin_cases j <;> simp [cycle, halfSwap]

/-- Four canonical cycle lifts induce precisely the original half-swap sign. -/
theorem cycle_four_phase (n : Charge) :
    gamma n + gamma (cycle n) + gamma (cycle (cycle n)) +
      gamma (cycle (cycle (cycle n))) = eta n := by
  simp [gamma, cycle, eta]
  ring_nf
  simp only [show (4 : Phase) = 0 by decide, mul_zero, add_zero]

theorem gamma_q : gamma q = 2 := by decide

/-- Under the ADDITIONAL scalar eigenoperator law C X C^-1 = nu X,
    X^2 proportional to E_q and the source cycle sign give nu^2=-1.
    General sector-dependent lifts need not satisfy this scalar law.
    This implication does not establish existence or physical selection of X. -/
theorem square_minus_one_fourth_one (nu : ℂ) (h : nu ^ 2 = -1) : nu ^ 4 = 1 := by
  calc
    nu ^ 4 = (nu ^ 2) ^ 2 := by ring
    _ = 1 := by rw [h]; norm_num

theorem same_fourth_lift_selects_plus (nu kappa : ℂ)
    (square : nu ^ 2 = -1) (sameLift : kappa = nu ^ 4) : kappa = 1 := by
  rw [sameLift, square_minus_one_fourth_one nu square]

abbrev QCharge := Fin 8 → ℚ
def cycleQ (p : QCharge) : QCharge := ![p 7, p 0, p 1, p 2, p 3, p 4, p 5, p 6]
def e0 : QCharge := ![1, 0, 0, 0, 0, 0, 0, 0]
def charge (p : QCharge) : ℚ := p 0 + p 1 + p 2 + p 3 + p 4 + p 5 + p 6 + p 7
def normSq (p : QCharge) : ℚ :=
  p 0 ^ 2 + p 1 ^ 2 + p 2 ^ 2 + p 3 ^ 2 + p 4 ^ 2 + p 5 ^ 2 + p 6 ^ 2 + p 7 ^ 2
def deltaDot (alpha : ℚ) (p : QCharge) : ℚ :=
  (alpha / 8 - 1 / 2) * p 0 + ((1 + alpha) / 8 - 1 / 2) * p 1 +
  ((2 + alpha) / 8 - 1 / 2) * p 2 + ((3 + alpha) / 8 - 1 / 2) * p 3 +
  ((4 + alpha) / 8 - 1 / 2) * p 4 + ((5 + alpha) / 8 - 1 / 2) * p 5 +
  ((6 + alpha) / 8 - 1 / 2) * p 6 + ((7 + alpha) / 8 - 1 / 2) * p 7
def energy (alpha N : ℚ) (p : QCharge) : ℚ := N + normSq p / 2 + deltaDot alpha p

theorem single_step_charge (p : QCharge) : charge (cycleQ p + e0) - charge p = 1 := by
  simp [charge, cycleQ, e0, Pi.add_apply]
  ring

/-- The original fixed quadratic time is retained for every rational charge,
    including both the integer and half-lattice sectors. -/
theorem single_step_time (alpha N : ℚ) (p : QCharge) :
    energy alpha N (cycleQ p + e0) - energy alpha N p = (charge p + alpha) / 8 := by
  simp [energy, normSq, deltaDot, charge, cycleQ, e0, Pi.add_apply]
  ring

#print axioms cycle_cocycle
#print axioms cycle_four_phase
#print axioms same_fourth_lift_selects_plus
#print axioms single_step_time

end OriginalWinding
