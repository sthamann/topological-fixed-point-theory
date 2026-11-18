import Mathlib

namespace TFPT.WeakSingletCut

def cap : Matrix (Fin 3) (Fin 3) ℚ := fun _ _ => 1/3

def relativeTurn : Matrix (Fin 3) (Fin 3) ℚ := Matrix.diagonal ![-1, -1, 1]

theorem capProjection : cap * cap = cap := by
  ext i j
  norm_num [cap, Matrix.mul_apply, Fin.sum_univ_three]

theorem relativeTurnInvolution : relativeTurn * relativeTurn = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [relativeTurn, cap, Matrix.mul_apply, Fin.sum_univ_three, Matrix.diagonal_apply,
      Matrix.one_apply, Matrix.smul_apply, Matrix.cons_val_two, Matrix.vecHead, Matrix.vecTail,
      Matrix.cons_val_succ, Matrix.cons_val_zero]

theorem oneLegCompressed : cap * relativeTurn * cap = (-1/3 : ℚ) • cap := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [relativeTurn, cap, Matrix.mul_apply, Fin.sum_univ_three, Matrix.diagonal_apply,
      Matrix.one_apply, Matrix.smul_apply, Matrix.cons_val_two, Matrix.vecHead, Matrix.vecTail,
      Matrix.cons_val_succ, Matrix.cons_val_zero]

theorem twoTurnsNoIntermediateCap :
    cap * (relativeTurn * relativeTurn) * cap = cap := by
  rw [relativeTurnInvolution, mul_one, capProjection]

theorem twoTurnsWithIntermediateCap :
    cap * relativeTurn * cap * relativeTurn * cap = (1/9 : ℚ) • cap := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [relativeTurn, cap, Matrix.mul_apply, Fin.sum_univ_three, Matrix.diagonal_apply,
      Matrix.one_apply, Matrix.smul_apply, Matrix.cons_val_two, Matrix.vecHead, Matrix.vecTail,
      Matrix.cons_val_succ, Matrix.cons_val_zero]

#print axioms oneLegCompressed
#print axioms twoTurnsNoIntermediateCap
#print axioms twoTurnsWithIntermediateCap

end TFPT.WeakSingletCut
