import Mathlib

namespace TFPT.Det16SelectorCensus

def bit (w p : Nat) : Nat := w / (2^p) % 2
def color (w : Nat) : Nat := bit w 4 + bit w 3 + bit w 2
def weak (w : Nat) : Nat := bit w 1 + bit w 0
def charge6 (w : Nat) : Nat := (4 * color w + 3 * weak w) % 6
def splitCharge (w : Nat) : Nat := (color w + bit w 1) % 2
def centralCharge (w : Nat) : Nat := (2 * color w + weak w) % 4
def evenWedgeMasks : List Nat :=
  (List.range 32).filter fun w => (color w + weak w) % 2 == 0

def address (a b c : Nat) : Nat := (a * 2 + b) * 4 + c

abbrev Moments := Nat × Nat × Nat

def next (d : List Moments) (w : Nat) : List Moments :=
  (List.range 48).map fun k =>
    let a := k / 8
    let b := k / 4 % 2
    let c := k % 4
    let previous := address ((a + 6 - charge6 w) % 6)
      ((b + 2 - splitCharge w) % 2) ((c + 4 - centralCharge w) % 4)
    let (n, m, s) := d[k]!
    let (n', m', s') := d[previous]!
    (n+n', m+m'+n', s+s'+2*m'+n')

def census : List Moments :=
  evenWedgeMasks.foldl next ((1,0,0) :: List.replicate 47 (0,0,0))
def splitRank : Nat := ((List.range 4).map (fun c => (census[address 0 0 c]!).1)).sum
def centralRank : Nat := ((List.range 2).map (fun b => (census[address 0 b 0]!).1)).sum
def intersectionRank : Nat := (census[address 0 0 0]!).1

def centeredSquare (x : Moments) : Int :=
  (x.2.2 : Int) - 16 * (x.2.1 : Int) + 64 * (x.1 : Int)
def splitSquare : Int := ((List.range 4).map
  (fun c => centeredSquare (census[address 0 0 c]!))).sum
def centralSquare : Int := ((List.range 2).map
  (fun b => centeredSquare (census[address 0 b 0]!))).sum

set_option maxRecDepth 8192 in
set_option maxHeartbeats 0 in
theorem exactChargeCensus :
    (splitRank, centralRank, intersectionRank, centralSquare - splitSquare) =
    (5464,5464,3072,624) := by
  decide

noncomputable def exchange (phase : ℂ) : Matrix (Fin 2) (Fin 2) ℂ := !![0, star phase; phase, 0]
def densityContrast : Matrix (Fin 2) (Fin 2) ℂ := !![-1, 0; 0, 1]

theorem exchangeContrast (phase : ℂ) :
    exchange phase * densityContrast = -(densityContrast * exchange phase) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [exchange, densityContrast, Matrix.mul_apply, Fin.sum_univ_two]

#print axioms exactChargeCensus
#print axioms exchangeContrast
end TFPT.Det16SelectorCensus
