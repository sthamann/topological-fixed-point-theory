"""Check the Fuchs address closure and the original two-sheet sewing tests."""
import json
import sys
from pathlib import Path

import sympy as sp

repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "verification"))
from v117_monodromy_weyl_a3 import M0_EXACT as M, U_EXACT as U
from v120_address_table import L_LEP

I3 = sp.eye(3)
Z3 = sp.zeros(3)
I6 = sp.eye(6)
Z6 = sp.zeros(6)
w = (-1+sp.I*sp.sqrt(3))/2
Cs = sp.Matrix([[0,0,1],[-1,0,0],[0,-1,0]])
G = sp.Matrix([[0,0,-1+sp.I],[1,sp.I,0],[sp.I,1,0]])/sp.sqrt(2)
Ts = Z3.row_join(Cs).col_join(I3.row_join(Z3))
J = I3.row_join(I3).col_join(Cs.row_join(-Cs))/sp.sqrt(2)
raw_sign = sp.diag(I3,-I3)
sign_swap = Z3.row_join(I3).col_join(I3.row_join(Z3))
assert sp.simplify(J.H * raw_sign * J) == sign_swap
assert sp.simplify(J.H * Ts**3 * J) == raw_sign
assert sp.simplify(raw_sign*Ts*raw_sign + Ts) == Z6
assert sp.simplify(G*Cs*G.H-M) == Z3
M6, U6 = sp.diag(M,-M),sp.diag(U,U)
assert sp.simplify(M6**3) == raw_sign


def line(z):
    return sp.simplify(sum((sp.conjugate(z)**j*M**j for j in range(3)),Z3)/3)


expected = (w,-w,-1)
for length, phase in zip(L_LEP.values(), expected):
    assert sp.simplify(sp.expand_complex(sp.exp(sp.I*sp.pi*(length%6)/3))-phase) == 0

P3 = sp.diag(line(w),line(w)+line(1))
Q3 = I6-P3
assert sp.simplify(P3*P3-P3) == Z6 and sp.trace(P3)==3
assert sp.simplify(P3*M6-M6*P3) == Z6
assert sp.trace(P3*raw_sign)==-1

P4 = sp.diag(line(w)+line(1),line(w)+line(1))
swap_loss = sp.simplify(P3*sign_swap*Q3*sign_swap*P3)
assert sp.trace(swap_loss)==1 and swap_loss.rank()==1
assert sp.simplify(swap_loss-sp.diag(Z3,line(1))) == Z6
assert sp.simplify(P4*sign_swap-sign_swap*P4) == Z6
assert sp.simplify(P4*P3-P3)==Z6 and sp.trace(P4)==4

# Exact composition defect with the actual quarter marker, not OS time.
leakage = sp.simplify(P3*U6.H*Q3*U6*P3)
two_word_defect = sp.simplify(P3*U6.H*U6*P3-(P3*U6.H*P3)*(P3*U6*P3))
assert sp.simplify(leakage-two_word_defect)==Z6
assert leakage.rank()==2

def closure(P, maps):
    columns = P.columnspace()
    ranks = [len(columns)]
    while True:
        whole = sp.Matrix.hstack(*columns)
        new = sp.simplify(sp.Matrix.hstack(whole, *(T*whole for T in maps))).columnspace()
        if len(new)==len(columns):
            return ranks
        ranks.append(len(new))
        columns=new

ranks = closure(P3,[M6,U6,U6.H])
assert ranks[-1]==6
trace_per_lepton = [sp.simplify(sp.trace(line(w)*leakage[:3,:3])),
                    sp.simplify(sp.trace(line(w)*leakage[3:,3:])),
                    sp.simplify(sp.trace(line(1)*leakage[3:,3:]))]

# Actual C6 signed kernel. Both sheets are held in the same marked metric.
y,t = sp.symbols("y t",real=True)
D = y*I6-t*M6
opposite_product = sp.simplify((y*I3-t*M)*(y*I3+t*M))
assert sp.simplify(opposite_product-(y*y*I3-t*t*M**2))==Z3
assert sp.simplify(sign_swap*D*sign_swap-D.subs(t,-t))==Z6
det_y = sp.factor(D.det())
assert sp.simplify(det_y-(y**6-t**6))==0
# y=1 and the actual charged-family t=1/2; compare amplitudes, not masses.
paired_return = sp.sqrt(((1-sp.Rational(1,2)*w)*(1-sp.Rational(1,2)*sp.conjugate(w)))
                        *((1+sp.Rational(1,2)*w)*(1+sp.Rational(1,2)*sp.conjugate(w))))
paired_return = sp.simplify(1/paired_return)
assert paired_return==4*sp.sqrt(21)/21
assert sp.simplify(D.subs({y:1,t:sp.Rational(1,2)}).H*D.subs({y:1,t:sp.Rational(1,2)})-I6)!=Z6

# The already established v487 six-step pair can be iterated, but its
# connected contraction bound is strictly gapped at a fixed time step.
R=sp.Matrix([[1,0,0],[0,sp.Rational(1,2),sp.Rational(1,6)],[0,sp.Rational(1,6),sp.Rational(1,2)]])
rates=sorted(R.__pow__(6).eigenvals())
assert rates==[sp.Rational(1,729),sp.Rational(64,729),sp.Integer(1)]
# The pure normalization can obey exact sewing too. A scalar Hamiltonian
# shift multiplies each interval, and survives the adjoint reflected norm.
t1,t2,E=sp.symbols("t1 t2 E",real=True)
assert sp.simplify(sp.exp(-E*t1)*sp.exp(-E*t2)-sp.exp(-E*(t1+t2)))==0
assert sp.simplify(sp.conjugate(sp.exp(-E*t1))*sp.exp(-E*t1))==sp.exp(-2*E*t1)

# v81 determinant cover in the v82 coordinate q=-3x.
z, q = sp.symbols("z q")
mu = sp.Rational(64,729)
lift_amplitude = sp.Rational(8,27)
rho = (q-2)/(5-q)
def base_of(z):
    return (2-5*z*z)/(1-z*z)


def fiber_of(z):
    return 3*z/(1-z*z)


assert sp.cancel(fiber_of(z)**2-(base_of(z)-2)*(base_of(z)-5))==0
assert sp.cancel(rho.subs(q,base_of(z))+z*z)==0
assert sp.cancel(base_of(-z)-base_of(z))==0
assert sp.cancel(fiber_of(-z)+fiber_of(z))==0
assert sp.cancel(base_of(-1/z)-(7-base_of(z)))==0
assert sp.cancel(fiber_of(-1/z)-fiber_of(z))==0
F = (2+5*mu*rho)/(1+mu*rho)
assert sp.cancel(F.subs(q,base_of(z))-base_of(lift_amplitude*z))==0
assert sp.cancel(base_of(lift_amplitude*z)-base_of(-lift_amplitude*z))==0
assert sp.cancel(rho.subs(q,7-q)*rho-1)==0
assert sp.cancel((-1/(-lift_amplitude/z))-z/lift_amplitude)==0
assert lift_amplitude**2==rates[1] and lift_amplitude**2!=rates[0]

result={
 "signed_path_parity": {"half_turn_trace_on_three_addresses":sp.trace(P3*raw_sign),
                        "fiber_sign_in_signed_frame":"exchanges the two signs"},
 "direct_sign_exchange": {"three_addresses":3,"first_invariant_space":4,
                          "loss_rank":swap_loss.rank(),"missing_address":"(+,1), partner of tau (-,1)"},
 "original_quarter_marker": {"reachable_ranks":ranks,"return_defect_rank":leakage.rank(),
                             "defect_trace":sp.trace(leakage),"loss_on_e_mu_tau":trace_per_lepton,
                             "nonzero_defect_eigenvalues":{str(x):n for x,n in leakage.eigenvals().items() if x!=0}},
 "signed_c6": {"determinant":det_y,"forward_and_opposite_green_absolute_product_at_omega":paired_return,
               "opposite_is_inverse":False},
 "v81_v82_lifts":{"coordinate":"z=Y/(5-q), z^2=-rho(q)",
                  "deck":"z -> -z, base unchanged",
                  "branch_exchange":"z -> -1/z, q -> 7-q",
                  "two_lifts":["(8/27) z","-(8/27) z"],
                  "branch_exchange_conjugates_lift_to_inverse":True,
                  "squared_multiplier":"64/729",
                  "faster_v487_mode_not_obtained_by_deck_sign":"1/729"},
 "v487_fixed_step_spectrum":rates,
 "continuous_scalar_sewing_ambiguity":"exp(-t E(g)), reflected factor exp(-2 t E(g))",
}
print(json.dumps(result,indent=2,default=str))
