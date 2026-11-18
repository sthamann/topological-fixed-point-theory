"""Conditional audit of the explicitly stated equations on PDF page 17.

This tests that given pulse, not unavailable external encoding matrices or
the microscopic TFPT origin. No repository source is changed or promoted.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.diag([1., -1.]).astype(complex)
I4 = np.eye(4, dtype=complex)
G = [np.kron(X, I2), np.kron(Y, I2), np.kron(Z, X),
     np.kron(Z, Y), np.kron(Z, Z)]

def require(cond, msg):
    if not cond:
        raise RuntimeError(msg)

for a in range(5):
    for b in range(5):
        require(np.max(np.abs(G[a]@G[b]+G[b]@G[a] -
                (2*I4 if a == b else np.zeros((4,4))))) < 1e-14,
                'Clifford identity')

# System-first ordering; |v>=|0>, |a>=|a+1>.
T = np.zeros((24, 24), complex)
for a, ga in enumerate(G):
    e = np.zeros((6,6)); e[a+1,0] = 1
    T += np.kron(ga,e)/np.sqrt(5)
P, Q = T.conj().T@T, T@T.conj().T
H = P+Q+T+T.conj().T
require(np.max(np.abs(H@H-2*H)) < 1e-14, 'H^2=2H')
require(np.linalg.eigvalsh(H).min() > -1e-14, 'H positive')

def M(x):
    return sum(ga@x@ga for ga in G)

def target(x):
    return 7*x/12+M(x)/12

def predicted(x, theta, p0):
    w = (1-p0)/5
    s = np.sin(theta)**2
    return ((1-s*(p0+2*w))*x+(p0+w)*s*M(x)/5 +
            w*s*M(M(x))/25)

def direct(x, theta, p0):
    env = np.diag([p0]+[(1-p0)/5]*5)
    u = expm(-1j*theta*H)
    out = u@np.kron(x,env)@u.conj().T
    return np.trace(out.reshape(4,6,4,6),axis1=1,axis2=3)

matrix_units = []
for a in range(4):
    for b in range(4):
        e = np.zeros((4,4),complex); e[a,b]=1
        matrix_units.append(e)

max_formula_error = max(np.max(np.abs(direct(e,t,p)-predicted(e,t,p)))
    for e in matrix_units for t in [0.,0.31,0.81,1.3]
    for p in [1.,7/12,1/6])
theta_native = np.arcsin(np.sqrt(5/12))
pure_target_error = max(np.max(np.abs(direct(e,theta_native,1)-target(e)))
    for e in matrix_units)
require(max_formula_error < 2e-14,'all 16 matrix units formula')
require(pure_target_error < 2e-14,'pure source reproduces native channel')

p0,w,s = sp.symbols('p0 w s',real=True)
lv = 1-s*(sp.Rational(8,5)*p0+sp.Rational(56,25)*w)
lb = 1-s*(sp.Rational(4,5)*p0+sp.Rational(44,25)*w)
ratio_obstruction = sp.factor((1-lv)-2*(1-lb))
require(ratio_obstruction == -sp.Rational(32,25)*s*w,
        'exact target-ratio obstruction')
thermal_lv = sp.simplify(lv.subs({p0:sp.Rational(7,12),w:sp.Rational(1,12)}))
thermal_lb = sp.simplify(lb.subs({p0:sp.Rational(7,12),w:sp.Rational(1,12)}))
sv = sp.solve(sp.Eq(thermal_lv,sp.Rational(1,3)),s)[0]
sb = sp.solve(sp.Eq(thermal_lb,sp.Rational(2,3)),s)[0]

# Minimax of absolute errors of the two contraction eigenvalues on s in [0,1].
best_s = sp.solve(sp.Eq(thermal_lv-sp.Rational(1,3),
               -(thermal_lb-sp.Rational(2,3))),s)[0]
best_err = sp.simplify(thermal_lv.subs(s,best_s)-sp.Rational(1,3))

results={
 'scope':'Given 24D H_sym on page 17 and Spin(5)-invariant diagonal bath; pair H switched off during pulse.',
 'source_pdf':'/Users/stefanhamann/Projekte/tfpt-theoryv4/TFPT-Kernessenz-Deutsch.pdf',
 'source_page':17,
 'verdict':'NO_EXACT_NATIVE_CHANNEL_FOR_NONPURE_INVARIANT_BATH_WITH_THIS_PULSE',
 'all_matrix_units_checked':16,
 'max_formula_error':float(max_formula_error),
 'pure_source_native_error':float(pure_target_error),
 'vector_contraction':str(lv),
 'bivector_contraction':str(lb),
 'target_vector':'1/3','target_bivector':'2/3',
 'obstruction':'(1-lambda_vector)-2*(1-lambda_bivector)=-32*s*w/25',
 'required_bath':'w=0, p0=1, s=5/12',
 'thermal_state':'p0=7/12, each quintet weight w=1/12',
 'thermal_vector':str(thermal_lv),'thermal_bivector':str(thermal_lb),
 'required_s_vector':str(sv),'required_s_bivector':str(sb),
 'minimax_s':str(best_s),'minimax_contraction_error':str(best_err),
 'minimax_error_float':float(best_err),
 'native_pulse_thermal_vector':str(thermal_lv.subs(s,sp.Rational(5,12))),
 'native_pulse_thermal_bivector':str(thermal_lb.subs(s,sp.Rational(5,12))),
 'not_derived':['a different pulse is impossible','all thermal source realizations fail',
     'physical TFPT source selection','external encoding/64x64 matrix identity'],
}

# Exact Pauli-channel diamond distance: a Bell input attains the l1 difference
# of 16 orthogonal unitary probabilities; triangle inequality is the upper bound.
thermal_diamond = (abs(sp.Rational(5,12)-sp.Rational(11,15)*s) +
                  5*abs(sp.Rational(2,15)*s-sp.Rational(1,12)) + s/15)
optimal_s = sp.Rational(25,44)
results.update({
 'thermal_optimal_s_for_diamond':str(optimal_s),
 'thermal_exact_optimal_unhalved_diamond':str(thermal_diamond.subs(s,optimal_s)),
 'thermal_exact_optimal_half_diamond':str(thermal_diamond.subs(s,optimal_s)/2),
 'thermal_native_pulse_half_diamond':str(thermal_diamond.subs(s,sp.Rational(5,12))/2),
 'thermal_cooling_variable':'r=exp(-8*beta*g)',
 'thermal_cooling_optimal_s':'25*(1+5*r)/(12*(5+9*r))',
 'thermal_cooling_optimal_half_diamond':'5*r/(3*(5+9*r))',
 'temperature_for_half_diamond_epsilon':'beta*g >= log((5-27*epsilon)/(15*epsilon))/8 for 0<epsilon<5/42',
})

# Positive construction retaining the pair-preparation Hamiltonian during pulse.
# In the antisymmetric bath 1+5 its energies are 0 and 8*g. Starting from |v>,
# only P+Q is reached: the effective bright pair is [[k,k],[k,k+8g]].
g, k = 1., 4.
omega = np.sqrt(k*k+16*g*g)
t_match = np.arcsin(np.sqrt((5/12)*(1+16*g*g/(k*k))))/omega
pair_bath = np.diag([0.]+[8*g]*5)
H_total = k*H+np.kron(I4,pair_bath)
U_total = expm(-1j*t_match*H_total)
env0 = np.diag([1.]+[0.]*5)
pair_on_error = 0.
for e in matrix_units:
    out = U_total@np.kron(e,env0)@U_total.conj().T
    reduced = np.trace(out.reshape(4,6,4,6),axis1=1,axis2=3)
    pair_on_error = max(pair_on_error,float(np.max(np.abs(reduced-target(e)))))
require(pair_on_error < 2e-14,'pair-on construction all16')
results.update({
 'pair_on_population':'kappa^2/(kappa^2+16*g^2) * sin(t*sqrt(kappa^2+16*g^2))^2',
 'pair_on_exact_match_iff':'kappa/g >= sqrt(80/7), for positive g,kappa and pure singlet bath',
 'pair_on_matching_time':'asin(sqrt((5/12)*(1+16*g^2/kappa^2)))/sqrt(kappa^2+16*g^2)',
 'pair_on_example_kappa_over_g':4.,
 'pair_on_example_g_times_t':float(g*t_match),
 'pair_on_max_matrix_unit_error':pair_on_error,
 'pair_on_scope':'Only the two displayed H_pair/H_sym equations and six-dimensional antisymmetric bath; pure ground prep and kappa/g choice remain inputs.',
})
Path(__file__).with_name('joint_bath_pulse_results.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
