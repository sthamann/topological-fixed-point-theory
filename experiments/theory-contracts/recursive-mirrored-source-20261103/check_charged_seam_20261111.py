"""Exact charged field checks in the original CAR chart and finite v252 carrier."""
from fractions import Fraction as F
import importlib.util
import json
from itertools import product
from pathlib import Path
import sys

import numpy as np


def add(a, b, factor=F(1)):
    out = a.copy()
    for key, value in b.items():
        out[key] = out.get(key, F(0)) + factor*value
    return {k: v for k, v in out.items() if v}


def field(v, r, create):
    i = r+8
    assert 0 <= i < 16
    out = {}
    for mask, value in v.items():
        if bool(mask & (1<<i)) != create:
            sign = (-1)**((mask & ((1<<i)-1)).bit_count())
            out[mask^(1<<i)] = sign*value
    return out


def vacuum(q):
    return {(1<<(8+q))-1: F(1)}


def energy(mask, a):
    return sum((F(i-8)+a)*((mask>>i & 1)-(i<8)) for i in range(16))


def mode(v, n, old, alpha=None):
    if n == 0:
        return {mask: value*e for mask, value in v.items()
                if (e := mask.bit_count()-8 if alpha is None else energy(mask, alpha))}
    out = {}
    for mask, coefficient in v.items():
        for i in range(max(0,n),min(16,16+n)):
            move = old.fermion_move(mask,i,i-n)
            if move is not None:
                target, sign = move
                weight = F(1) if alpha is None else F(i-8)+alpha-F(n,2)
                out[target] = out.get(target,F(0))+sign*coefficient*weight
    return {k:v for k,v in out.items() if v}


def charged(repo):
    entry = (repo/"experiments/theory-contracts/original-admissible-source-closure-20261025"
             /"check_w8_stress_cyclic.py")
    spec = importlib.util.spec_from_file_location("original_stress",entry)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    def original_quarter_current(v,n):
        upper = old.original_stress_action(v,n,F(3,4))
        lower = old.original_stress_action(v,n,F(1,4))
        difference = {k:2*c for k,c in add(upper,lower,-1).items()}
        assert difference == mode(v,-n,old)
        return difference

    coefficients, wards, flows = 0, 0, 0
    for q in range(-3,4):
        psi_coefficients = [field(vacuum(q),q,True)]
        for degree in range(5):
            if degree:
                out = {}
                for k in range(1,degree+1):
                    out = add(out,original_quarter_current(psi_coefficients[degree-k],k),F(1,degree))
                psi_coefficients.append(out)
            assert psi_coefficients[degree] == field(vacuum(q),q+degree,True)
            coefficients += 1
        for offset in (1,2):
            v = mode(vacuum(q),-offset,old)
            # Charge spectral flow on the entire tested descendant, not only the ground.
            shift = lambda a: {(mask<<1)|1: value for mask,value in a.items()}
            assert all((mask>>15) == 0 for mask in v)
            assert shift(v) == mode(vacuum(q+1),-offset,old)
            for alpha in (F(1,4),F(3,4)):
                for mask in v:
                    new = (mask<<1)|1
                    assert energy(new,alpha)-energy(mask,alpha) == q+alpha
                    neutral_energy = energy(mask,F(1,2))
                    assert (neutral_energy-F(q,2)).denominator == 1
                    assert (energy(mask,alpha)-alpha*q).denominator == 1
                    flows += 1
            for n,r,create in product(range(-2,3),range(-4,5),(True,False)):
                fv = field(v,r,create)
                r_out = r-n if create else r+n
                sign = 1 if create else -1
                expected = field(v,r_out,create)
                assert add(mode(fv,n,old),field(mode(v,n,old),r,create),-1) == add({},expected,sign)
                for alpha in (F(1,4),F(3,4)):
                    lfv = mode(fv,n,old,alpha)
                    flv = field(mode(v,n,old,alpha),r,create)
                    assert add(lfv,flv,-1) == add({},expected,sign*(r+alpha-sign*F(n,2)))
                    # Subtract the actual quarter shift using that same current.
                    corrected = add(add(lfv,flv,-1),expected,-sign*(alpha-F(1,2)))
                    assert corrected == add({},expected,sign*(r+F(1,2)-sign*F(n,2)))
                    wards += 1
    return {"charged_vertex_coefficients":coefficients,"original_and_corrected_stress_ward_checks":wards,
            "relative_energy_and_spin_checks":flows}


def finite_matter(repo):
    # The identical finite source matrices and Higgs chart used in QTheta.
    pkg = (repo/"experiments/theory-contracts/original-paper-source-reconciliation-20261102"
           /"contribution_20261106Q")
    sys.path.insert(0,str(pkg/"runtime"))
    from assemble import Assembler
    a = Assembler(pkg/"domain/native-diagnostic-v1.json")
    v = a.original
    zeros = np.zeros((3,3),complex)
    K = a.K
    gamma = a.gamma
    y_dict = {"nu":{"L":-3,"R":0},"e":{"L":-3,"R":-6},
              "u":{"L":1,"R":4},"d":{"L":1,"R":-2}}
    y = np.array([y_dict[f][side] for f,side,_,_ in v.P])
    w = np.array([{"nu":1,"e":-1,"u":1,"d":-1}[f] if side=="L" else 0 for f,side,_,_ in v.P])
    labels = v.P
    y = np.r_[y,-y]
    w = np.r_[w,-w]
    sample = [v.rep(lam,q,m) for lam,q,m in (
        (1,np.zeros((2,2)),zeros), (1j,np.zeros((2,2)),zeros),
        (0,np.eye(2),zeros), (0,v.quaternion(0,1),zeros),
        (0,v.quaternion(1j,0),zeros),(0,v.quaternion(0,1j),zeros))]
    for i,j in product(range(3),repeat=2):
        for coeff in (1,1j):
            m = zeros.copy()
            m[i,j] = coeff
            sample.append(v.rep(0,np.zeros((2,2)),m))
    right = [K@p.conj()@K for p in sample]
    proj = v.rep(1,np.zeros((2,2)),zeros)

    def inner(D,h):
        B = v.rep(0,v.quaternion(h[1].conjugate(),h[0]),zeros) @ (D@proj-proj@D)
        B = B+B.conj().T
        return B+K@B.conj()@K

    weak = (v.quaternion(0,1),v.quaternion(1j,0))
    color = (np.diag([1j,-1j,1]), np.roll(np.eye(3),1,axis=0))
    higgs_tests = [np.array(h,dtype=complex) for h in ((1,0),(1j,0),(0,1),(0,1j))]
    charged_checks = 0
    gauge_checks = 0
    for scale in (1,2):
        Y = list(a.Y)
        Y[2] = scale*Y[2]
        DY = v.full_D(*Y,zeros)
        T = -1j*gamma@(inner(DY,higgs_tests[0])-1j*inner(DY,higgs_tests[1]))/2
        for g,c,h in product((*weak,np.eye(2)),(*color,np.eye(3)),higgs_tests):
            left = v.rep(1,g,c)
            action = left@K@left.conj()@K
            np.testing.assert_array_equal(action@inner(DY,h)@action.conj().T,inner(DY,g@h))
            gauge_checks += 1
        # The quaternion convention selects the same positive charged arm.
        assert np.count_nonzero((y[:,None]-y[None,:]-3)*T)==0
        assert np.count_nonzero((w[:,None]-w[None,:]-1)*T)==0
        assert np.count_nonzero(gamma@T+T@gamma)==0
        D = v.full_D(*Y,a.MR)
        assert np.array_equal(D,D.conj().T)
        assert np.array_equal(K@D.conj()@K,D)
        assert np.count_nonzero(D@gamma+gamma@D)==0
        for L,R in product(sample,right):
            DL = D@L-L@D
            assert np.count_nonzero(DL@R-R@DL)==0
            charged_checks+=1
    before = float(np.linalg.norm(a.Y[2])/np.linalg.norm(a.Y[0]))
    return {"same_first_order_checks":charged_checks,"same_color_and_weak_covariance_checks":gauge_checks,"electron_over_up_yukawa_norm":[before,2*before]}


if __name__=="__main__":
    repo=Path(sys.argv[1]).resolve()
    print(json.dumps({"original_car":charged(repo),"finite_matter":finite_matter(repo)},indent=2))
