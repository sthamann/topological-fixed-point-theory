"""Original quarter folding, native clock charges and line projection of v1000 TT."""
import ast
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as ss
import sympy as sp

from v1_e8_glue import e8_roots

repo = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[3]
contract = repo/"experiments/theory-contracts/chat-triangulation-20261108"
spec = importlib.util.spec_from_file_location("quarter",contract/"check_quarter_reference.py")
quarter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quarter)


def folding():
    G, DY, DM = quarter.grading, quarter.higgs_matrix, quarter.mass_matrix
    gamma = ss.csr_matrix(G)
    majorana = ss.csr_matrix(-1j*G@DM)
    yuk = ss.csr_matrix(-1j*G@DY)
    modes = np.arange(-8,8)
    order = sorted(range(len(modes)), key=lambda i:(modes[i]%8,modes[i]//8))
    new = modes[order]
    k, n = new%8,new//8
    one = ss.eye(len(modes),format="csr")
    edges = ss.diags([np.ones(len(modes)-1)]*2,[-1,1],format="csr")
    # Small smooth external profile in the actual original Higgs direction.
    amplitude = .23
    radius = 1
    perm = np.array([i*len(G)+j for i in order for j in range(len(G))])
    candidate = ss.csr_matrix(amplitude/2*(abs(8*(n[:,None]-n)+k[:,None]-k)==1))
    inter = (k[:,None]!=k)
    assert np.all(inter[candidate.nonzero()])
    rows=[]
    for alpha in (.25,.75):
        signed = (modes+alpha)/radius
        folded_signed = (8/radius)*(n+(k+alpha)/8)
        assert np.array_equal(signed[order],folded_signed)
        assert np.array_equal(signed[order]<0,n<0)
        assert np.array_equal(signed[order]>0,n>=0)
        original = ss.kron(ss.diags(signed),gamma)+ss.kron(one,majorana)+ss.kron(amplitude/2*edges,yuk)
        folded = ss.kron(ss.diags(folded_signed),gamma)+ss.kron(one,majorana)+ss.kron(candidate,yuk)
        residual = np.max(np.abs((original.tocsr()[perm,:][:,perm]-folded).data),initial=0)
        assert residual == 0
        rows.append({"quarter":alpha,"original_folded_matrix_max_error":float(residual),
                     "cross_residue_mode_entries":int(candidate.nnz),
                     "yukawa_coefficient_rank":int(np.linalg.matrix_rank(yuk.toarray()))})
    return rows


def clock():
    M=sp.Matrix(json.loads((repo/"experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json").read_text())["explicit_Weyl_matrix"])
    tree=ast.parse((repo/"verification/v634_st31_structure.py").read_text())
    idx=next(ast.literal_eval(node.value) for node in tree.body if isinstance(node,ast.Assign)
             and any(isinstance(target,ast.Name) and target.id=="SIGMA_IDX" for target in node.targets))
    S=sp.eye(8)[idx,:]
    R=M.T*S*M
    p=sp.eye(8)[:,0]+sp.eye(8)[:,2]
    orbit=[R**j*p for j in range(3)]

    def sector(v):
        if all(x.q==1 for x in v):
            return "0" if int(sum(v))%2==0 else "v"
        return "s" if int(sum(v-sp.ones(8,1)/2))%2==0 else "c"

    assert [sector(q) for q in orbit]==["0","s","0"]
    # The actual previously tested 0->2->0 glue word in fixed target coordinates.
    beta=sp.eye(8)[:,0]+sp.eye(8)[:,3]
    previous_word=[M.T*S**j*beta for j in range(3)]
    word_types=[sector(q) for q in previous_word]
    assert word_types==["s","s","s"]
    # Find a witness in the same fixed class orbit which starts and ends in
    # the even CAR lattice, without choosing a different clock or frame.
    def grade(v):
        d=v[:5]
        if all(x.q==1 for x in d):
            return 2*(int(sum(d))%2)
        return 1+2*(int(sum(x-sp.Rational(1,2) for x in d))%2)

    counts={}
    witness=None
    for v in e8_roots():
        v=sp.Matrix([sp.Rational(float(x)) for x in v])
        path=[S**j*v for j in range(3)]
        if [grade(x) for x in path] != [0,2,0]:
            continue
        kinds=tuple(sector(M.T*x) for x in path)
        counts[kinds]=counts.get(kinds,0)+1
        if kinds==("0","s","0") and witness is None:
            witness=v
    assert witness is not None
    assert counts=={("s","s","s"):4,("0","s","0"):4,("0","0","s"):4}
    energies=[]
    for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
        delta=sp.Matrix([(i+alpha)/8-sp.Rational(1,2) for i in range(8)])
        K=lambda q:q.dot(q)/2+delta.dot(q)
        assert K(R*p)-K(p)==sp.Rational(3,4)
        multiples={str(n):str(K(n*R*p)-K(n*p)) for n in (1,2,8,32)}
        assert multiples=={"1":"3/4","2":"3/2","8":"6","32":"24"}
        target=[S**j*witness for j in range(3)]
        path=[M.T*x for x in target]
        phase1=sp.I**(grade(target[1])-grade(target[0]))
        phase2=sp.I**(grade(target[2])-grade(target[0]))
        assert (phase1,phase2)==(-1,1)
        assert sector(path[1]) not in ("0","v") and sector(path[0])==sector(path[2])=="0"
        energies.append({"quarter":str(alpha),"original_clock_delta_K_by_charge_multiple":multiples,
                         "even_spinor_even_witness_energy":[str(K(q)) for q in path],
                         "witness_word_decay":str(K(path[1])+K(path[2])-2*K(path[0])),
                         "early_and_late_phase_factors":[int(phase1),int(phase2)],
                         "word_with_bare_CAR_projection_at_middle":0})
    return {"R100_CAR_sectors":[sector(q) for q in orbit],"previous_word_CAR_sectors":word_types,
            "fixed_glue_020_paths_by_CAR_sectors":{" ".join(k):v for k,v in counts.items()},
            "even_spinor_even_target_root":list(map(str,witness)), "cases":energies}


def line_spectral():
    # v1000 unit-spacing half-link symbol, k_parallel=pi/3 hence q_parallel=1.
    # The two sources below are Frobenius unit symmetric matrices E11 and E23.
    # Physical linearized tidal readout in the free TT vacuum has weight omega^3/8.
    def sum_band(size,eps):
        k=2*np.pi*np.fft.fftfreq(size)
        q=2*np.sin(k/2)
        # Only these transverse modes can enter a threshold band of width eps.
        q=q[abs(q)<=np.sqrt(2*eps+eps**2)]
        qy,qz=np.meshgrid(q,q)
        s=1+qy*qy+qz*qz
        w=np.sqrt(s)
        inside=(w<=1+eps)
        wt=w**3/(8*size**2)
        g11=.5*((s-1)/s)**2
        # ||P_TT E23||^2 = theta22*theta33 for the symmetric off-diagonal unit.
        g23=(1-qy*qy/s)*(1-qz*qz/s)
        return float((wt[inside]*g11[inside]).sum()),float((wt[inside]*g23[inside]).sum())
    n=3072
    rows=[]
    for eps in (.04,.02,.01,.005):
        a,b=sum_band(n,eps)
        rows.append([eps,a/eps**3,b/eps])
    target=[1/(24*np.pi),1/(16*np.pi)]
    np.testing.assert_allclose(rows[-1][1:],target,rtol=.025)
    atoms=[{"transverse_size":s,"tangential_atom":0,
            "transverse_tidal_atom":1/(8*s*s)} for s in (96,192,384)]
    return {"threshold_band_ratios":rows,"analytic_limits":target,"finite_threshold_atoms":atoms}


print(json.dumps({"original_R032_mode_fold":folding(),"native_clock":clock(),
                  "v1000_fixed_line_TT":line_spectral()},indent=2))
