"""One M3 source -> traceless eight -> CAR16 -> candidate E8 phase counts.

Proves algebraic functorial bridge and a new NO-GO: simply lifting the
unique covariant M3 damping to the 128 even spinor states DOES NOT
supply a common weight-one E8 affine current grading.
"""
from itertools import combinations
from math import comb
from collections import Counter
import sympy as sp


def e(i,j):return sp.eye(3)[:,i]*sp.eye(3)[j,:]


def wedge_op(mask,i,creation):
    occupied=bool(mask & (1<<i))
    if occupied==creation:return None
    sign=(-1)**(int(mask&((1<<i)-1)).bit_count())
    return (mask^(1<<i),sign)


def apply_pair(mask,i,j,types):
    second=wedge_op(mask,j,types[1])
    if second is None:return {}
    first=wedge_op(second[0],i,types[0])
    if first is None:return {}
    return {first[0]:first[1]*second[1]}


def anti_comp(mask,i,j,types):
    a=apply_pair(mask,i,j,types)
    b=apply_pair(mask,j,i,types[::-1])
    for k,v in b.items():a[k]=a.get(k,0)+v
    return {k:v for k,v in a.items() if v}


def controls():
    herm=[sp.I*(e(0,1)-e(1,0)),sp.I*(e(0,2)-e(2,0)),
          sp.I*(e(1,2)-e(2,1)),e(0,1)+e(1,0),
          e(0,2)+e(2,0),e(1,2)+e(2,1),
          sp.diag(1,-1,0),sp.diag(1,1,-2)]
    metric=sp.Matrix(8,8,lambda i,j:sp.trace(herm[i].H*herm[j]))
    assert metric==sp.diag(*([2]*7+[6]))
    b=[x/sp.sqrt(metric[i,i]) for i,x in enumerate(herm)]
    assert sp.Matrix(8,8,lambda i,j:sp.simplify(sp.trace(b[i].H*b[j])))==sp.eye(8)
    H=sp.Matrix.hstack(*[v.reshape(9,1) for v in herm])
    assert H.rank()==8 and all(x.H==x and sp.trace(x)==0 for x in herm)
    rot=sp.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    a8=sp.Matrix.hstack(*[H.gauss_jordan_solve((rot*x*rot.T).reshape(9,1))[0]
                           for x in herm])
    spec=a8.eigenvals()
    assert spec=={sp.S.One:2,sp.S.NegativeOne:2,sp.I:2,-sp.I:2}
    assert a8**4==sp.eye(8)
    # Direct 8-complex-mode canonical anticommutation algebra on 256 masks.
    for mask in range(256):
      for i in range(8):
       for j in range(8):
        for types in ((False,False),(True,True),(False,True)):
         observed=anti_comp(mask,i,j,types)
         expected={mask:1} if types==(False,True) and i==j else {}
         if observed!=expected:
          raise ArithmeticError(('CAR',mask,i,j,types,observed,expected))
    degrees=(0,0,1,1,2,2,3,3)
    even=Counter(sum(degrees[i] for i in range(8) if mask & (1<<i))%4
                 for mask in range(256) if mask.bit_count()%2==0)
    vector=degrees+tuple((-d)%4 for d in degrees)
    bilinear=Counter((vector[i]+vector[j])%4
                     for i,j in combinations(range(16),2))
    assert [even[k] for k in range(4)]==[32]*4
    assert [bilinear[k] for k in range(4)]==[28,32,28,32]
    assert len([mask for mask in range(256) if mask.bit_count()%2==0])==128
    # E8=so16(120)+S+16(128) is a KNOWN abstract extension, NOT yet
    # affine E8 on this zero-dimensional Fock carrier.
    assert [even[k]+bilinear[k] for k in range(4)]==[60,64,60,64]
    # Second quantization of qI_8 has eigenvalue q^k on degree-k wedges.
    # Its vacuum exterior matrix element stays EXACTLY one.
    q=sp.Rational(1,2)
    spectral=Counter(mask.bit_count() for mask in range(256))
    assert spectral==Counter({k:comb(8,k) for k in range(9)})
    trace=sum(comb(8,k)*q**k for k in range(9))
    supertrace=sum((-1)**k*comb(8,k)*q**k for k in range(9))
    assert trace==(1+q)**8 and supertrace==(1-q)**8
    assert trace==sp.Rational(6561,256) and supertrace==sp.Rational(1,256)
    # Even spinor states live at k=0,2,4,6,8 for a lifted H=N.
    # No nonzero common affine weight-one assignment from H=kappa*N+c.
    weights=tuple(range(0,9,2))
    assert weights==(0,2,4,6,8)
    kappa,c=sp.symbols('kappa c')
    assert sp.solve((c-1,2*kappa+c-1),(kappa,c))=={kappa:0,c:1}
    return {
      'status':'EXACT_M3_TO_CAR16_SPINOR128_WITH_AFFINE_WEIGHT_ONE_NO_GO',
      'source_one_particle_hilbert':'traceless M3(C), Hilbert-Schmidt metric, dimension 8',
      'positive_metric_diagonal':([2]*7+[6]),
      'quarter_rotation_phase_counts':[2,2,2,2],
      'all_CAR_pairs_on_256_Fock_masks_verified':256*8*8*3,
      'complex_CAR_modes':8,
      'real_Majorana_generators':16,
      'full_Fock_dimension':256,
      'positive_even_spinor_dimension':128,
      'SO16_bilinear_dimension':120,
      'E8_abstract_extension_dimension':248,
      'E8_abstract_quarter_phase_counts':[60,64,60,64],
      'lifted_transfer':'Gamma(q I8) acts by q^k on wedge degree k',
      'lifted_transfer_trace_at_q_half':'6561/256',
      'lifted_transfer_supertrace_at_q_half':'1/256',
      'original_zero_Fock_cap_of_lift':'1 (no determinant)',
      'different_even_spinor_source_energies':[0,2,4,6,8],
      'source_fermion_number_dynamics_gives_all_E8_spinor_currents_weight_one':False,
      'nontrivial_common_weight_one_possible_by_affine_rescaling':False,
      'original_physical_CAR_or_E8_current_OPE_derived':False,
      'original_TFPT_Q0_and_3plus1_parent_derived':False,
      'F1_solved':False
    }


if __name__=='__main__':
 import json
 print(json.dumps(controls(),sort_keys=True,indent=2))
