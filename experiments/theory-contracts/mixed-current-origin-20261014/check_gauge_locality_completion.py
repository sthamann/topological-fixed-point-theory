"""Exact control of a CONDITIONAL locality change, not its P1 derivation.

Question: can a fermion-parity Gauss completion make the existing half flux
endpoint-local while preserving the full even physical Hilbert space and
the original register Hamiltonian? The proof is all-size; this small exact
CAR control checks its signs with the actual v1033 hopping coefficients.
The first new premise is the Gauss/local-observable rule, not a new energy.
"""
from pathlib import Path
import hashlib
import importlib
import json
import sys
import numpy as np
import sympy as s

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(ROOT/'verification'))
v1033=importlib.import_module('v1033_charged_disorder')


def need(value,label):
    if not value:raise ValueError(label)


def zero(matrix):
    return all(s.expand(x)==0 for x in matrix.todok().values())


def annihilator(m,k):
    out=s.MutableSparseMatrix(2**m,2**m,{})
    for bits in range(2**m):
        if bits&(1<<k):
            out[bits^(1<<k),bits]=(-1)**((bits&((1<<k)-1)).bit_count())
    return out


def exact_matrix(np_matrix):
    return s.SparseMatrix([[s.nsimplify(x) for x in row] for row in np_matrix])


def second_quantize(h,ops):
    out=s.zeros(ops[0].rows)
    for j in range(h.rows):
        for k in range(h.cols):
            if h[j,k]:out+=h[j,k]*ops[j].T*ops[k]
    return s.SparseMatrix(out)


# Three columns, one transverse site, both actual QWZ orbital components.
# This is solely an exact group/CAR/sign control, not the physical width8
# source or a new source model. The analytic isometry works for arbitrary
# column Hilbert spaces, onsite terms and original hopping coefficients.
L,ny=(int(sys.argv[1]) if len(sys.argv)>1 else 3),1
m=2*L*ny
dim=2**m
ops=[annihilator(m,k) for k in range(m)]
even=[n for n in range(dim) if n.bit_count()%2==0]
states=[(r,n) for r in range(4) for n in even]
pos={a:i for i,a in enumerate(states)}


def col_parity(n,j):
    mask=((1<<(2*ny))-1)<<(2*ny*j)
    return (n&mask).bit_count()%2


def ext_index(b,z,n):
    return (b*2**L+z)*dim+n


extended=2*2**L*dim
E=s.MutableSparseMatrix(extended,len(states),{})
for (r,n),col in pos.items():
    b,w=r%2,r//2
    for tbits in range(2**(L-1)):
        t=[0]+[(tbits>>(j-1))&1 for j in range(1,L)]
        z=sum((((w if j==L-1 else 0)+t[j]+t[(j+1)%L])%2)<<j for j in range(L))
        phase=(-1)**sum(t[j]*col_parity(n,j) for j in range(1,L))
        E[ext_index(b,z,n),col]=phase/s.sqrt(2**(L-1))
E=s.SparseMatrix(E)
need(E.T*E==s.eye(len(states)),'isometry on FULL even register-matter carrier')

gauss=[]
electric=[]
for j in range(L):
    G=s.MutableSparseMatrix(extended,extended,{})
    X=s.MutableSparseMatrix(extended,extended,{})
    for b in range(2):
        for z in range(2**L):
            for n in range(dim):
                col=ext_index(b,z,n)
                G[ext_index(b,z^(1<<j)^(1<<((j-1)%L)),n),col]=(-1)**col_parity(n,j)
                X[ext_index(b,z^(1<<j),n),col]=1
    G,X=s.SparseMatrix(G),s.SparseMatrix(X)
    need(G*E==E,'all Gauss constraints, including anchor, on even carrier')
    gauss.append(G);electric.append(X)

# Assemble the exact original dGamma(h_r) from the current builder.
Hblocks=[second_quantize(exact_matrix(v1033.qwz_cylinder(L,ny,1.,r)),ops) for r in range(4)]
H_old=s.MutableSparseMatrix(len(states),len(states),{})
for r in range(4):
    for n in even:
        for nprime in even:
            value=Hblocks[r][nprime,n]
            if value:H_old[pos[r,nprime],pos[r,n]]=value
H_old=s.SparseMatrix(H_old)

# The gauged hopping is obtained from the SAME coefficients: link signs z_j
# multiply the forward QWZ hopping at every j; i^b remains at the seam.
H_new=s.MutableSparseMatrix(extended,extended,{})
TX=exact_matrix(v1033.TX)
for b in range(2):
    for z in range(2**L):
        h=exact_matrix(v1033.qwz_cylinder(L,ny,1.,b))
        for x in range(L):
            amp=(s.I**b if x==L-1 else 1)
            change=((-1)**((z>>x)&1)-1)*amp
            if change:
                for a in range(2):
                    for c in range(2):
                        i,j=2*((x+1)%L)+a,2*x+c
                        h[i,j]+=change*TX[a,c]
                        h[j,i]+=s.conjugate(change*TX[a,c])
        block=second_quantize(h,ops)
        for n in range(dim):
            for nprime in range(dim):
                if block[nprime,n]:H_new[ext_index(b,z,nprime),ext_index(b,z,n)]=block[nprime,n]
H_new=s.SparseMatrix(H_new)
need(zero(H_new*E-E*H_old),'EXACT original physical Hamiltonian and all time words')
need(all(H_new*G==G*H_new for G in gauss),'local Gauss constraints commute with Hamiltonian')

W=s.diag(*[s.I**b*(-1)**z.bit_count() for b in range(2) for z in range(2**L) for n in range(dim)])
Zold=s.diag(*[s.I**r for r,n in states])
need(W*E==E*Zold,'old local register Z becomes global Wilson loop')

for endpoint in range(L):
    F=s.MutableSparseMatrix(len(states),len(states),{})
    for (r,n),col in pos.items():
        parity=sum(col_parity(n,x) for x in range(endpoint+1))%2
        F[pos[(r+2)%4,n],col]=(-1)**parity
    need(electric[endpoint]*E==E*F,'existing half string becomes LOCAL electric endpoint')

# Any nontrivial register-grade q continues to obstruct old local inclusion.
X4=s.zeros(4)
for r in range(4):X4[(r+1)%4,r]=1
Z4=s.diag(1,s.I,-1,-s.I)
need(all(Z4*X4**q*Z4.conjugate().T==s.I**q*X4**q for q in range(4)),
     'all register grades survive arbitrary matter dressing and homogeneous words')

# Independently check the one-particle gauge covariance at the ACTUAL width8.
# All values are binary rationals times powers of i, so these identities are
# exact in floating representation, not tolerance-based spectral checks.
def signed_cylinder(nx,ny,b,z):
    h=v1033.qwz_cylinder(nx,ny,1.,b).copy()
    for x in range(nx):
        if (z>>x)&1:
            change=-2*(1j**b if x==nx-1 else 1)
            for y in range(ny):
                i=2*((x+1)%nx*ny+y)
                j=2*(x*ny+y)
                h[i:i+2,j:j+2]+=change*v1033.TX
                h[j:j+2,i:i+2]+=np.conj(change)*v1033.TX.conj().T
    return h

width8_cases=0
for b in range(2):
    for z in range(8):
        h=signed_cylinder(3,8,b,z)
        for x in range(3):
            u=np.ones(48)
            u[16*x:16*(x+1)]=-1
            znew=z^(1<<x)^(1<<((x-1)%3))
            need(np.array_equal(u[:,None]*h*u[None,:],signed_cylinder(3,8,b,znew)),
                 'actual width8 QWZ gauge covariance with SAME coefficients')
            width8_cases+=1

result={
    'date':'2026-10-22',
    'verdict':'EXACT_GAUSS_LOCALITY_COMPLETION_CONDITIONAL_NOT_P1_DERIVED',
    'control_columns':L,'control_transverse_sites':ny,
    'original_even_physical_dimension':len(states),
    'extended_redundant_dimension':extended,
    'isometry_and_all_Gauss_constraints_exact':True,
    'same_original_QWZ_hopping_coefficients_and_physical_H_exact':True,
    'half_flux_maps_to_local_electric_link_operator_exact':True,
    'old_register_Z_maps_to_GLOBAL_Wilson_loop_exact':True,
    'actual_width8_one_body_gauge_covariance_cases':width8_cases,
    'no_new_physical_dimension_or_energy_coefficient':True,
    'first_new_assumption':'Local Gauss constraints and gauge-invariant local observable assignment; not supplied by original P1.',
    'full_Dirac_chiral_measure_selected':False,
    'second_edge_removed':False,
    'E8_spinor_origin_derived':False,
    'source_sha256':hashlib.sha256((ROOT/'verification/v1033_charged_disorder.py').read_bytes()).hexdigest(),
}
Path(__file__).with_name('gauge-locality-completion-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
