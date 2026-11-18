"""Decide the actual character-valued lattice lift by finite GF2 equations."""
from collections import deque
from pathlib import Path
import itertools as it
import json
import pickle
import numpy as np
import sympy as S

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'equivariant-monodromy-20261023/imported/work'

def key(x):
    return tuple(map(int,x[0].flat))+tuple(map(int,x[1].flat))

def mul(x,y):
    a,b=x;c,d=y
    r,s=a@c-b@d,a@d+b@c
    if np.any(r%2) or np.any(s%2): raise RuntimeError('denominator')
    return r//2,s//2

def real_twice(x):
    a,b=x;R=np.zeros((8,8),dtype=np.int64)
    R[::2,::2]=a;R[::2,1::2]=-b;R[1::2,::2]=b;R[1::2,1::2]=a
    return R

def act(cols,x):
    a=0
    for i in range(8):
        if (x>>i)&1:a^=cols[i]
    return a

def bilinear(C,x,y):
    return sum(int(C[i,j]) for i,j in it.product(range(8),repeat=2)
               if (x>>i)&1 and (y>>j)&1)%2

def phase_table(M,C):
    D=(M.T@C@M+C)%2
    if np.any(np.diag(D)) or np.any(D!=D.T):raise RuntimeError('quadratic refinement')
    return [sum(int(D[i,j]) for i in range(8) for j in range(i+1,8)
                if (x>>i)&1 and (x>>j)&1)%2 for x in range(256)]

def solve(gens,B,C,stop_first=True):
    n=len(gens);Id=(2*np.eye(4,dtype=np.int64),np.zeros((4,4),dtype=np.int64))
    Bi=np.array((2*S.Matrix(B.tolist()).inv()).tolist(),dtype=np.int64)
    def columns(g):
        N=Bi@real_twice(g)@B
        if np.any(N%4):raise RuntimeError('not integral lattice isometry')
        M=(N//4)%2
        return [sum(int(M[i,j])<<i for i in range(8)) for j in range(8)]
    gs=[]
    for g in gens:
        cols=columns(g)
        M=np.array([[(cols[j]>>i)&1 for j in range(8)] for i in range(8)],dtype=np.int64)
        gs.append((cols,phase_table(M,C)))
    # At g, each coefficient row is a linear expression in generator basis signs.
    data={key(Id):(Id,[1<<i for i in range(8)],[0]*8,0,[])}
    queue=deque([key(Id)]);basis={};pivots=[];equations=0
    while queue:
        gk=queue.popleft();g,gcols,forms,const,word=data[gk]
        for si,s in enumerate(gens):
            y=mul(s,g);yk=key(y)
            yforms=[forms[i] ^ (gcols[i] << (8*si)) for i in range(8)]
            yc=const ^ sum(gs[si][1][gcols[i]]<<i for i in range(8))
            yw=[si]+word
            if yk not in data:
                cols=[act(gs[si][0],x) for x in gcols]
                data[yk]=(y,cols,yforms,yc,yw);queue.append(yk)
            else:
                _,_,oldforms,oldc,oldword=data[yk]
                for i in range(8):
                    co=yforms[i]^oldforms[i];rhs=((yc^oldc)>>i)&1
                    equations+=1;row=co;value=rhs;cert=0
                    while row:
                        p=row.bit_length()-1
                        if p not in basis:break
                        z,t,sup=basis[p];row^=z;value^=t;cert^=sup
                    desc={'left_word':yw,'right_word':oldword,'coordinate':i,
                          'coefficients':co,'rhs':rhs}
                    if row:
                        ordinal=len(pivots);pivots.append(desc)
                        basis[row.bit_length()-1]=(row,value,cert^(1<<ordinal))
                    elif value:
                        support=[pivots[j] for j in range(len(pivots)) if (cert>>j)&1]+[desc]
                        out={'split':False,'generators':n,'group_elements_reached':len(data),
                             'rank':len(basis),'equations_processed':equations,
                             'inconsistency_certificate':support}
                        return out
    # Solve pivot system with free variables zero; all solutions are affine.
    sol=0
    for p in sorted(basis):
        row,value,_=basis[p]
        if ((row&sol).bit_count()&1)!=value:sol^=1<<p
    return {'split':True,'generators':n,'group_order':len(data),'rank':len(basis),
            'dimension_of_splittings':8*n-len(basis),'one_solution':sol,
            'equations_processed':equations}

def main():
    refs,perms,pg,chosen,orbits,full,orbit_id,rays=pickle.loads((OLD/'state.pkl').read_bytes())
    old=json.loads((OLD/'results.json').read_text());B=np.array(old['integral_lattice_basis'],dtype=np.int64)
    Gram=np.array(old['integral_gram'],dtype=np.int64)
    C=np.tril(Gram,-1)%2
    for i in range(8):C[i,i]=(Gram[i,i]//2)%2
    results=[]
    for k in range(1,len(chosen)+1):
        out=solve([refs[j] for j in chosen[:k]],B,C)
        out['generator_reflection_indices']=chosen[:k];results.append(out)
        print(json.dumps({a:b for a,b in out.items() if a!='inconsistency_certificate'}),flush=True)
        if not out['split']:break
    (HERE/'probe-results.json').write_text(json.dumps({'Gram':Gram.tolist(),'C':C.tolist(),'prefix_results':results},indent=2)+'\n')

if __name__=='__main__':main()
