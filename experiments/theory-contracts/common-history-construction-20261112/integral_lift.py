from pathlib import Path
import importlib.util,itertools as it,json
import numpy as np
ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
HERE=Path(__file__).resolve().parent
p=ROOT/'experiments/theory-contracts/compiler-correlated-event-clock-20260919/source_channel.py'
spec=importlib.util.spec_from_file_location('native',p)
n=importlib.util.module_from_spec(spec);spec.loader.exec_module(n)
zs=n.source_rays();I=np.eye(4,dtype=complex)
rs=[I-np.outer(z,z.conj())/2 for z in zs]
def key(A):
    B=2*A
    assert np.array_equal(B.real,np.rint(B.real)) and np.array_equal(B.imag,np.rint(B.imag))
    return tuple(int(x) for z in B.ravel() for x in (z.real,z.imag))
def group(gens,bound):
    found={key(I):I};todo=[I]
    while todo:
        A=todo.pop()
        for g in gens:
            B=g@A;k=key(B)
            if k not in found:
                found[k]=B;todo.append(B)
                if len(found)>bound:return None
    return list(found.values())
gs=[(i,j,1j*rs[i]@rs[j]) for i,j in it.combinations(range(60),2) if np.vdot(zs[i],zs[j])==0]
hs=[(i,j,-rs[i]@rs[j]) for i,j in it.combinations(range(60),2) if abs(np.vdot(zs[i],zs[j]))**2==4]
answer=None
for i,j,g in gs:
    for k,l,h in hs:
        assert np.array_equal(g@g,-I) and np.array_equal(np.linalg.matrix_power(h,3),-I)
        product5=np.linalg.matrix_power(g@h,5)
        if np.array_equal(product5,-I) or np.array_equal(product5,I):
            if np.array_equal(product5,I):g=-g
            B=group([g,h],120)
            if B is not None and len(B)==120:
                G=group([g,h,1j*I],240)
                if G is not None and len(G)==240:
                    answer=(i,j,k,l,g,h,B,G);break
    if answer:break
assert answer is not None
i,j,k,l,g,h,B,G=answer
assert (i,j,k,l)==(0,46,1,21)
assert np.array_equal(g,-1j*rs[0]@rs[46])
assert np.array_equal(h,-rs[1]@rs[21])
roots=[phase*z for z in zs for phase in [1,1j,-1,-1j]]
def rkey(z):return tuple(z)
lookup={rkey(z):a for a,z in enumerate(roots)}
rootperms=[]
for A in G:
    assert np.array_equal(A.conj().T@A,I)
    perm=[lookup[rkey(A@z)] for z in roots]
    rootperms.append(perm)
remaining=set(range(240));orbits=[]
while remaining:
    a=min(remaining);orb={p[a] for p in rootperms}
    remaining-=orb;orbits.append(sorted(orb))
center=[A for A in G if np.array_equal(A@g,g@A) and np.array_equal(A@h,h@A)]
assert len(center)==4
# Binary tetrahedral subgroup as normalizer of an actual quaternion subgroup.
b=next(A for A in B if np.array_equal(A@A,-I) and np.array_equal(A@g,-g@A))
Q=group([g,b],8);assert Q is not None and len(Q)==8
Qkeys={key(A) for A in Q}
H=[A for A in G if all(key(A@x@A.conj().T) in Qkeys for x in Q)]
assert len(H)==48
Hkeys={key(A) for A in H}
cosets=[];reps=[];covered=set()
for A in G:
    if key(A) in covered:continue
    orbit={key(A@x) for x in H};cosets.append(orbit);reps.append(A);covered|=orbit
assert len(cosets)==5 and len(covered)==240
def coset_action(A):
    return tuple(next(i for i,c in enumerate(cosets) if key(A@r) in c) for r in reps)
acts={coset_action(A) for A in G}
assert len(acts)==60
assert all(sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0 for p in acts)
assert coset_action(1j*I)==tuple(range(5))
quarter=next(A for A in B if np.array_equal(A@A,-I) and sum(i==x for i,x in enumerate(coset_action(A)))==1)
qperm=coset_action(quarter)
assert all(qperm[qperm[i]]==i for i in range(5))
_,raw,_=n.reflection_actions(zs)
outer,_=n.outer_mark_action(raw)
def compose(p,q):return tuple(p[q[a]] for a in range(len(p)))
def pgroup(gens):
    ident=tuple(range(len(gens[0])));found={ident};todo=[ident]
    while todo:
        p=todo.pop()
        for q in gens:
            r=compose(q,p)
            if r not in found:found.add(r);todo.append(r)
    return found
def orbit_sizes(group):
    remaining=set(range(len(next(iter(group)))));out=[]
    while remaining:
        a=min(remaining);orb={p[a] for p in group};remaining-=orb;out.append(len(orb))
    return sorted(out)
pg=pgroup([compose(outer[i],outer[j]),compose(outer[k],outer[l])])
pr=pgroup([compose(raw[i],raw[j]),compose(raw[k],raw[l])])
assert len(pg)==len(pr)==60
assert orbit_sizes(pg)==[1,5] and orbit_sizes(pr)==[6]
assert all(p[0]==0 for p in pg)
quartic_trace=0
for A in G:
    t1=np.trace(A);t2=np.trace(A@A);t3=np.trace(A@A@A);t4=np.trace(A@A@A@A)
    quartic_trace+=(t1**4+6*t1*t1*t2+3*t2*t2+8*t1*t3+6*t4)/24
assert quartic_trace==480
assert sum(sum(p[a]==a for a in range(6))-1 for p in pg)==60
# Equivariant source-event dictionary, including the actual mark endpoint.
raylookup={n.canonical_mu4(z):a for a,z in enumerate(zs)}
raygens=[tuple(raylookup[n.canonical_mu4(A@z)] for z in zs) for A in (g,h)]
markgens=[compose(outer[i],outer[j]),compose(outer[k],outer[l])]
fivegens=[tuple(p[a]-1 for a in range(1,6)) for p in markgens]
raygroup=pgroup(raygens)
remaining=set(range(60));rayorbits=[]
while remaining:
    a=min(remaining);orb={p[a] for p in raygroup};remaining-=orb;rayorbits.append(orb)
assert sorted(map(len,rayorbits))==[10,20,30]
selected={a for a,p in enumerate(outer) if p[0]!=0}
assert selected==next(o for o in rayorbits if len(o)==20)
targets={20:list(it.permutations(range(5),2)),10:list(it.combinations(range(5),2)),30:[(a,b,c) for a,b in it.combinations(range(5),2) for c in range(5) if c not in (a,b)]}
def target_action(t,p,size):
    if size==20:return tuple(p[a] for a in t)
    if size==10:return tuple(sorted(p[a] for a in t))
    a,b,c=t;return tuple(sorted((p[a],p[b])))+(p[c],)
event_dictionary={}
for orb in rayorbits:
    size=len(orb);seed=min(orb);mapping=None
    for target in targets[size]:
        if size==20 and target[1]!=outer[seed][0]-1:continue
        trial={seed:target};todo=[seed];ok=True
        while todo and ok:
            r=todo.pop()
            for rg,fg in zip(raygens,fivegens):
                r2=rg[r];t2=target_action(trial[r],fg,size)
                if r2 in trial and trial[r2]!=t2:ok=False;break
                if r2 not in trial:trial[r2]=t2;todo.append(r2)
        if ok and len(trial)==size and len(set(trial.values()))==size:
            mapping=trial;break
    assert mapping is not None
    if size==20:assert all(t[1]==outer[r][0]-1 for r,t in mapping.items())
    event_dictionary[str(size)]={str(r):t for r,t in sorted(mapping.items())}
(HERE/'event_dictionary.json').write_text(json.dumps(event_dictionary,indent=2)+'\n')
def serialize(A):return [[[int((2*z).real),int((2*z).imag)] for z in row] for row in A]
result={'source_reflection_indices':[i,j,k,l],'generators_twice_gaussian_coordinates':{'g':serialize(g),'h':serialize(h),'J':serialize(1j*I)},'binary_subgroup_order':len(B),'extended_group_order':len(G),'center_order':len(center),'root_orbits':sorted(map(len,orbits),reverse=True),'relations':{'g2':'-I','h3':'-I','gh5':'-I'},'tetrahedral_normalizer_order':len(H),'coset_count':len(cosets),'coset_image_order':len(acts),'coset_image':'A5','noncentral_order4_coset_action':qperm,'coset_eigenspace_split':[3,2],'central_J_coset_action':'identity','scope':'Explicit native Gaussian E8 lattice subgroup; equality with a particular period-lattice Picard representation requires a marked intertwiner. Coset5 and its noncentral order4 lift are not identified with Picard fixed-point5 or its central quarter.'}
result['native_outer_chart_orbits']=[1,5]
result['native_outer_fixed_mark']=0
result['native_raw_odd_forms_orbits']=[6]
result['full_Sym4_fixed_dimension']=2
result['native_V5_fixed_dimension']=1
result['ray_orbits']=[10,20,30]
result['ray20_is_exactly_outer_mark_moving_events']=True
result['ray20_ordered_pair_second_coordinate_is_mark_endpoint']=True
result['equivariant_event_dictionary']='event_dictionary.json'
(HERE/'integral_lift.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
