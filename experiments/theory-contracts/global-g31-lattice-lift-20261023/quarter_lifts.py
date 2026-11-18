"""All strict character-valued lifts of the original 62 clocked objects."""
from collections import Counter, deque
from pathlib import Path
import json
import pickle
import runpy
import numpy as np
import sympy as S

HERE=Path(__file__).resolve().parent

def main():
    p=runpy.run_path(str(HERE/'probe.py'))
    refs,perms,pg,chosen,orbits,full,orbit_id,rays=pickle.loads((p['OLD']/'state.pkl').read_bytes())
    old=json.loads((p['OLD']/'results.json').read_text())
    C=np.array(json.loads((HERE/'probe-results.json').read_text())['C'],dtype=np.int64)
    B=np.array(old['integral_lattice_basis'],dtype=np.int64)
    Bi=np.array((2*S.Matrix(B.tolist()).inv()).tolist(),dtype=np.int64)
    I=(2*np.eye(4,dtype=np.int64),np.zeros((4,4),dtype=np.int64));ki=p['key'](I)
    cache={}
    def data(g):
        k=p['key'](g)
        if k not in cache:
            v=Bi@p['real_twice'](g)@B
            if np.any(v%4):raise RuntimeError('nonintegral original action')
            M=(v//4)%2;cols=[sum(int(M[i,j])<<i for i in range(8)) for j in range(8)]
            h=p['phase_table'](M,C)
            action=[p['act'](cols,x) for x in range(256)];inv=[0]*256
            for x,y in enumerate(action):inv[y]=x
            invcols=[inv[1<<i] for i in range(8)]
            chars=[sum(((a&v).bit_count()%2)<<i for i,v in enumerate(invcols)) for a in range(256)]
            cache[k]=(cols,h,chars)
        return cache[k]
    def commute(a,b):return p['key'](p['mul'](a,b))==p['key'](p['mul'](b,a))
    def row_basis(rows):
        bs={}
        for x in rows:
            while x:
                k=x.bit_length()-1
                if k in bs:x^=bs[k]
                else:bs[k]=x;break
        return bs
    def reduce(v,bs):
        for k in sorted(bs,reverse=True):
            if (v>>k)&1:v^=bs[k]
        return v
    basecent={}
    for oi in (0,12,15):
        rr=[refs[j] for j in set(old['orbits'][oi]['representative'])]
        basecent[oi]=[g for g in full.values() if all(commute(g,r) for r in rr)]
    result=[]
    for oi,o in enumerate(old['orbits']):
        for ci,c in enumerate(o.get('clock_classes',[])):
            if c['order']!=4:continue
            a=(np.array(c['matrix_real_twice'],dtype=np.int64),np.array(c['matrix_imag_twice'],dtype=np.int64))
            b=p['mul'](refs[o['representative'][0]],a);gs=[a,b]
            gd=[data(g) for g in gs];rows=[]
            for si,(cols,phase,_) in enumerate(gd):
                for i in range(8):
                    x=1<<i;norm=value=0
                    for _ in range(4):norm^=x;value^=phase[x];x=p['act'](cols,x)
                    if x!=1<<i:raise RuntimeError('wrong fourth power')
                    rows.append((norm<<(8*si),value))
            solutions=[v for v in range(65536) if all(((v&r).bit_count()%2)==t for r,t in rows)]
            need=len(solutions) in (16384,65536)
            if not need:raise RuntimeError('unexpected based lift count')
            cob=[]
            for j in range(8):
                mask=0
                for si,(cols,_,_) in enumerate(gd):
                    mask|=sum(((((1<<i)^cols[i])>>j)&1)<<(8*si+i) for i in range(8))
                cob.append(mask)
            bs=row_basis(cob);quot={reduce(v,bs) for v in solutions}
            if len(quot)*(1<<len(bs))!=len(solutions):raise RuntimeError('H1 quotient size')
            cent=[g for g in basecent[oi] if commute(g,a)]
            if len(cent)!=c['automorphism_order']:raise RuntimeError('wrong actual centralizer')
            # Actual matrix generators of this centralizer.
            generated={ki:I};generators=[]
            for g in cent:
                if p['key'](g) in generated:continue
                generators.append(g);todo=list(generated.values())
                while todo:
                    x=todo.pop()
                    for h in generators:
                        z=p['mul'](h,x);k=p['key'](z)
                        if k not in generated:generated[k]=z;todo.append(z)
            if len(generated)!=len(cent):raise RuntimeError('centralizer generator closure')
            actions=[]
            for h in generators:
                hc,hp,hchars=data(h);shifts=[]
                for gc,gp,_ in gd:
                    shifts.append(sum((hp[gc[i]]^gp[hc[i]])<<i for i in range(8)))
                actions.append((hchars,shifts))
            def transform(v,action):
                lookup,shift=action
                return reduce(lookup[(v&255)^shift[0]]|(lookup[(v>>8)^shift[1]]<<8),bs)
            remaining=set(quot);classes=[]
            while remaining:
                root=min(remaining);remaining.remove(root);todo=[root];orb={root}
                while todo:
                    x=todo.pop()
                    for ac in actions:
                        y=transform(x,ac)
                        if y not in quot:raise RuntimeError('lift space not invariant')
                        if y not in orb:
                            orb.add(y);todo.append(y);remaining.remove(y)
                classes.append({'phase_pair':root,'orbit_size_mod_kernel':len(orb),
                                'full_automorphism_order':256*len(cent)//((1<<len(bs))*len(orb))})
            if sum(q['orbit_size_mod_kernel'] for q in classes)!=len(quot):raise RuntimeError('orbit exhaustion')
            item={'monodromy_id':oi,'clock_id':ci,'based_lift_count':len(solutions),
                  'kernel_coboundary_rank':len(bs),'kernel_conjugacy_classes':len(quot),
                  'base_automorphism_order':len(cent),'full_lift_isomorphism_classes':len(classes),
                  'normal_forms':classes,'power_equations':rows,'kernel_coboundary_basis':list(bs.values())}
            result.append(item)
            print(oi,ci,'based',len(solutions),'H1',len(quot),'unbased',len(classes),flush=True)
    out={'contract':'MATH.SEAM.G31.LATTICE.03','original_clock_classes':62,
         'all_strict_character_lifts_exist':len(result)==62,
         'based_count_distribution':dict(Counter(v['based_lift_count'] for v in result)),
         'total_unbased_lift_classes':sum(v['full_lift_isomorphism_classes'] for v in result),
         'classes':result}
    (HERE/'quarter-lifts.json').write_text(json.dumps(out,indent=2)+'\n')
    print('TOTAL',out['total_unbased_lift_classes'],out['based_count_distribution'])

if __name__=='__main__':main()
