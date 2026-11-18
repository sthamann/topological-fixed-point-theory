"""Independent Burnside counts using affine fixed-point equations."""
from collections import Counter
from pathlib import Path
import json
import pickle
import runpy
import numpy as np
import sympy as S

HERE=Path(__file__).resolve().parent

def linear_count(rows,n):
    basis={}
    for co,rhs in rows:
        while co:
            i=co.bit_length()-1
            if i in basis:
                a,b=basis[i];co^=a;rhs^=b
            else:basis[i]=(co,rhs);break
        if not co and rhs:return 0
    return 1<<(n-len(basis))

def main():
    p=runpy.run_path(str(HERE/'probe.py'))
    refs,perms,pg,chosen,orbits,full,orbit_id,rays=pickle.loads((p['OLD']/'state.pkl').read_bytes())
    old=json.loads((p['OLD']/'results.json').read_text())
    quarter=json.loads((HERE/'quarter-lifts.json').read_text())
    C=np.array(json.loads((HERE/'probe-results.json').read_text())['C'],dtype=np.int64)
    B=np.array(old['integral_lattice_basis'],dtype=np.int64)
    Bi=np.array((2*S.Matrix(B.tolist()).inv()).tolist(),dtype=np.int64)
    cache={}
    def data(g):
        key=p['key'](g)
        if key not in cache:
            M=(Bi@p['real_twice'](g)@B//4)%2
            cols=[sum(int(M[i,j])<<i for i in range(8)) for j in range(8)]
            cache[key]=(cols,p['phase_table'](M,C))
        return cache[key]
    def commutes(a,b):return p['key'](p['mul'](a,b))==p['key'](p['mul'](b,a))
    # Every canonical tuple starts with the original r0.
    cent0=[h for h in full.values() if commutes(h,refs[0])]
    original=[]
    for oi,o in enumerate(old['orbits']):
        gs=[refs[j] for j in o['representative'][:3]];gd=[data(g) for g in gs]
        cent=[h for h in cent0 if all(commutes(h,g) for g in gs)]
        total=0
        for h in cent:
            hc,hp=data(h);rows=[]
            for si,(gc,gp) in enumerate(gd):
                for i in range(8):
                    coeff=(((1<<i)^hc[i])<<(8*si)) ^ (((1<<i)^gc[i])<<24)
                    rows.append((coeff,hp[gc[i]]^gp[hc[i]]))
            total+=linear_count(rows,32)
        divisor=256*len(cent)
        if total%divisor:raise RuntimeError('nonintegral monodromy Burnside count')
        original.append({'monodromy_id':oi,'representative':o['representative'],
                         'based_lifts':1<<24,'base_centralizer_order':len(cent),
                         'burnside_fixed_pair_sum':total,'isomorphism_classes':total//divisor})
        print('F3',oi,'classes',total//divisor,flush=True)
    audited=[];direct_checks=0
    for case in quarter['classes']:
        oi,ci=case['monodromy_id'],case['clock_id'];o=old['orbits'][oi];c=o['clock_classes'][ci]
        a=(np.array(c['matrix_real_twice'],dtype=np.int64),np.array(c['matrix_imag_twice'],dtype=np.int64))
        b=p['mul'](refs[o['representative'][0]],a);gs=[a,b];gd=[data(g) for g in gs]
        cent=[h for h in cent0 if all(commutes(h,g) for g in gs)]
        total=0
        for h in cent:
            hc,hp=data(h);rows=[tuple(x) for x in case['power_equations']]
            for si,(gc,gp) in enumerate(gd):
                for i in range(8):
                    coeff=(((1<<i)^hc[i])<<(8*si)) ^ (((1<<i)^gc[i])<<16)
                    rows.append((coeff,hp[gc[i]]^gp[hc[i]]))
            total+=linear_count(rows,24)
        divisor=256*len(cent)
        if total%divisor or total//divisor!=case['full_lift_isomorphism_classes']:
            raise RuntimeError('independent clock Burnside count disagrees')
        # Check the affine conjugation formula directly on all eight generators,
        # using inversion/composition of full 256-point phase functions.
        for h in cent[:2]:
            hc,hp=data(h);ma=[p['act'](hc,x) for x in range(256)];inv=[0]*256
            for x,y in enumerate(ma):inv[y]=x
            for si,(gc,gp) in enumerate(gd):
                v=case['normal_forms'][0]['phase_pair']>>(8*si)&255
                shift=sum((hp[gc[i]]^gp[hc[i]])<<i for i in range(8))
                formula=sum((((v^shift)&inv[1<<i]).bit_count()%2)<<i for i in range(8))
                direct=0
                for i in range(8):
                    y=inv[1<<i];z=p['act'](gc,y)
                    val=hp[y]^gp[y]^((v&y).bit_count()%2)^hp[z]
                    direct|=val<<i
                if formula!=direct:raise RuntimeError('direct affine conjugation disagrees')
                direct_checks+=1
        audited.append({'monodromy_id':oi,'clock_id':ci,'burnside_fixed_pair_sum':total,
                        'divisor':divisor,'isomorphism_classes':total//divisor})
    out={'contract':'MATH.SEAM.G31.LATTICE.03','all62_orbit_counts_independently_confirmed':True,
         'direct_conjugation_checks':direct_checks,'quarter_total':sum(c['isomorphism_classes'] for c in audited),
         'quarter_classes':audited,'ordinary_total':sum(c['isomorphism_classes'] for c in original),
         'ordinary_classes':original}
    (HERE/'burnside-audit.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'ORDINARY_LIFTS.tsv').write_text('monodromy_id\treflection_indices\tbase_centralizer\tlift_isomorphism_classes\n'+
        ''.join(f"{c['monodromy_id']}\t{c['representative']}\t{c['base_centralizer_order']}\t{c['isomorphism_classes']}\n" for c in original))
    (HERE/'CLOCK_LIFTS.tsv').write_text('monodromy_id\tclock_id\tbased_lifts\tkernel_conjugacy_classes\tfull_isomorphism_classes\n'+
        ''.join(f"{c['monodromy_id']}\t{c['clock_id']}\t{c['based_lift_count']}\t{c['kernel_conjugacy_classes']}\t{c['full_lift_isomorphism_classes']}\n" for c in quarter['classes']))
    print('TOTALS',out['quarter_total'],out['ordinary_total'])

if __name__=='__main__':main()
