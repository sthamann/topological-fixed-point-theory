"""Exact first type test on existing E8 C/X fields; no chosen Higgs map."""
from pathlib import Path
from fractions import Fraction
import json
from tfpt_explorer.current_block_geometry import _actual_current_source
from tfpt_explorer.source_program import _c_currents

s=_actual_current_source()
rows=[]
for a in range(4):
    ci,cp=next(iter(_c_currents()[a].items()))
    c=next(root for root,index in s['chevalley'].ridx.items() if index==ci)
    for i in range(5):
        x=tuple(s['w20'][i,a])
        q=tuple(cj+xj for cj,xj in zip(c,x))
        dot=lambda p,r:sum(Fraction(u*v,4) for u,v in zip(p,r))
        rows.append({'i':i,'a':a,'C_doubled_charge':list(c),'X_doubled_charge':list(x),
                     'original_C_phase':str(cp),'original_X_phase':str(s['phases'][i,a]),
                     'composite_doubled_charge':list(q),'C_dot_X':str(dot(c,x)),
                     'norm_squared':str(dot(q,q)),'conformal_weight':str(dot(q,q)/2),
                     'D5_vector_charge':q[:5]==tuple(-2 if j==i else 0 for j in range(5)),
                     'family_doubled_charge':list(q[5:]),
                     'physical_Higgs_identified':False})
assert all(r['C_dot_X']=='0' and r['norm_squared']=='4' and r['D5_vector_charge'] for r in rows)
out={'scope':'original chiral E8 source composite :C_a X_ia:, not a new physical source',
     'exact_pairs_checked':len(rows),'all_pass':True,'selected_family_line':None,
     'same_original_contact_time_Higgs_identification':None,'rows':rows}
p=Path(__file__).with_name('native-composite-first-type-test.json')
p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
