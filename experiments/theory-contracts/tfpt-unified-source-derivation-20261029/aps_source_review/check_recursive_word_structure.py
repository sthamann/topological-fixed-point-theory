"""Exact structure of the already derived original word realization."""
from pathlib import Path
import json
import sympy as s
p=Path(__file__).parent
d=json.loads((p/'RECURSIVE-WORD-SOURCE.json').read_text())
A=s.Matrix(d['operators']['L']);J=s.Matrix(d['operators']['J'])
G=s.diag(*map(s.sympify,d['metric_diagonal']));N=G[:5,:5]
C=J[5:,5:];U=A[:5,5:];B=A[5:,:5]
powers=[s.eye(63),C,C*C,C**3]
ws=[s.Integer(1),s.I,s.Integer(-1),-s.I];forms=[];dims=[]
for w in ws:
    E=sum((w**(-r)*powers[r] for r in range(4)),s.zeros(63))/4
    forms.append(U*E*B);dims.append(s.trace(E))
def mat(M):return [[str(a) for a in row] for row in M.tolist()]
defect=s.Matrix(d['clock_ordered_word_defect'])
minors=[(-N*defect)[:k,:k].det() for k in range(1,6)]
checks={
    'memory_dimension63':sum(dims)==63,
    'sector_sum_is_original_memory_defect':sum(forms,s.zeros(5))==(A*A)[:5,:5]-A[:5,:5]**2,
    'weighted_sector_sum_is_clock_defect':sum(((w-1)*F for w,F in zip(ws,forms)),s.zeros(5))==defect,
    'clock_defect_metric_Hermitian':(N*defect).H==N*defect,
    'negative_clock_defect_strictly_positive':all(x>0 for x in minors),
    'opposite_imaginary_sector_input_forms_equal':forms[1]==forms[3],
}
out={'contract':d['contract'],'all_passed':all(checks.values()),'checks':checks,
     'memory_clock_multiplicities':dict(zip(map(str,ws),map(int,dims))),
     'clock_sector_input_Gram_forms':{str(w):mat(F) for w,F in zip(ws,forms)},
     'input_form_ranks':[int(F.rank()) for F in forms],
     'code_metric':list(map(str,N.diagonal())),
     'negative_clock_defect_leading_minors':list(map(str,minors)),
     'original_unit_f4_response_amplitudes':{
         'L_L':str((A*A)[4,4]),'L_J_L':str((A*J*A)[4,4]),
         'incorrect_code_only_L_L':str((A[:5,:5]**2)[4,4])},
     'scope':'Same original single-beta L and native J alphabet only. Sector dimensions are not physical spatial dimensions or particle counts.'}
(p/'RECURSIVE-WORD-SOURCE-STRUCTURE.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='clock_sector_input_Gram_forms'},indent=2))
if not out['all_passed']:raise RuntimeError('Exact original word structure check failed')
