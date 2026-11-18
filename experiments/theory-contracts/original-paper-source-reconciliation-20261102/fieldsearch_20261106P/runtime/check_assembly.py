"""Small original assembly controls, not a new general test harness."""
from pathlib import Path
import json
import numpy as np
from assemble import Assembler
BASE=Path(__file__).resolve().parents[1]
A=Assembler(BASE/'domain/native-diagnostic-v1.json')
def candidate(genes):return {'schema':'tfpt.field-configuration.v1','domain_sha256':A.domain['_sha256'],'genes':genes}
a=A.assemble(candidate([1,1,1]));b=A.assemble(candidate([2,2,2]))
radial=A.original.full_D(*A.Y,A.MR)
checks={
 'baseline_equals_original_full_D':float(np.max(np.abs(a['phi_fourier'][0]-radial)))<1e-12,
 'all96_in_each_native_mode':b['h_N'].shape==(384,384),
 'actual_nonconstant_mode_mixing_retained':b['diagnostics']['mode_mixing_frobenius_norm']>0,
 'original_MR_retained_at_zero_mode':float(np.max(np.abs(a['phi_fourier'][0]-A.DM-sum(A.V[i]*x for i,x in enumerate([0,0,1,0])))))<1e-12,
 'all_matrix_constraints_pass':max(b['kinematic_residuals'].values())<1e-10,
 'spatial_not_OS_time':b['diagnostics']['spatial_coordinate_is_OS_time']is False}
assert all(checks.values()),checks
out={'status':'PASS_ORIGINAL_FULL96_SPATIAL_ASSEMBLY_ONLY','checks':checks,'nonconstant':b['diagnostics'],
 'residuals':b['kinematic_residuals'],'full_original_source_evaluated':False}
print(json.dumps(out,indent=2))
