"""Reproduce an EXISTING original diagnostic fixture, never physical inputs."""
from pathlib import Path
import hashlib,json
import numpy as np
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'runtime'))
from repo_paths import repo_root
R=repo_root()
paths=[R/'verification/v252_full_finite_triple.py',R/'experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/check_original_higgs_fiber.py',R/'experiments/theory-contracts/original-local-source-variable-background-20261102/check_variable_background.py']
rng=np.random.default_rng(252)
Ys=[rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))for _ in range(4)]
M=rng.normal(size=(3,3));MR=M+M.T
def enc(a):return np.stack([a.real,a.imag],axis=-1).tolist()
d={'schema':'tfpt.spatial-field-search-domain.v1','domain_id':'original-variable-H-quarter-full96-diagnostic-v1',
 'input_role':'existing_original_diagnostic_fixture_NOT_physical_Y_selection','physical_source_identity':False,
 'coordinate_role':'spatial_tangential_circle_NOT_OS_time','quarter':'1/4','modes':[-2,-1,0,1],
 'couplings':dict(zip(['Yu','Yd','Ye','Ynu','MR'],[enc(a)for a in Ys+[MR]])),
 'source_originals':[{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in paths],
 'fixture_origin':{'Y_MR':'v252.run rng252 fixture unchanged; not physical supplied Yukawas',
 'variable_H':'check_variable_background.py potential[1]=.065*T0-.045j*T1, potential[2]=.035*T3, baseT2',
 'extension':'same completePhi(H) family and originalquarter; all96 retained instead of reducingquark36'},
 'field_proposal':{'role':'finite proposal study, NOT original integration measure or physical amplitude bound',
 'base_constant':['0','0','1','0'],'allowed_multipliers':['-1','0','1'],
 'slots':[{'kind':'cos','harmonic':1,'component':0,'unit':'13/100'},
 {'kind':'sin','harmonic':1,'component':1,'unit':'9/100'},
 {'kind':'cos','harmonic':2,'component':3,'unit':'7/100'}]},
 'matrix_residual_tolerance':1e-10,'whole_source_bundle':None,
 'ensemble_requirements':['original_joint_weighted_Z0_ZU_ZY_ZUY','original_OS_history_RP','actual_quantum_Ward_SD_contacts','whole_physical_P_Ref_Eta_domains','original_regulator_limit_order']}
p=Path(__file__).with_name('native-diagnostic-v1.json');p.write_text(json.dumps(d,indent=2)+'\n');print(p)
