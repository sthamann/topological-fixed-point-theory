"""Reuse the actual full96 Higgs function; spatial Fourier cut, never OS time."""
from pathlib import Path
import ast, hashlib, importlib.util, json
from fractions import Fraction
import numpy as np
from repo_paths import repo_root,original_path
ROOT=repo_root()
V252=ROOT/'verification/v252_full_finite_triple.py'
HIGGS=ROOT/'experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/check_original_higgs_fiber.py'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_original():
    import sys
    sys.path.insert(0,str(ROOT/'verification'))
    spec=importlib.util.spec_from_file_location('field_original_v252',V252)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    # Extract the unchanged owning function only: importing its checker would
    # execute all controls and overwrite its canonical result file.
    tree=ast.parse(HIGGS.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='higgs')
    z=np.zeros((3,3),complex);ns={'original':mod,'z3':z,'K':mod.Kswap(),'P':mod.rep(1,np.zeros((2,2)),z)}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(HIGGS),'exec'),ns)
    return mod,ns['higgs']
def decode_matrix(a):
    v=np.asarray(a,dtype=float)
    if v.shape!=(3,3,2)or not np.isfinite(v).all():raise ValueError('Expected finite3x3 complex matrix encoded as [real,imag] pairs')
    return v[:,:,0]+1j*v[:,:,1]
def load_domain(path):
    p=Path(path).resolve();d=json.loads(p.read_text());d['_path']=str(p);d['_sha256']=sha(p);d['domain_directory']=str(p.parent)
    if d.get('schema')!='tfpt.spatial-field-search-domain.v1':raise ValueError('Unknown domain schema')
    if d.get('coordinate_role')!='spatial_tangential_circle_NOT_OS_time':raise ValueError('Spatial/time type mismatch')
    for key in ['source_originals','couplings','quarter','modes','field_proposal']:
        if key not in d:raise ValueError('Missing original domain input '+key)
    d['_resolved_source_originals']=[]
    for item in d['source_originals']:
        resolved=original_path(item['path'],ROOT)
        if sha(resolved)!=item['sha256']:raise ValueError('Original source changed: '+item['path'])
        d['_resolved_source_originals'].append(str(resolved))
    modes=d['modes']
    if not modes or any(type(n)is not int for n in modes)or sorted(set(modes))!=modes:raise ValueError('Actual distinct ordered integer modes required')
    if len(modes)%2 or modes!=list(range(-len(modes)//2,len(modes)//2)):raise ValueError('Use the actual symmetric original quarter cut [-N,N)')
    prop=d['field_proposal'];multipliers=[Fraction(x)for x in prop['allowed_multipliers']]
    if not multipliers or len(set(multipliers))!=len(multipliers):raise ValueError('Distinct finite rational proposal multipliers required')
    if len(prop['base_constant'])!=4:raise ValueError('All four real H components required')
    for x in prop['base_constant']:Fraction(x)
    if not prop['slots']:raise ValueError('No field search coordinates supplied')
    seen=set()
    for slot in prop['slots']:
        if set(slot)!={'kind','harmonic','component','unit'}:raise ValueError('Unexpected field coordinate key')
        key=(slot['kind'],slot['harmonic'],slot['component'])
        if key in seen:raise ValueError('Duplicate Fourier coordinate would duplicate physical genotypes')
        seen.add(key)
        if type(slot['harmonic'])is not int or not 0<=slot['harmonic']<=modes[-1]-modes[0]:raise ValueError('Fourier harmonic invisible at this original cutoff')
        if type(slot['component'])is not int or not 0<=slot['component']<4:raise ValueError('H component out of range')
        if slot['kind']not in ('constant','cos','sin')or(slot['harmonic']==0)!=(slot['kind']=='constant'):raise ValueError('Fourier coordinate type mismatch')
        if Fraction(slot['unit'])==0:raise ValueError('Zero Fourier unit would duplicate genotypes')
    if d['quarter']not in ['1/4','3/4']:raise ValueError('Only the two declared original quarter domains supported')
    Y=[decode_matrix(d['couplings'][k])for k in ['Yu','Yd','Ye','Ynu']];MR=decode_matrix(d['couplings']['MR'])
    if not np.allclose(MR,MR.T,rtol=0,atol=1e-13):raise ValueError('Original Majorana symmetry failed')
    return d,Y,MR
class Assembler:
    def __init__(self,domain_path):
        self.domain,self.Y,self.MR=load_domain(domain_path);self.original,self.higgs=load_original()
        self.gamma=self.original.grading();self.K=self.original.Kswap()
        self.DM=self.original.full_D(*[np.zeros((3,3),complex)for _ in range(4)],self.MR)
        Q=[np.array([[0,1],[-1,0]],complex),np.array([[0,1j],[1j,0]],complex),np.eye(2),np.diag([-1j,1j])]
        self.V=[self.higgs(q,self.Y,self.MR)-self.DM for q in Q]
    def coefficients(self,candidate):
        prop=self.domain['field_proposal'];values=[Fraction(v)for v in prop['allowed_multipliers']]
        base=[Fraction(x)for x in prop['base_constant']];coeff={0:np.array([float(x)for x in base],complex)}
        for gene,slot in zip(candidate['genes'],prop['slots']):
            amp=values[gene]*Fraction(slot['unit']);a=slot['component'];harm=slot['harmonic'];kind=slot['kind']
            if harm==0:
                coeff[0][a]+=float(amp)
            else:
                coeff.setdefault(harm,np.zeros(4,complex));coeff.setdefault(-harm,np.zeros(4,complex))
                plus=float(amp)/2*(1 if kind=='cos' else -1j)
                coeff[harm][a]+=plus;coeff[-harm][a]+=plus.conjugate()
        return coeff
    def assemble(self,candidate):
        coeff=self.coefficients(candidate);phi={k:sum((v*x for v,x in zip(self.V,xs)),np.zeros((96,96),complex))for k,xs in coeff.items()}
        phi[0]=phi.get(0,np.zeros((96,96),complex))+self.DM
        ns=self.domain['modes'];alpha=float(Fraction(self.domain['quarter']));z=np.zeros((96,96),complex)
        # Same original Fourier multiplication rule mult(coeff,ns), now all96.
        pot={k:-1j*self.gamma@v for k,v in phi.items()}
        h=np.block([[pot.get(n-m,z)+(n+alpha)*self.gamma if n==m else pot.get(n-m,z)for m in ns]for n in ns])
        zero=lambda a:float(np.max(np.abs(a)))
        residuals={'h_N_selfadjoint':zero(h-h.conj().T),'Majorana_symmetric':zero(self.MR-self.MR.T),
          'Phi_fourier_adjoint':max(zero(v-phi.get(-k,z).conj().T)for k,v in phi.items()),
          'Phi_fourier_grading':max(zero(self.gamma@v+v@self.gamma)for v in phi.values()),
          'Phi_fourier_real_structure':max(zero(self.K@v.conj()@self.K-phi.get(-k,z))for k,v in phi.items())}
        coupled=h.copy()
        for i in range(len(ns)):coupled[96*i:96*(i+1),96*i:96*(i+1)]=0
        return {'candidate_id':hashlib.sha256(json.dumps(candidate,sort_keys=True).encode()).hexdigest(),
          'domain_id':self.domain['domain_id'],'native_inputs_sha256':self.domain['_sha256'],
          'h_N_role':'original_tangential_one_particle_NOT_joint_OS_generator',
          'h_N':h,'phi_fourier':phi,'higgs_fourier':coeff,'gamma96':self.gamma,'K96':self.K,
          'modes':ns,'alpha':alpha,'kinematic_residuals':residuals,
          'diagnostics':{'full_internal_dimension':96,'matrix_dimension':len(h),'mode_mixing_frobenius_norm':float(np.linalg.norm(coupled)),
          'spatial_coordinate_is_OS_time':False,'numerical_rounding_enclosed':False}}
