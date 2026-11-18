"""Exact neutral bilocal reconstruction in the retained quarter source.

Finite-energy matrix elements of an infinite-mode operator identity.
These matrices are not an executed projection or a new dynamics.
"""
from pathlib import Path
import sys,json,argparse,hashlib
import sympy as s
REPO=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
ORIGINAL=REPO/'experiments/theory-contracts/joint-native-algebra-selection-20261016'
sys.path[:0]=[str(ORIGINAL),str(REPO)]
import check_original_current_schur_kernel as original

def need(ok,msg):
 if not bool(ok):raise RuntimeError(msg)

def partitions(n,top=None):
 if n==0:
  yield ();return
 for a in range(min(n,n if top is None else top),0,-1):
  for tail in partitions(n-a,a):yield(a,)+tail

def occupancy(lam,depth):
 return frozenset((lam[i-1] if i<=len(lam) else 0)-i for i in range(1,depth+1))

def truncated_exp(A,degree):
 out=s.eye(A.rows);term=s.eye(A.rows)
 for k in range(1,degree+1):
  term=term*A/k;out+=term
 need(term*A==s.zeros(A.rows),'Strict graded nilpotence')
 return out

def one_case(degree,z,w):
 basis=[lam for n in range(degree+1) for lam in partitions(n)]
 index={lam:i for i,lam in enumerate(basis)}
 depth=degree+3;occ=[occupancy(lam,depth) for lam in basis]
 occindex={v:i for i,v in enumerate(occ)};vac=occ[0]
 raising=s.zeros(len(basis));lowering=s.zeros(len(basis))
 for mode in range(1,degree+1):
  R=s.zeros(len(basis))
  for col,lam in enumerate(basis):
   for mu,sign in original.add_ribbons(lam,mode):
    if mu in index:R[index[mu],col]+=sign
  raising+=(z**mode-w**mode)*R/mode
  lowering+=(w**(-mode)-z**(-mode))*R.T/mode
 C=w/(z-w)
 boson=C*truncated_exp(raising,degree)*truncated_exp(lowering,degree)
 # Independent exact CAR occupancy computation with the infinite diagonal
 # sea summed analytically; no finite-sea replacement.
 car=s.zeros(len(basis))
 for col,O in enumerate(occ):
  car[col,col]=C+sum((z/w)**r for r in O-vac)-sum((z/w)**r for r in vac-O)
  for remove in O:
   for add in range(-depth,degree+1):
    if add in O:continue
    target=O-{remove}|{add}
    if target not in occindex:continue
    sign=(-1)**sum(min(remove,add)<r<max(remove,add) for r in O)
    car[occindex[target],col]+=sign*z**add*w**(-remove)
 need(car==boson,'Full finite-energy matrix equality of original bilocal and current vertex')
 need(car!=C*s.eye(len(basis)),'Reconstruction is not a scalar vacuum contraction')
 twist=[]
 for alpha in (s.Rational(1,4),s.Rational(3,4)):
  factor=(z/w)**alpha
  need((factor*car)==factor*boson,'Same actual quarter factor')
  need(factor*car!=car,'Omitting the original quarter phase is detected')
  twist.append({'alpha':str(alpha),'factor':str(factor)})
 return {'degree':degree,'basis_size':len(basis),'matrix_entries_checked':len(basis)**2,'z':str(z),'w':str(w),'quarters':twist,'equal':True}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=Path(__file__).with_name('bilocal_current_results.json'));args=parser.parse_args()
 cases=[one_case(d,z,w) for d in (4,8) for z,w in [(s.Rational(1,2),s.Rational(1,32)),(s.Rational(3,4),s.Rational(3,64))]]
 cauchy=[]
 for n in (1,2,4,8):
  zs=[s.Rational(1,2)+s.Rational(k,64) for k in range(n)];ws=[s.Rational(k,128) for k in range(1,n+1)]
  K=s.Matrix([[w/(z-w) for w in ws] for z in zs])
  numerator=s.prod(ws)*s.prod((zs[i]-zs[j])*(ws[j]-ws[i]) for i in range(n) for j in range(i+1,n))
  denominator=s.prod(z-w for z in zs for w in ws)
  need(K.det()==numerator/denominator,'Full neutral multi-fermion Cauchy/vertex contraction')
  cauchy.append({'pairs':n,'exact_match':True,'nonzero':K.det()!=0})
 result={'verdict':'EXACT_SCOPED_CURRENT_VERTEX_RECONSTRUCTION; ORIGINAL_P1_CURRENT_IMAGE_OPEN','cases':cases,'neutral_multi_fermion_checks':cauchy,'source_sha256':{str((ORIGINAL/'check_original_current_schur_kernel.py').relative_to(REPO)):hashlib.sha256((ORIGINAL/'check_original_current_schur_kernel.py').read_bytes()).hexdigest()},'scope':'Original infinite sea kept analytically. Finite-energy matrix elements checked, not a finite physical Hilbert space. Current strings require connected containing intervals or a specified relative-commutant extension, not the additive eight-arc subnet. No primitive field-dependent coupling or P1 map supplied.'}
 args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
