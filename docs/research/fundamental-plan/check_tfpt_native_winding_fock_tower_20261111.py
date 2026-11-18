"""Full longitudinal tower for the given v252/R032 profile and its free reference."""
from pathlib import Path
from contextlib import redirect_stdout
import io
import sys
import numpy as np
from scipy.linalg import eigh, svdvals
from scipy.special import polygamma

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
source_path=ROOT/'experiments/theory-contracts/chat-triangulation-20261108/check_native_quarter_winding_and_geometry_20261111.py'
saved=sys.argv
sys.argv=[str(source_path),str(ROOT)]
original={'__file__':str(source_path)}
with redirect_stdout(io.StringIO()):
    exec(compile(source_path.read_text().split('baseline, _ = relative_cross()',1)[0],
                 str(source_path),'exec'),original)
sys.argv=saved
G=original['grading']
Id=original['identity']
source=original['source']
dt=np.array(source['durations'])
mass=np.array([-1j*G@(source['inner_form'](profile)+source['majorana'])
               for profile in source['history']])
assert len(G)==96 and np.linalg.matrix_rank((Id-G)/2)==48
assert np.max(np.abs(mass-mass.conj().transpose(0,2,1)))<1e-13
assert np.max(np.abs(mass@G+G@mass))<1e-13
def largest_logs(steps):
    prod=Id.copy()
    for x in steps:
        prod=x@prod
    return np.log(svdvals(prod)[:48])

def data(p):
    forward,backward=[],[]
    sea=top=bottom=0.
    for m,t in zip(mass,dt):
        vals,U=eigh(p*G+m)
        forward.append((U*np.exp(-t*(vals-vals[0])))@U.conj().T)
        backward.append((U*np.exp(t*(vals-vals[-1])))@U.conj().T)
        sea+=t*vals[:48].sum()
        top+=t*vals[-1]
        bottom+=t*vals[0]
    lead_forward=largest_logs(forward)
    logs=np.r_[2*(lead_forward-bottom), -2*(largest_logs(reversed(backward))+top)]
    if np.min(logs[:48])+1e-9<np.max(logs[48:]):
        raise ArithmeticError('Overlapping halves of the singular spectrum')
    # Use exact filled-sea convention. Cancel large linear parts before summing.
    free=96*np.log1p(np.exp(-2*dt.sum()*abs(p)))
    filled=2*(sea-48*bottom)+2*lead_forward.sum()
    log_empty=filled+np.logaddexp(0, -logs[:48]).sum()+np.logaddexp(0,logs[48:]).sum()-free
    return logs,log_empty

def twist(logs, v):
    a=np.exp(-np.abs(logs))
    positive=(logs>=0)
    correction=np.where(positive,v.conjugate(),v)
    return (positive.sum()%4)*np.log(v)+np.log1p(correction*a).sum()-np.log1p(a).sum()

cutoff=256
values={}
for n in range(-cutoff,cutoff):
    p=n+.25
    logs,ell=data(abs(p))
    if p<0:
        logs=-logs
    values[n]=(ell,twist(logs,1j),logs)
# Mirror momentum on alpha=3/4 reverses singular-log signs and twist.
for n in [-1,0,1,4]:
    p=n+.25
    forward,ell=data(p)
    assert abs(ell-values[n][0])<1e-9
    assert np.max(np.abs(np.sort(forward)-np.sort(values[n][2])))<1e-9
old, block=original['relative_cross']()
probe=data(.25)[1]+data(.75)[1]
assert abs(probe-old['relative_log_empty'])<1e-9
ratio0=np.exp(twist(data(.25)[0],1j)+twist(data(.75)[0],-1j))
print('Original two-mode log norm:',probe)
print('Original two-mode full CAR empty winding cross / diagonal:',ratio0)
for N in [1,2,4,8,16,32,64,128,cutoff]:
    a=2*sum(values[n][0] for n in range(-N,N))
    b=2*sum(values[n][1] for n in range(-N,N))
    print('Symmetric tower cutoff N',N,'log diagonal',a,'log |cross/diagonal|',b.real,
          'phase principal',np.angle(np.exp(1j*b.imag)))
lead=sum(np.linalg.norm(b-a,'fro')**2 for a,b in zip(mass[:-1],mass[1:]))/8
tail=polygamma(1,cutoff+.25)+polygamma(1,cutoff+.75)
print('Leading 1/p^2 coefficient of log norm:',-lead)
print('Log diagonal with leading tail added:',2*sum(x[0] for x in values.values())-2*lead*tail)
print('Durations:',dt.tolist(),'minimum:',dt.min(),'max mass norm:',max(np.linalg.norm(x,2) for x in mass))

# Bounds on the omitted modes. Anticommutation gives P(p) in closed form.
# These bound the analytic tail; finite linear algebra remains numerical.
p0=cutoff+.25
near=0.
fourth=0.
for A,B in zip(mass[:-1],mass[1:]):
    a=np.linalg.norm(A-B,'fro')/2
    c=(np.linalg.norm(A,2)+np.linalg.norm(B,2))/2
    F=(np.linalg.norm(A,'fro')+np.linalg.norm(B,'fro'))/2
    D=(np.linalg.norm(A@A,'fro')+np.linalg.norm(B@B,'fro'))/4
    E=(np.linalg.norm(A@A@A,'fro')+np.linalg.norm(B@B@B,'fro'))/4
    near+=F**2/2
    fourth+=(2*a*E+D**2+E**2/p0**2)/2+F**2*c**2/(4*(1-c**2/p0**2))
dmin=dt.min()
C=6*np.sqrt(96*np.exp(96*np.exp(-2*dmin*p0)))
eps0=C*np.exp(-dmin*p0)
baseline=np.sqrt(1-near/p0**2)
low=baseline-eps0
assert low>.9 and 2*eps0**2/low**2<.5
fourth_sum=(polygamma(3,cutoff+.25)+polygamma(3,cutoff+.75))/6
def geometric(rate):
    return (np.exp(-rate*(cutoff+.25))+np.exp(-rate*(cutoff+.75)))/(-np.expm1(-rate))
diag_err=2*(fourth*fourth_sum+2*C/low*geometric(dmin)+96*geometric(2*dt.sum()))
cross_err=2*(4*C**2/low**2)*geometric(2*dmin)
log_diag=2*sum(x[0] for x in values.values())-2*lead*tail
log_cross_rel=2*sum(x[1] for x in values.values())
print('Analytic omitted-mode error bound, log diagonal after leading tail:',diag_err)
print('Analytic omitted-mode error bound, complex log cross/diagonal:',cross_err)
print('Computed fixed-tower relative diagonal:',np.exp(log_diag))
print('Computed fixed-tower cross/diagonal:',np.exp(log_cross_rel))
print('Computed fixed-tower relative absolute cross:',np.exp(log_diag+log_cross_rel))

# Original one-way weak pulse on the lowest pair of quarter modes.
# On a mixed Fock overlap, the same normal-ordered two-point derivative
# is d(I-S^-1) = S^-1 dS S^-1, S=I+K_left^* K_right.
from scipy.linalg import solve
a,b=[original['arm'](p,False)[0] for p in (.25,.75)]
marker=original['raising']
ix=np.ix_(original['u'],original['d'])
def marked(A,B,M,va,vb):
    Sa=Id+va*A.conj().T@A
    Sb=Id+vb*B.conj().T@B
    Y=va*B.conj().T@M@A
    X=solve(Sa.T,solve(Sb,Y).T).T
    zero=np.zeros_like(Id)
    S=np.block([[Sa,zero],[zero,Sb]])
    D=np.block([[zero,zero],[Y,zero]])
    dd=(np.linalg.inv(S)-np.linalg.inv(S+D))[96:,:96]
    assert np.max(np.abs(dd-X))<1e-12
    return X
A,B=[original['direct'](x) for x in (a,b)]
assert np.max(np.abs(marked(A,B,marker,1,1)[ix]-block))<1e-13
mix_marked=marked(A,B,marker,1j,-1j)
for cuts in [(2,),(4,),(2,4)]:
    Ao,Bo=[original['sew'](x,cuts) for x in (a,b)]
    frames=original['frames']
    Mo=frames[-1].conj().T@marker@frames[-1]
    transformed=frames[0]@marked(Ao,Bo,Mo,1j,-1j)@frames[0].conj().T
    assert np.max(np.abs(mix_marked-transformed))<1e-12
block_m=mix_marked[ix]
print('Mixed full Fock minor: weakest-to-up block singular values:',svdvals(block_m).tolist())
print('Max entry of the lowest marked mixed minor relative to full free tower:',
      np.exp(log_diag+log_cross_rel.real)*np.max(np.abs(block_m)))

# Same native radius probe already used by the source, p -> p exp(-r).
# At fixed dt and mass profiles, the first UV coefficient is -A exp(2r).
for p in [64.25,128.25,256.25]:
    step=.001
    value=(data(p*np.exp(-step))[1]-data(p*np.exp(step))[1])/(2*step)
    assert abs(p*p*value+2*lead)<.4
    print('Existing radius response, p and p^2 times diagonal log derivative:',p,p*p*value)
print('Limiting p^2 times radius derivative from filled overlaps:',-2*lead)
