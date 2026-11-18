"""Independent diagnostics for R066. Analytic proof covers infinite smooth class."""
from pathlib import Path
import json
import numpy as np
import mpmath as mp
from scipy.linalg import expm, eigh

mp.mp.dps = 100
b = 4 * mp.pi
checks = []
def check(name, error, tolerance):
    error, tolerance = mp.mpf(str(abs(error))), mp.mpf(str(tolerance))
    if not error < tolerance:
        raise RuntimeError(f'{name}: {error} >= {tolerance}')
    checks.append(dict(name=name, error=str(error), tolerance=str(tolerance)))

def theta(a, t=0):
    return sum(mp.exp(-b*(mp.mpf(q*q)/2+(a-mp.mpf('.5'))*q)+1j*t*q)
               for q in range(-30,31))

def moments(a):
    z = theta(a).real
    mean = sum(q*mp.exp(-b*(mp.mpf(q*q)/2+(a-mp.mpf('.5'))*q))
               for q in range(-30,31))/z
    second = sum(q*q*mp.exp(-b*(mp.mpf(q*q)/2+(a-mp.mpf('.5'))*q))
                 for q in range(-30,31))/z
    return mean, second-mean*mean

quarters = [mp.mpf(1)/4, mp.mpf(3)/4]
quarter_data = []
for a in quarters:
    mean, var = moments(a)
    fermi_mean = sum(1/(1+mp.exp(b*(r+a)))-(1 if r<0 else 0)
                     for r in range(-80,81))
    fermi_var = sum(1/(4*mp.cosh(b*(r+a)/2)**2) for r in range(-80,81))
    check(f'original Fermi charge a={a}', fermi_mean-mean, '1e-95')
    check(f'original Fermi variance a={a}', fermi_var-var, '1e-95')
    check(f'exact quarter variance a={a}', var-1/(8*mp.pi), '1e-95')
    check(f'original exterior half a={a}', 1-b*var-mp.mpf('.5'), '1e-95')
    quarter_data.append(dict(a=str(a), mean=str(mean), var=str(var)))

def square(t):
    return sum(mp.exp(-mp.pi*p*p+1j*t*p) for p in range(-30,31))
for t in [mp.mpf(0), mp.mpf('.3'), mp.mpf('1.1'), mp.mpf('.3')+.2j]:
    paired = theta(quarters[0],t)*theta(quarters[1],t)/(theta(quarters[0])*theta(quarters[1]))
    check(f'full pair square theta t={t}', paired-square(t)/square(0), '1e-95')
for v in [mp.mpf('.03'), mp.mpf('.2'), mp.mpf('.4')]:
    phi=theta(quarters[0],2*mp.pi*v)*theta(quarters[1],2*mp.pi*v)/(theta(quarters[0])*theta(quarters[1]))
    original_norm=(mp.jtheta(3,mp.pi*v,mp.exp(-mp.pi))/mp.jtheta(3,0,mp.exp(-mp.pi)))**2
    check(f'original real Quillen norm v={v}',abs(phi)**2-original_norm,'1e-95')
check('original Quillen real Hessian',mp.diff(lambda v:2*mp.log(square(2*mp.pi*v)/square(0)),0,2)+2*mp.pi,'1e-90')
for scale in [mp.mpf('.7'), mp.mpf(1), mp.mpf('1.4')]:
    z = lambda s: sum(mp.exp(-mp.pi*s*p*p) for p in range(-30,31))
    check(f'Poisson modular scale={scale}', z(scale)-z(1/scale)/mp.sqrt(scale), '1e-95')

# Actual finite oscillator traces, not just evaluation of the claimed formula.
cut = 32
ann = np.diag(np.sqrt(np.arange(1,cut,dtype=float)),1).astype(complex)
eps = .23
finite_errors = []
for n in [1,2,3]:
    jn = np.sqrt(n)*ann
    jm = jn.conj().T
    hn = jm@jn
    rho = np.diag(np.exp(-float(b)*n*np.arange(cut)))
    rho /= np.trace(rho)
    f, g = .31+.19j, -.24+.17j
    jf = f*jm+f.conjugate()*jn
    jg = g*jm+g.conjugate()*jn
    h = hn+eps*jf
    ev, vectors = eigh(h)
    weights = np.exp(-float(b)*(ev-ev.min()))
    re = (vectors*weights)@vectors.conj().T/weights.sum()
    measured = np.trace(re@jg).real
    expected = -2*eps*(f*g.conjugate()).real
    check(f'forward perturbed Gibbs mean mode{n}', measured-expected, '4e-13')
    logz = -float(b)*ev.min()+np.log(weights.sum())-np.log(np.diag(rho)[0]**-1)
    check(f'original partition scalar mode{n}', logz-float(b)*eps**2*abs(f)**2, '4e-13')
    k = 1j*eps*f/n
    wk = expm(1j*(k*jm+k.conjugate()*jn))
    transported = wk@rho@wk.conj().T
    check(f'whole oscillator state transport mode{n}', np.linalg.norm(re-transported), '2e-13')
    wg = expm(1j*jg)
    analytic = np.exp(-2j*eps*(f*g.conjugate()).real-.5*n/np.tanh(float(b)*n/2)*abs(g)**2)
    check(f'whole Gibbs Weyl mode{n}', np.trace(re@wg)-analytic, '3e-13')
    pulse = expm(-1j*eps*jf)
    # Independent real time observable, pulse on unchanged original rho.
    t = .41
    gt = g*np.exp(1j*n*t)
    jgt = gt*jm+gt.conjugate()*jn
    direct = np.trace(rho@pulse.conj().T@jgt@pulse)
    comm = n*(f.conjugate()*gt-f*gt.conjugate())
    check(f'causal pulse commutator mode{n}', direct-1j*eps*comm, '3e-13')
    # Three ordered pulses vs exact central Magnus phase, including SK crossphase.
    plus = [.18+.07j, -.09+.12j, .05-.03j]
    minus = [-.04+.13j, .11-.05j]
    def path(fs):
        u = np.eye(cut,dtype=complex)
        central = 0j
        for i, v in enumerate(fs):
            u = expm(-1j*(v*jm+v.conjugate()*jn))@u
            for earlier in fs[:i]:
                central -= .5*n*(v.conjugate()*earlier-v*earlier.conjugate())
        return u, sum(fs), central
    up, fp, bp = path(plus)
    um, fm, bm = path(minus)
    c = n*(fm.conjugate()*fp-fm*fp.conjugate())
    predicted = np.exp(bp-bm+c/2-.5*n/np.tanh(float(b)*n/2)*abs(fm-fp)**2)
    check(f'closed contour unequal paths mode{n}', np.trace(rho@um.conj().T@up)-predicted, '3e-13')
    check(f'closed contour equal paths mode{n}', np.trace(rho@up.conj().T@up)-1, '3e-13')
    finite_errors.append(dict(mode=n, gibbs_mean=float(measured), expected=float(expected)))

# Duhamel integral from independent original thermal ordered mode contractions.
for n in [1,2,3]:
    f, g = mp.mpc('.23','.14'), mp.mpc('-.11','.2')
    integrand = lambda tau: n*(f*mp.conj(g)*mp.exp(n*(tau-b))/(1-mp.exp(-b*n))
                             +mp.conj(f)*g*mp.exp(-n*tau)/(1-mp.exp(-b*n)))
    integral = mp.quad(integrand,[0,b/2,b])
    check(f'original Duhamel ordered integral n={n}', integral-2*mp.re(f*mp.conj(g)), '1e-90')

# Direct finite-dimensional forward one-particle matrix perturbation and
# divided-difference derivative; matrix eigenvalues do not assume a gauge answer.
rr = np.arange(-6,7,dtype=float)
fmat = np.zeros((len(rr),len(rr)),complex)
coeffs = {0:.13,1:.27+.08j,2:-.09+.05j}
for i,r in enumerate(rr):
    for j,s in enumerate(rr):
        d=int(r-s)
        if abs(d) in coeffs:
            fmat[i,j]=coeffs[abs(d)] if d>=0 else np.conj(coeffs[abs(d)])
aa=.25
def fermi_matrix(mat):
    val, vec = eigh(mat)
    occ=1/(1+np.exp(float(b)*val))
    return (vec*occ)@vec.conj().T
step=1e-5
numeric=(fermi_matrix(np.diag(rr+aa)+step*fmat)-fermi_matrix(np.diag(rr+aa)-step*fmat))/(2*step)
c=1/(1+np.exp(float(b)*(rr+aa)))
pred=np.empty_like(fmat)
for i,r in enumerate(rr):
    for j,s in enumerate(rr):
        pred[i,j]=(-float(b)*c[i]*(1-c[i]) if i==j else (c[i]-c[j])/(r-s))*fmat[i,j]
check('actual forward Fermi Frechet response',np.linalg.norm(numeric-pred),'2e-9')

# Genuine smooth disjoint supports; Fourier truncation is only a numerical
# diagnostic. Exact support/vanishing follows from the proof, not this grid.
size=32768
x=np.arange(size)*2*np.pi/size
def bump(center,width):
    dx=(x-center+np.pi)%(2*np.pi)-np.pi
    z=dx/width
    out=np.zeros_like(x)
    inside=abs(z)<1
    out[inside]=np.exp(-1/(1-z[inside]**2))
    return out
f=bump(1.,.35);g=bump(3.,.4)
ff=np.fft.fft(f)/size;gg=np.fft.fft(g)/size
osc=2*np.real(np.vdot(gg[1:size//2],ff[1:size//2]))
mean0f,mean0g=f.mean(),g.mean()
check('proper disjoint support integral',(f*g).mean(),'1e-30')
check('full Fourier Parseval outside cancellation',osc+mean0f*mean0g,'1e-16')
eq=-osc-.5*mean0f*mean0g
check('actual equilibrium exterior half',eq-.5*mean0f*mean0g,'1e-16')
if not eq>0: raise RuntimeError('outside response must be nonzero')
dx=2*np.pi/size
freq=np.fft.fftfreq(size,1/size)
derivative=np.fft.ifft(1j*freq*np.fft.fft(f)).real
check('proper exterior real time response at t0',(derivative*g).mean(),'1e-14')
# Fixed holonomy insertion f=bump derivative has supported primitive k=-epsilon*bump.
check('fixed holonomy derivative mean',derivative.mean(),'1e-14')
check('fixed holonomy exterior response',np.mean((-.23*f)*g),'1e-30')

result=dict(verdict='EXACT_SCOPED_LOCAL_INSERTION_AND_CAUSAL_RESPONSE',
            checks=len(checks),diagnostics=checks,quarters=quarter_data,
            exterior_example=dict(f0=float(mean0f),g0=float(mean0g),
                 normalized_equilibrium_derivative=float(eq),pulse_response_t0=0),
            oscillator_diagnostics=finite_errors,
            scope='Given original free elliptic quantization; all smooth-class proof separate; not primitive P1, physical4D or spin2 closure')
out=Path(__file__).with_name('ELLIPTIC-LOCAL-INSERTION.json')
out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(verdict=result['verdict'],checks=len(checks),output=str(out))))
