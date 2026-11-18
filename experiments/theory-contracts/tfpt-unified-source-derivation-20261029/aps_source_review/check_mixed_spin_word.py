"""R071: mixed original color-spin/CAR word on the retained thermal source.

No physical source selection or boson measure is made. Infinite sea and
oscillator products are analytically retained; finite sums check formulas.
"""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
mp.mp.dps = 70


def theta_one(T, delta, X, shift=0, N=10):
    value = mp.fsum(mp.exp(-T * ((n + shift)**2 / 2 + delta * (n + shift)) + (n + shift) * X)
                    for n in range(-N, N + 1))
    h = mp.re(X) - T * delta
    kp, km = N + 1 + shift, -N - 1 + shift
    rp, rm = mp.exp(-T * (kp + mp.mpf('.5')) + h), mp.exp(T * (km - mp.mpf('.5')) - h)
    tail = mp.exp(-T * kp**2 / 2 + h * kp) / (1 - rp) + mp.exp(-T * km**2 / 2 + h * km) / (1 - rm)
    return value, tail


def theta(T, delta, X, shift=0):
    return mp.fprod(theta_one(T, d, x, shift)[0] for d, x in zip(delta, X))


def prime(T, x, M=8):
    return 2 * mp.sinh(x / 2) * mp.fprod(
        (1 - mp.exp(-T*m + x)) * (1 - mp.exp(-T*m - x)) / (1 - mp.exp(-T*m))**2
        for m in range(1, M + 1))


def original_odd_pair(T, d, x, N=10):
    # Entire negative sea summed exactly, only exponentially small thermal
    # corrections truncated. This is independent of the theta reader.
    sea = mp.exp((d - mp.mpf('.5')) * x) / (1 - mp.exp(-x))
    correction = mp.fsum(
        mp.exp(x * (n + mp.mpf('.5') + d)) / (1 + mp.exp(T * (n + mp.mpf('.5') + d)))
        - mp.exp(-x * (n + mp.mpf('.5') - d)) / (1 + mp.exp(T * (n + mp.mpf('.5') - d)))
        for n in range(N))
    return sea + correction


def main():
    checks = []

    def need(name, condition):
        if not bool(condition):
            raise AssertionError(name)
        checks.append(name)

    def close(name, x, y, tol=mp.mpf('1e-55')):
        need(name, abs(x-y) < tol * max(1, abs(x), abs(y)))

    spec = importlib.util.spec_from_file_location('original_sector_cocycle', ROOT / 'experiments/theory-contracts/canonical-seam-fermionization-20261025/check_same_source_affine_defect.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    m = sp.Matrix(json.loads((ROOT / 'experiments/theory-contracts/v16-operator-pullback-20261014/fold-results.json').read_text())['explicit_Weyl_matrix'])
    e = [sp.eye(8)[:, j] for j in range(8)]
    charges = [m.T * q for q in [e[1]-e[2], e[2]-e[1], e[0]-e[2], e[2]-e[0]]]
    beta = 4 * mp.pi
    T = 8 * beta
    xs = [mp.mpf(3), mp.mpf(2), mp.mpf(1), mp.mpf(0)]
    wF = [1, -1, 0, -1, -1, -1, 0, -1]
    zero_spin = sp.Matrix([sp.Rational(1, 2)]*4 + [-sp.Rational(1, 2)]*4)
    perm = sp.zeros(8)
    for j in range(8):
        perm[j, (j+4) % 8] = 1

    def label(q):
        c = 0 if all(x.q == 1 for x in q) else 1
        return tuple(int(x-sp.Rational(c, 2)) for x in q), c

    def group_word(ps):
        acc, phase = ((0,)*8, 0), 1+0j
        for p in ps:
            lb = label(p)
            phase *= old.eps_ext(acc, lb)
            acc = old.add_labels(acc, lb)
        need('balanced original charge product', acc == ((0,)*8, 0))
        return phase

    records = []
    for alpha in [mp.mpf(1)/4, mp.mpf(3)/4]:
        ds = [(j+alpha)/8-mp.mpf('.5') for j in range(8)]
        zNS = theta(T, ds, [0]*8)
        zR = theta(T, ds, [0]*8, mp.mpf('.5'))
        d_exact = sp.Matrix([(j + sp.Rational(str(alpha))) / 8 - sp.Rational(1, 2) for j in range(8)])
        need('original NS/R affine energy intertwiner', perm.T * (zero_spin+d_exact) == d_exact and zero_spin.dot(zero_spin)/2 + d_exact.dot(zero_spin) == 0)
        close('exact affine NS/R partition identity', zNS, zR)

        def color_theta(X, shift=0, parity=None):
            terms = []
            for k in range(3):
                xx = [x+2j*mp.pi*k*w/3 for x, w in zip(X, wF)]
                ordinary = theta(T, ds, xx, shift)
                if parity is not None:
                    marked = theta(T, ds, [x+1j*mp.pi for x in xx], shift)
                    ordinary = (ordinary + parity*marked)/2
                terms.append(ordinary)
            return mp.fsum(terms)/3

        zNS0 = color_theta([0]*8)
        zR0 = color_theta([0]*8, mp.mpf('.5'))
        zE80 = color_theta([0]*8, parity=1) + color_theta([0]*8, mp.mpf('.5'), parity=1)
        need('center restricted NS/R states have positive weights', mp.re(zNS0)>0 and mp.re(zR0)>0)
        need('affine NS/R identification fails fixed center intertwining', perm.T*sp.Matrix(wF) != sp.Matrix(wF))
        need('center restricted partitions differ', abs(zNS0-zR0)>mp.mpf('1e-12'))
        a_bound = (sp.Rational(1, 20)+sp.Rational(1, 20)**3)/(1-sp.Rational(1, 20)**4)
        need('exact rational bound for whole source partition', 1/(1-a_bound)<sp.Rational(53, 50))
        osc = mp.fprod((1-mp.exp(-T*n))**-8 for n in range(1, 9))
        zfermion = mp.fprod((1+mp.exp(-beta*(n+alpha)))*(1+mp.exp(-beta*(n+1-alpha))) for n in range(50))
        close('same whole one-source partition after actual eight-fold', zNS*osc, zfermion)
        for j in range(8):
            close('literal frequency transport', 8*(mp.mpf('.5')+ds[j]), j+alpha)
            tj, bound = theta_one(T, ds[j], 1)
            need('explicit theta tail bound', bound < mp.mpf('1e-1000'))
            two = mp.exp(ds[j]) * tj / theta_one(T, ds[j], 0)[0] / prime(T, 1)
            close('original entire sea equals theta odd two-point', two, original_odd_pair(T, ds[j], 1))
        for k, q in enumerate(charges):
            need('actual missing color root has h1', q.dot(q) == 2)
            need('actual missing color root belongs to spinor coset s', label(q)[1] == 1 and sum(label(q)[0]) % 2 == 0)
            for j, v in enumerate(e):
                ps = [q, v, -v, -q]
                # The old charge cocycle and actual adjoints fix the paired
                # coefficient; no arbitrary square-root phase is inserted.
                adj = old.eps_ext(label(v), label(-v)).conjugate() * old.eps_ext(label(q), label(-q)).conjugate()
                need('original paired adjoint cocycle is exactly one', group_word(ps) * adj == 1)
                rho = mp.mpf(str(q[j]))
                X = [mp.fsum(mp.mpf(str(p[a]))*x for p, x in zip(ps, xs)) for a in range(8)]
                thermal = mp.exp(mp.fdot(ds, X)) * theta(T, ds, X) / zNS
                cross = prime(T, 1)**2 / prime(T, 2)**2
                reduced = thermal * cross**rho / prime(T, 3)**2 / prime(T, 1)
                direct = thermal
                for a in range(4):
                    for b in range(a+1, 4):
                        direct *= prime(T, xs[a]-xs[b])**mp.mpf(str(ps[a].dot(ps[b])))
                close('whole mixed word agrees with cross-ratio reduction', direct, reduced)
                # Actual rectangular trace cyclicity: moving the color
                # spinor through the thermal cut changes NS to R. The
                # shifted insertion remains in a real ordered strip.
                rotated_ps = ps[1:]+ps[:1]
                rotated_xs = xs[1:]+[xs[0]-T]
                need('thermal cyclic word keeps exact original adjoint scalar', group_word(rotated_ps)*adj == 1)
                Y = [mp.fsum(mp.mpf(str(p[a]))*x for p, x in zip(rotated_ps, rotated_xs)) for a in range(8)]
                cyclic_kernel = mp.mpf(1)
                for a in range(4):
                    for b in range(a+1, 4):
                        cyclic_kernel *= prime(T, rotated_xs[a]-rotated_xs[b])**mp.mpf(str(rotated_ps[a].dot(rotated_ps[b])))
                cycled = mp.exp(mp.fdot(ds, Y))*theta(T, ds, Y, mp.mpf('.5'))/zR*cyclic_kernel
                close('complete thermal NS-to-R sewing of actual mixed word', direct*zNS, cycled*zR)
                need('nonzero mixed answer', abs(direct) > mp.mpf('1e-10'))
                need('full local monodromy remains minus one', sp.exp(2*sp.pi*sp.I*q[j]) == -1)
                if k == 0 and j in [0, 1]:
                    records.append({'alpha':str(alpha),'root':'M^T(e2-e3)','odd_channel':j,'rho':str(rho),'mixed_word_value':mp.nstr(direct,25)})
                if k == 0 and j == 0:
                    projected = mp.exp(mp.fdot(ds, X))*color_theta(X)/zNS0*cross**rho/prime(T, 3)**2/prime(T, 1)
                    need('whole original color restriction retains actual mixed word', mp.re(projected)>0 and abs(mp.im(projected))<mp.mpf('1e-55'))
                    cycled_color = mp.exp(mp.fdot(ds, Y))*color_theta(Y, mp.mpf('.5'))/zR0*cyclic_kernel
                    close('center-restricted full cyclic sewing retains absolute block normalizers', projected*zNS0, cycled_color*zR0)
                    records[-1]['center_restricted_mixed_value'] = mp.nstr(mp.re(projected),25)
        # A single cyclic transfer by a spinor charge changes the base
        # charge coset. These exact Fourier-period signs separate the blocks.
        close('NS charge theta is 2pi-periodic', theta_one(T, ds[0], 2j*mp.pi)[0], theta_one(T, ds[0], 0)[0])
        close('R charge theta is 2pi-antiperiodic', theta_one(T, ds[0], 2j*mp.pi, mp.mpf('.5'))[0], -theta_one(T, ds[0], 0, mp.mpf('.5'))[0])
        records.append({'alpha':str(alpha),'original_NS_Z0':mp.nstr(zNS*osc,25),'whole_R_Z0':mp.nstr(zR*osc,25),'conditional_full_parent_NS_weight':mp.nstr(zNS/(zNS+zR),25),
                        'center_NS_Z0':mp.nstr(mp.re(zNS0*osc),25),'center_R_Z0':mp.nstr(mp.re(zR0*osc),25),'center_E8_Z0':mp.nstr(mp.re(zE80*osc),25),
                        'conditional_full_center_NS_weight':mp.nstr(mp.re(zNS0/(zNS0+zR0)),25)})

    out = {'result':'R071','verdict':'GIVEN_SOURCE_MIXED_WORD_AND_SECTOR_SEWING_EXACT_PHYSICAL_F1_OPEN',
           'checks_passed':len(checks),'checks':checks,'examples':records,
           'scope':'Actual fixed four color roots and all eight original folded odd channels, both given quarters, beta4pi. Full original NS source and its original cocycle. Positive full-parent Gibbs extension only conditional on the entire B(H) algebra with same given H and beta. No source/measure/96-matter/spacetime selection.'}
    (HERE/'MIXED-SPIN-WORD.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'checks_passed':len(checks),'verdict':out['verdict'],'examples':records}))


if __name__ == '__main__':
    main()
