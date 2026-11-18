"""Symbolic controls only; physical source and reconstruction are not derived."""
import hashlib
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def main():
    checks = []
    def require(name, condition):
        passed = bool(condition)
        checks.append({"name": name, "passed": passed})
        if not passed:
            raise ValueError(name)

    x, e, n, z = sp.symbols('x epsilon n z', real=True)
    f, u, potential = sp.Function('f')(x), sp.Function('u')(x), sp.Function('V')(x)
    scale = sp.exp(-f / 2)
    df_action = scale * (-sp.I * sp.diff(scale * sp.exp(f / 2) * u, x)
                         + potential * scale * sp.exp(f / 2) * u)
    d_action = -sp.I * sp.diff(u, x) + potential * u
    require('symmetric source action maps arbitrary zero solution',
            sp.simplify(df_action - scale * d_action) == 0)
    require('nonzero spectral equation contains interior weight exp(f)',
            sp.simplify(scale * z * sp.exp(f) * u - z * sp.exp(f / 2) * u) == 0)
    require('fibre conformal law plus L2 density has left exponent -1/2',
            sp.simplify(n / 2 - (n + 1) / 2 + sp.Rational(1, 2)) == 0)
    require('fibre conformal law plus L2 density has right exponent -1/2',
            sp.simplify((n - 1) / 2 - n / 2 + sp.Rational(1, 2)) == 0)
    h, R, grad2 = sp.symbols('h R grad_h_squared', real=True)
    coefficient_after_ibp = 2 * (n - 1) * (n - 2) - (n - 1) * (n - 2)
    require('integrated conformal gradient coefficient in 4D is +6',
            coefficient_after_ibp.subs(n, 4) == 6)
    density = sp.exp(2 * e * h) * (R + 6 * e**2 * grad2)
    require('first curvature variation is 2 h R',
            sp.simplify(sp.diff(density, e).subs(e, 0) - 2 * h * R) == 0)
    require('second curvature variation is 4 h^2 R + 12 grad_h_squared',
            sp.simplify(sp.diff(density, e, 2).subs(e, 0)
                        - (4 * h**2 * R + 12 * grad2)) == 0)
    require('Gram weight from two solution factors is exp(f)',
            sp.simplify(sp.exp(f / 2) * sp.exp(f / 2) - sp.exp(f)) == 0)
    require('first Gram profile variation is h',
            sp.diff(sp.exp(e * h), e).subs(e, 0) == h)
    weighted_coefficient = 2 * (n - 1) * 4 - (n - 1) * (n - 2)
    require('chi^2-weighted conformal curvature in seed frame has coefficient +18',
            weighted_coefficient.subs(n, 4) == 18)
    sigma, grad_sigma2, lap_sigma = sp.symbols('sigma grad_sigma_squared lap_sigma', real=True)
    local_L = sp.exp(4 * sigma) * (R + 18 * grad_sigma2)
    local_EL = sp.diff(local_L, sigma) - 36 * sp.exp(4 * sigma) * (lap_sigma + 4 * grad_sigma2)
    require('weighted curvature subterm EL is R-9 lap_sigma-18 grad_sigma_squared',
            sp.simplify(local_EL / (4 * sp.exp(4 * sigma))
                        - (R - 9 * lap_sigma - 18 * grad_sigma2)) == 0)
    c = sp.symbols('c', real=True)
    require('both constant scale factors jointly give exp(4c)',
            sp.simplify(sp.exp(2 * c) * sp.exp(2 * c) - sp.exp(4 * c)) == 0)
    s = sp.symbols('s', real=True)
    r1, r0, regular0, alpha, area, kappa, chi0, L = sp.symbols(
        'r1 r0 regular0 alpha area kappa chi0 Lambda', real=True)
    scale_diff = sp.exp(2*s*c) - 1
    require('even-reference scale multiplier vanishes at baseline',
            scale_diff.subs(c, 0) == 0)
    require('relative Einstein residue retains the reference subtraction',
            sp.simplify(sp.residue(scale_diff * r1/(s-1), s, 1)
                        -(sp.exp(2*c)-1)*r1) == 0)
    require('regular order-zero zeta has zero relative constant-scale finite part',
            sp.limit(scale_diff * regular0, s, 0) == 0)
    require('simple zero-pole contributes 2c times its residue',
            sp.limit(scale_diff * (r0/s + regular0), s, 0) == 2*c*r0)
    require('weighted relative residue is exp(4c)-exp(2c)',
            sp.simplify(-6*L**2*sp.exp(2*c)*(sp.exp(2*c)-1)*r1
                        +6*L**2*(sp.exp(4*c)-sp.exp(2*c))*r1) == 0)
    W, Wref = sp.Function('W')(alpha), sp.Function('Wref')(alpha)
    relative_moment = W-Wref-alpha*(W.subs(alpha,1)-Wref.subs(alpha,1))
    expected_entropy_difference = -sp.diff(W,alpha).subs(alpha,1)+W.subs(alpha,1) \
        +sp.diff(Wref,alpha).subs(alpha,1)-Wref.subs(alpha,1)
    require('normalized replica quotient gives entropy difference',
            sp.simplify(-sp.diff(relative_moment,alpha).subs(alpha,1)
                        -expected_entropy_difference) == 0)
    cone_W = kappa*(alpha*regular0+4*sp.pi*(1-alpha)*area)
    cone_S = -sp.diff(cone_W-alpha*cone_W.subs(alpha,1),alpha).subs(alpha,1)
    require('conditional original EH and entropy coefficients coincide',
            sp.simplify(cone_S.subs(kappa,chi0**2/(4*sp.pi**2))
                        -chi0**2*area/sp.pi) == 0)
    require('conditional modular coefficient equals inverse reduced Planck mass squared',
            sp.simplify(2*sp.pi/(4*sp.pi*kappa)-1/(2*kappa)) == 0)
    # Minimal exact symbol control: opposite nonzero boundary covectors.
    cn = sp.I * sp.Matrix([[1, 0], [0, -1]])
    ct = sp.I * sp.Matrix([[0, 1], [1, 0]])
    bpos = sp.I * cn * ct
    bneg = sp.I * cn * (-ct)
    require('Calderon principal symbol is a Hermitian involution',
            bpos.H == bpos and bpos*bpos == sp.eye(2))
    require('opposite boundary covectors require opposite nonzero involutions',
            bneg == -bpos and bneg != bpos)
    ppos, pneg = (sp.eye(2)+bpos)/2, (sp.eye(2)+bneg)/2
    require('Calderon symbol projectors are complementary for opposite covectors',
            ppos*ppos == ppos and pneg*pneg == pneg
            and ppos+pneg == sp.eye(2) and ppos*pneg == sp.zeros(2))
    require('local matrix matching one covector fails on its opposite',
            bpos-bneg == 2*bpos and bpos-bneg != sp.zeros(2))
    # Finite symbol/domain algebra of the analytic double, not a TFPT model.
    F = sp.I * cn
    Q = bpos * F
    require('normal interchange anticommutes with Calderon polarization',
            F.H == F and F*F == sp.eye(2)
            and F*bpos+bpos*F == sp.zeros(2))
    require('boundary order-four product has only plus/minus i phases',
            Q*Q == -sp.eye(2) and Q**4 == sp.eye(2)
            and set(Q.eigenvals()) == {sp.I, -sp.I})
    I2, Z2 = sp.eye(2), sp.zeros(2)
    K = sp.BlockMatrix([[Z2, sp.I*I2], [-sp.I*I2, Z2]]).as_explicit()
    double_symbol = sp.diag(ct, -ct)
    glued_trace = I2.col_join(cn)
    domain_constraint = (-cn).row_join(I2)
    require('pair reflection is a Hermitian unit involution and Dirac odd',
            K.H == K and K*K == sp.eye(4)
            and K*double_symbol+double_symbol*K == sp.zeros(4))
    require('correct pair reflection preserves clutching domain',
            domain_constraint*K*glued_trace == sp.zeros(2))
    naive_lift = sp.BlockMatrix([[Z2, sp.I*cn], [sp.I*cn, Z2]]).as_explicit()
    require('normal-times-swap shorthand fails clutching domain',
            domain_constraint*naive_lift*glued_trace != sp.zeros(2))
    gamma = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Gdouble = sp.diag(gamma, -gamma)
    U0 = sp.I*Gdouble*K
    require('double chirality preserves domain and anticommutes with reflection',
            domain_constraint*Gdouble*glued_trace == sp.zeros(2)
            and Gdouble*K+K*Gdouble == sp.zeros(4))
    require('twisted even lift is Hermitian involution commuting with Dirac',
            U0.H == U0 and U0*U0 == sp.eye(4)
            and U0*double_symbol-double_symbol*U0 == sp.zeros(4))
    Eof, Ebase = sp.Function('E')(e), sp.symbols('E_baseline')
    require('frozen reference zero at baseline does not kill derivative',
            sp.diff(Eof-Ebase,e) == sp.diff(Eof,e)
            and sp.diff(Eof-Eof,e) == 0)
    d1, Zn, Z1 = sp.symbols('d1 Zn Z1', positive=True)
    require('tensor-compatible fixed denominator cancels normalized replica moments',
            sp.simplify((Zn/d1**n)/(Z1/d1)**n-Zn/Z1**n) == 0)
    r_eff, a = sp.symbols('r_eff a', real=True)
    a2_density = r_eff*(-sp.Rational(1,4)+sp.Rational(1,6))/(4*sp.pi)**2
    require('ordinary Dirac curvature residue is minus rank over 192 pi squared',
            sp.simplify(a2_density+r_eff/(192*sp.pi**2)) == 0)
    require('asserted Einstein residue selects effective rank eight',
            sp.solve(sp.Eq(a2_density,-1/(24*sp.pi**2)),r_eff) == [8])
    require('rank eight gives asserted reduced Planck mass',
            sp.simplify((-12*chi0**2*a2_density).subs(r_eff,8)
                        -chi0**2/(2*sp.pi**2)) == 0)
    # Twirling is a conditional expectation, not a rank division.
    clock = sp.diag(1, sp.I, -1, -sp.I)
    entries = sp.symbols('a0:16')
    A = sp.Matrix(4,4,entries)
    clock_CE = sum((clock**j*A*clock**(-j) for j in range(4)),sp.zeros(4))/4
    require('clock conditional expectation preserves ordinary trace',
            sp.simplify(sp.trace(clock_CE)-sp.trace(A)) == 0)
    require('clock conditional expectation does not divide rank',
            sum((clock**j*sp.eye(4)*clock**(-j) for j in range(4)),sp.zeros(4))/4
            == sp.eye(4))
    # Joint source coefficients: trace independent ratio and a repair control.
    H, p = sp.symbols('H p', real=True)
    trace_E = -R*r_eff/4-H*p
    trace_E2 = R**2*r_eff/16+R*H*p/2+H**2*p
    A2 = sp.expand((trace_E+R*r_eff/6)/(16*sp.pi**2))
    A4 = sp.expand((180*trace_E2+60*R*trace_E+5*R**2*r_eff)
                   /(360*16*sp.pi**2))
    require('common displayed endomorphism gives a2 Higgs coefficient',
            A2.coeff(H) == -p/(16*sp.pi**2))
    require('common displayed endomorphism gives a4 Higgs quartic',
            A4.coeff(H,2) == p/(32*sp.pi**2))
    require('common displayed endomorphism gives a4 mixed curvature Higgs',
            A4.coeff(H).coeff(R) == p/(192*sp.pi**2))
    require('mixed-to-mass ratio is trace independent one twelfth',
            sp.cancel(A4.coeff(H).coeff(R)/(-A2.coeff(H))) == sp.Rational(1,12))
    require('printed mixed-to-mass ratio is six times common operator ratio',
            (sp.Rational(1,96)/sp.Rational(1,48))/sp.Rational(1,12) == 6)
    equations = [sp.Eq(A2.coeff(H),-1/(48*sp.pi**2)),
                 sp.Eq(A4.coeff(H,2),1/(16*sp.pi**2)),
                 sp.Eq(A4.coeff(H).coeff(R),1/(96*sp.pi**2))]
    require('one common trace parameter cannot match all printed scalar coefficients',
            sp.solve(equations,[p]) == [])
    printed_A2 = -R/(24*sp.pi**2)-H/(48*sp.pi**2)
    require('conditional ordinary rank48 rank2 contraction differs by exactly six at a2',
            sp.simplify(A2.subs({r_eff:48,p:2})-6*printed_A2) == 0)
    require('uniform sqrt6 operator scale changes residue but not regular zero value',
            sp.residue(6**(-s)*r1/(s-1),s,1) == r1/6
            and sp.limit(6**(-s)*regular0,s,0) == regular0)
    require('common spectral unit rescaling preserves chi squared residue',
            sp.simplify((sp.sqrt(6)*chi0)**2*r1/6-chi0**2*r1) == 0)
    norm_omega2 = sp.symbols('norm_omega_squared', nonnegative=True)
    p_source = sp.Rational(1,3)
    quartic_with_unitary_connection = p_source/(32*sp.pi**2)-norm_omega2/(192*sp.pi**2)
    require('positive missing quartic exceeds maximal affine unitary connection coefficient',
            sp.simplify(1/(16*sp.pi**2)-quartic_with_unitary_connection)
            == (10+norm_omega2)/(192*sp.pi**2))
    source_paths = [
        '_archive/tfpt-45/source_extracts/01_boundary_kernel_source.tex',
        '_archive/tfpt-45/source_extracts/02_carrier_source.tex',
        '_archive/tfpt-45/source_extracts/05_metrology_source.tex',
        'experiments/theory-contracts/source-conformal-kernel-selection-20260921/README.md',
        'experiments/theory-contracts/original-source-functional-audit-20261013/METRIC_SOURCE_AUDIT.txt',
        'experiments/theory-contracts/source-operator-origin-20260920/ERGEBNIS.md',
    ]
    sources = []
    for relative in source_paths:
        data = (ROOT / relative).read_bytes()
        sources.append({"path": relative, "sha256": hashlib.sha256(data).hexdigest()})
    result = {
        'verdict': 'PARTIAL', 'symbolic_checks': checks,
        'passed': sum(c['passed'] for c in checks), 'total': len(checks),
        'analytic_review': 'independent scoped reviews incorporated for sections1-17, including double/frozen reference/replica/rank formulas, standard-connection joint coefficient mismatch and the uniform-scale dictionary alternative and affine connection sign; physical source not validated', 'physical_gates_closed': [],
        'not_executed': ['PDE spectral chart', 'unique continuation', 'inverse reconstruction',
                         'TFPT source selection', 'replica functional', 'full verification suite'],
        'sources': sources,
    }
    (HERE / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f"Symbolic controls: {result['passed']}/{result['total']}; verdict PARTIAL; physical gates closed: 0")


if __name__ == '__main__':
    main()
