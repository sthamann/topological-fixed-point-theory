"""Original TFPT Stieltjes spin-two spectral pole theorem counterexample.

TFPT 4.5 says F(z)=a+int dmu(t)/(t+z) with a,mu>=0 implies the
only spin-two pole of inverse z*F(z) is massless at z=0.
It does not: zeros on negative real axis produce extra physical poles.
"""
import sympy as sp


def controls():
    z=sp.symbols('z')
    a,m,mu=sp.symbols('a m mu',positive=True)
    F=a+mu/(z+m)
    inverse=sp.factor(1/(z*F))
    zero=sp.solve(sp.together(F),z)[0]
    residue=sp.factor(sp.limit((z-zero)*inverse,z,zero))
    assert sp.simplify(zero+m+mu/a)==0
    assert sp.simplify(residue-mu/(a*(a*m+mu)))==0
    assert residue.is_positive
    concrete=sp.factor(inverse.subs({a:1,m:1,mu:1}))
    split=sp.apart(concrete,z)
    assert concrete==(z+1)/(z*(z+2))
    assert split==sp.Rational(1,2)/z+sp.Rational(1,2)/(z+2)
    # F(i*y) is Herglotz-negative, as in the source paper, but that
    # zero-free upper-half-plane statement misses the physical axis.
    y=sp.symbols('y',positive=True)
    Fy=sp.simplify(sp.im(F.subs(z,sp.I*y)))
    assert sp.factor(Fy+mu*y/(m*m+y*y))==0
    return {
      'status':'EXACT_ORIGINAL_TFPT_SPIN2_STIELTJES_POLE_COUNTEREXAMPLE',
      'original_claim':'positive Stieltjes F forbids any TT poles beyond z=0',
      'valid_positive_stieltjes_F':'F(z)=a+mu/(z+m), a,m,mu>0',
      'zero_of_F':'z=-(m+mu/a) on physical negative real axis',
      'new_TT_propagator_pole':'z=-(m+mu/a)',
      'new_pole_residue':'mu/(a*(a*m+mu)) >0',
      'concrete_example_F':'F(z)=1+1/(z+1)',
      'concrete_TT_propagator':'(z+1)/(z*(z+2))=1/(2*z)+1/(2*(z+2))',
      'concrete_additional_spin2_mass_squared':2,
      'zero_free_in_open_complex_plane_implies_no_poles_on_negative_axis':False,
      'positive_stieltjes_precludes_ghost_sign_for_this_example':True,
      'all_general_spin2_ghosts_excluded_by_this_example':False,
      'original_stieltjes_factor_from_source_Q0_derived':False,
      'original_3plus1D_gravity_spin2_proved':False,
      'F1_solved':False
    }


if __name__=='__main__':
    import json
    print(json.dumps(controls(),indent=2,sort_keys=True))
