"""Exact arithmetic PROGRAM controls; no primitive TFPT source identification."""
from pathlib import Path
import json
import sympy as s
import native_constructors as n


def main():
    u,y=s.symbols('jU jY');b,q=s.symbols('body quartic',real=True)
    chart={'modes':2,'mode_chart_id':'ARITHMETIC_RETAINED_TWO_MODE_CHART'}
    p=n.quartic_source_program(chart,{'expression':b},{'expression':q},(u,y))
    op=n.nil_action(p)
    bindings={k:{'role':'arithmetic_control','description':'explicit arithmetic dictionary; not physical P1 marks, state or time'}
              for k in ('P','cap','reference','absolute','time','U','Y')}
    d={'P':s.eye(4),'cap':s.eye(4),'reference_scalar':s.Integer(3),'outer_scalar':s.Integer(2),'bindings':bindings}
    got=n.source_jets(op,d)
    checks={
        'body_and_nonGaussian_quartic_construct_a_full_kernel':op['matrix']==s.diag(b,b*(1+u),b*(1+y),b*(q+(1+u)*(1+y))),
        'actual_four_Berezin_trace_values':all(s.simplify(got[k]-v)==0 for k,v in
            {'Z0':2*b*(q+4)/3,'ZU':4*b/3,'ZY':4*b/3,'ZUY':2*b/3}.items()),
        'connected_from_same_unormalized_source':s.simplify(got['normalized']['connected_UY']-q/(q+4)**2)==0,
    }
    # Boundary projector and cap mix the even states; do not commute them through K.
    v=s.Matrix([1,0,0,s.I])/s.sqrt(2);P=v*v.H
    cap=s.Matrix([[2,0,0,1+s.I],[0,1,0,0],[0,0,1,0],[3-s.I,0,0,4]])
    mixed=n.source_jets(op,{**d,'P':P,'cap':cap})
    checks['whole_outer_P_and_noncommuting_cap_retained']=all(s.simplify(mixed[k]-s.Rational(2,3)*s.trace(cap*P*Q*P))==0 for k,Q in mixed['ordered_source_forms'].items())
    # Direct kernel chart includes zero-body sources; no inverse/body log is used.
    zero=n.even_kernel({'modes':1,'mode_chart_id':'ARITHMETIC_ONE_MODE','source_symbols':(u,y),
                        'coefficients':{3:1},'hypotheses':['zero-body arithmetic source']})
    checks['zero_body_kernel_constructor_supported']=zero['matrix']==s.diag(0,1)
    sewed=n.sew(op,op)
    checks['exact_original_Berezin_sewing_is_ordered_product']=sewed['matrix']==op['matrix']*op['matrix']
    # A separate explicit endpoint-word program gives real source derivatives.
    base=n.even_kernel({'modes':2,'mode_chart_id':chart['mode_chart_id'],'source_symbols':(u,y),
                        'coefficients':n.native_module().kernel(s.diag(2,3,5,7))})
    U=s.Matrix([[0,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,0]])
    Y=s.diag(1,2,3,4)
    marked=n.endpoint_marks(base,{'U':U,'Y':Y,'placement':'U_initial_Y_final',
                     'time_binding':{'role':'arithmetic_control','description':'specified initial and final endpoint'}})
    ep=n.source_jets(marked,d)
    checks['endpoint_whole_word_derivative_no_semigroup_assumption']=ep['ordered_source_forms']['ZUY']==Y*base['matrix']*U
    try:n.source_jets(op,{'bindings':bindings})
    except n.NativeInputMissing:checks['no_default_P_cap_ref_or_absolute_scalar']=True
    else:checks['no_default_P_cap_ref_or_absolute_scalar']=False
    # These are ORIGINAL scalar leaves used by a NEW declared source program.
    actual_program=n.quartic_source_program(chart,n.scalar_c3(),n.scalar_c6_ratio(),(u,y))
    actual_values=n.source_jets(n.nil_action(actual_program),d)
    result={'scope':'arithmetic program-constructor controls only; no original U/Y or time identity claimed',
            'checks':{k:bool(v) for k,v in checks.items()},'pass':all(checks.values()),
            'symbolic_program_values':{k:str(got[k]) for k in ('Z0','ZU','ZY','ZUY')},
            'original_scalar_seed_candidate_values':{k:str(actual_values[k]) for k in ('Z0','ZU','ZY','ZUY')},
            'candidate_connected':str(actual_values['normalized']['connected_UY']),
            'new_source_hypotheses':actual_values['hypotheses'],
            'native_version':n.version()}
    Path(__file__).with_name('native-constructor-controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['pass']:raise SystemExit(1)

if __name__=='__main__':main()
