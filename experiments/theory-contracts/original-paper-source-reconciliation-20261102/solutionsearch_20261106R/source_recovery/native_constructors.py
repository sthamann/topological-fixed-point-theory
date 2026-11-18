"""Typed PROGRAM constructors backed by the original coherent-kernel algebra.

Candidates are weighted kernel constructions, not numerical Higgs fields.
No physical H, state, time, projector, cap or reference is supplied by default.
The output has no physical-completion flag: the injected validator decides.
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import importlib.util
import sys
import sympy as sp

ROOT=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
NATIVE=ROOT/'experiments/theory-contracts/seam-collar-multifield-20261016'


class NativeInputMissing(ValueError):
    pass


@lru_cache(None)
def native_module():
    # Load the exact original functions; preserve unrelated module namespaces.
    saved={k:sys.modules.get(k) for k in ('joint_source','source_sewing')}
    try:
        for name in saved:
            spec=importlib.util.spec_from_file_location(name,NATIVE/(name+'.py'))
            mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
        spec=importlib.util.spec_from_file_location('_r_original_kernel_functional',NATIVE/'kernel_functional.py')
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        return mod
    finally:
        for name,old in saved.items():
            if old is None:sys.modules.pop(name,None)
            else:sys.modules[name]=old


@lru_cache(None)
def version():
    files=[NATIVE/'kernel_functional.py',NATIVE/'joint_source.py',NATIVE/'source_sewing.py',
           ROOT/'verification/v330_qgamb_admissible_measure.py',
           ROOT/'documents/archive--tfpt-45--source_extracts--01_boundary_kernel_source/sections/002-section-intrinsic-closed-branch-and-canonical-rigid-obje.tex']
    return hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()


def _metadata(value):
    return list(value.get('hypotheses',[])),list(value.get('provenance',[]))


def _mode(value):
    n=value.get('modes')
    if type(n) is not int or not 1<=n<=5:
        raise NativeInputMissing('the existing retained kernel interface requires an explicitly bound original1..5-mode chart')
    if not value.get('mode_chart_id'):
        raise NativeInputMissing('missing actual retained mode chart identification')
    return n


def _symbols(value):
    symbols=value.get('source_symbols')
    if not isinstance(symbols,(list,tuple)) or len(symbols)!=2 or not all(isinstance(x,sp.Symbol) for x in symbols) or symbols[0]==symbols[1]:
        raise NativeInputMissing('two distinct formal source symbols must be supplied explicitly')
    return tuple(symbols)


def _poly(value,n,zero_body=False):
    poly={}
    for mask,coefficient in value.items():
        if type(mask) is not int or not 0<=mask<(1<<(2*n)) or mask.bit_count()%2:
            raise ValueError('kernel/action masks must be even in the ORIGINAL2n Grassmann ordering')
        if zero_body and mask==0 and coefficient!=0:
            raise ValueError('nilpotent action must have zero body; absolute body is a separate searched program')
        poly[mask]=sp.sympify(coefficient)
    return poly


def even_kernel(program):
    """Actual forward constructor, including zero-body kernels and all vertices."""
    n=_mode(program);symbols=_symbols(program)
    poly=_poly(program['coefficients'],n)
    T=native_module().matrix_from_kernel(poly,n)
    # Check parity with the actual original constructor's checker.
    native_module().kernel(T)
    h,p=_metadata(program)
    return {'matrix':T,'modes':n,'mode_chart_id':program['mode_chart_id'],
            'source_symbols':symbols,'hypotheses':h,'provenance':p+[{'operation':'original.matrix_from_kernel','version':version()}],
            'scope':'candidate whole weighted even kernel on the explicitly bound retained source fibre'}


def nil_action(program):
    """Searchable scalar BODY plus arbitrary even nilpotent interactions."""
    n=_mode(program);poly=_poly(program['nil_action'],n,True)
    body=sp.sympify(program['body'])
    k=native_module()
    weighted=k.scale(k.exponential(k.scale(poly,-1),2*n),body)
    return even_kernel({**program,'coefficients':weighted})


def sew(left,right):
    """Exact original Berezin sewing, preserving absolute kernel coefficients."""
    if left['mode_chart_id']!=right['mode_chart_id'] or left['modes']!=right['modes'] or left['source_symbols']!=right['source_symbols']:
        raise ValueError('sewing needs the SAME actual chart, source symbols and labelled endpoint identification')
    k=native_module();poly=k.sew(left['matrix'],right['matrix'])
    return even_kernel({**left,'coefficients':poly,
                        'hypotheses':left['hypotheses']+right['hypotheses'],
                        'provenance':left['provenance']+right['provenance']+[{'operation':'original.sew','version':version()}]})


def adjoint(operator):
    return {**operator,'matrix':operator['matrix'].H,
            'provenance':operator['provenance']+[{'operation':'actual matrix adjoint'}]}


def endpoint_marks(operator,mark_pair):
    """Exact WHOLE endpoint word, never an inferred semigroup or chosen H."""
    if mark_pair.get('placement')!='U_initial_Y_final' or not mark_pair.get('time_binding'):
        raise NativeInputMissing('endpoint marks need explicit original/hypothesis order and time binding')
    T=operator['matrix'];U=sp.Matrix(mark_pair['U']);Y=sp.Matrix(mark_pair['Y'])
    if U.shape!=T.shape or Y.shape!=T.shape:
        raise ValueError('the actual endpoint marks must share the original whole source space')
    u,y=operator['source_symbols'];I=sp.eye(T.rows)
    marked=(I+y*Y)*T*(I+u*U)
    native_module().kernel(marked)
    h,p=_metadata(mark_pair)
    return {**operator,'matrix':marked,
            'hypotheses':operator['hypotheses']+h,
            'provenance':operator['provenance']+p+[{'operation':'whole endpoint source (I+jY Y) K (I+jU U)',
                                                  'time_binding':mark_pair['time_binding']}],
            'endpoint_marks':{'U':U,'Y':Y,'time_binding':mark_pair['time_binding']}}


def source_jets(operator,dictionary):
    """Compute four WHOLE unnormalized values by exact marked-source derivatives.

    No power of a retained transfer is substituted for a whole ordered word.
    Source coefficients belong to the proposed whole marked kernel K(jU,jY).
    Their mark/time identity remains an explicit validator obligation.
    """
    required=('P','cap','reference_scalar','outer_scalar','bindings')
    if any(key not in dictionary for key in required):
        raise NativeInputMissing('source_jets requires explicit P, cap, absolute/reference scalars and mark/time bindings')
    bindings=dictionary['bindings']
    for key in ('P','cap','reference','absolute','time','U','Y'):
        row=bindings.get(key)
        if not isinstance(row,dict) or row.get('role') not in ('original','hypothesis','arithmetic_control') or not row.get('description'):
            raise NativeInputMissing('each P/cap/reference/absolute/time/U/Y dictionary binding must declare its actual or hypothesis role')
    T=sp.Matrix(operator['matrix']);P=sp.Matrix(dictionary['P']);cap=sp.Matrix(dictionary['cap'])
    if P.shape!=T.shape or cap.shape!=T.shape:
        raise ValueError('whole P and cap must act on the SAME retained source space')
    ref=sp.sympify(dictionary['reference_scalar']);outer=sp.sympify(dictionary['outer_scalar'])
    if ref==0:raise ValueError('actual reference quotient has a zero denominator')
    u,y=operator['source_symbols'];zero={u:0,y:0}
    if any(value.has(u,y) for value in (ref,outer,*P,*cap)):
        raise ValueError('source-dependent P/cap/rest must be INSIDE the whole marked kernel before differentiation')
    forms={'Z0':T.subs(zero),'ZU':T.diff(u).subs(zero),'ZY':T.diff(y).subs(zero),
           'ZUY':T.diff(u,y).subs(zero)}
    k=native_module();values={}
    for key,Q in forms.items():
        # The original coherent trace performs the actual Berezin contraction.
        values[key]=sp.simplify(outer*k.coherent_trace(cap*P*Q*P)/ref)
        if values[key].has(sp.zoo,sp.nan,sp.oo,-sp.oo):
            raise ValueError('candidate marked kernel has a singular or nonfinite zero-source jet')
    z=values['Z0'];normalized={}
    if z!=0:
        normalized={key:sp.simplify(values[zkey]/z) for key,zkey in [('U','ZU'),('Y','ZY'),('UY','ZUY')]}
        normalized['connected_UY']=sp.simplify(normalized['UY']-normalized['U']*normalized['Y'])
    hypotheses=list(operator['hypotheses'])+[f'dictionary:{key}: {row["description"]}' for key,row in bindings.items() if row['role']=='hypothesis']
    return {'matrix':T,'modes':operator['modes'],'mode_chart_id':operator['mode_chart_id'],
            'source_symbols':operator['source_symbols'],'ordered_source_forms':forms,
            **values,'normalized':normalized,'bindings':bindings,'hypotheses':hypotheses,
            'dictionary':{**dictionary,'P':P,'cap':cap,'reference_scalar':ref,'outer_scalar':outer},
            'provenance':operator['provenance']+[{'operation':'original.coherent_trace of whole P Q P with cap','version':version()}],
            'scope':'computed common source values of a candidate marked kernel; original field/time/source identity is a validator obligation'}


def scalar_c3():
    return {'expression':1/(8*sp.pi),'hypotheses':[],
            'provenance':[{'source':'Boundary primitive c3=1/(8pi)','role':'actual original scalar'}]}


def scalar_c6_ratio():
    return {'expression':sp.Rational(64,729),'hypotheses':[],
            'provenance':[{'source':'verification/v330_qgamb_admissible_measure.py symmetric_transfer','role':'actual C6 ratio, not a full-field weight'}]}


def scalar_add(a,b):
    return {'expression':sp.simplify(a['expression']+b['expression']),
            'hypotheses':a.get('hypotheses',[])+b.get('hypotheses',[]),
            'provenance':a.get('provenance',[])+b.get('provenance',[])}


def scalar_multiply(a,b):
    return {**scalar_add(a,b),'expression':sp.simplify(a['expression']*b['expression'])}


def scalar_one():
    return {'expression':sp.Integer(1),'hypotheses':[],
            'provenance':[{'role':'operator coefficient-ring unit; not a physical normalization'}]}


def scalar_i():
    return {**scalar_one(),'expression':sp.I}


def scalar_negate(a):
    return {**a,'expression':-a['expression']}


def scalar_divide(a,b):
    if b['expression']==0:raise ValueError('candidate coefficient division by zero')
    return {**scalar_add(a,b),'expression':sp.simplify(a['expression']/b['expression'])}


def source_scalar_U(symbol_pair):
    return {'expression':_symbols({'source_symbols':symbol_pair})[0],
            'hypotheses':[],'provenance':[{'role':'formal marked-source variable U; mark identity is not assumed'}]}


def source_scalar_Y(symbol_pair):
    return {'expression':_symbols({'source_symbols':symbol_pair})[1],'hypotheses':[],
            'provenance':[{'role':'formal marked-source variable Y; mark identity is not assumed'}]}


def monomial(chart,mask,coefficient,symbol_pair):
    n=_mode(chart);symbols=_symbols({'source_symbols':symbol_pair})
    poly=_poly({mask:coefficient['expression']},n)
    h,p=_metadata(coefficient)
    ch,cp=_metadata(chart)
    return {**chart,'coefficients':poly,'source_symbols':symbols,'hypotheses':ch+h,'provenance':cp+p}


def kernel_add(left,right):
    if left['mode_chart_id']!=right['mode_chart_id'] or left['modes']!=right['modes'] or left['source_symbols']!=right['source_symbols']:
        raise ValueError('kernel coefficient programs need the same source chart')
    coefficients=dict(left['coefficients'])
    for mask,value in right['coefficients'].items():
        coefficients[mask]=sp.simplify(coefficients.get(mask,0)+value)
    return {**left,'coefficients':coefficients,'hypotheses':left['hypotheses']+right['hypotheses'],
            'provenance':left['provenance']+right['provenance']}


def action_from_polynomial(poly,body):
    _poly(poly['coefficients'],_mode(poly),True)
    return {**poly,'nil_action':poly['coefficients'],'body':body['expression'],
            'hypotheses':poly['hypotheses']+body.get('hypotheses',[]),
            'provenance':poly['provenance']+body.get('provenance',[])}


def quartic_source_program(chart,body,quartic,symbol_pair):
    """Explicit SEARCH HYPOTHESIS, not the original physical U/Y dictionary.

    A concrete nonGaussian source-program family whose four values are really
    evaluated. Its coherent-bilinear action marks are stated hypotheses.
    """
    n=_mode(chart)
    if n<2:raise ValueError('this quartic program template needs two actually retained modes')
    symbols=_symbols({'source_symbols':symbol_pair});u,y=symbols
    action={((1<<j)|(1<<(n+j))):-sp.Integer(1) for j in range(n)}
    action[(1<<0)|(1<<(n+0))]-=u
    action[(1<<1)|(1<<(n+1))]-=y
    action[(1<<0)|(1<<1)|(1<<(n+0))|(1<<(n+1))]=quartic['expression']
    return {'modes':n,'mode_chart_id':chart['mode_chart_id'],'body':body['expression'],
            'nil_action':action,'source_symbols':symbols,
            'hypotheses':chart.get('hypotheses',[])+body.get('hypotheses',[])+quartic.get('hypotheses',[])+['candidate scalar body and quartic SOURCE relation',
                          'jU marks first coherent bilinear in the kernel action, jY second; no primitive U/Y or endpoint-occupation identity is assumed',
                          'whole ordered-time identity must be checked; no semigroup is inferred'],
            'provenance':chart.get('provenance',[])+body.get('provenance',[])+quartic.get('provenance',[])}


def registry_operations():
    specs=[
        ('native.even_kernel',('EvenGrassmannKernel',),'EvenSourceOperator',even_kernel),
        ('native.nil_action',('KernelActionProgram',),'EvenSourceOperator',nil_action),
        ('native.sew',('EvenSourceOperator','EvenSourceOperator'),'EvenSourceOperator',sew),
        ('native.adjoint',('EvenSourceOperator',),'EvenSourceOperator',adjoint),
        ('native.endpoint_marks',('EvenSourceOperator','SourceMarkPair'),'EvenSourceOperator',endpoint_marks),
        ('native.source_jets',('EvenSourceOperator','WholeSourceDictionary'),'WholeSourceCandidate',source_jets),
        ('native.c3',(),'SourceScalar',scalar_c3),
        ('native.c6_ratio',(),'SourceScalar',scalar_c6_ratio),
        ('native.scalar_add',('SourceScalar','SourceScalar'),'SourceScalar',scalar_add),
        ('native.scalar_multiply',('SourceScalar','SourceScalar'),'SourceScalar',scalar_multiply),
        ('native.scalar_one',(),'SourceScalar',scalar_one),
        ('native.scalar_i',(),'SourceScalar',scalar_i),
        ('native.scalar_negate',('SourceScalar',),'SourceScalar',scalar_negate),
        ('native.scalar_divide',('SourceScalar','SourceScalar'),'SourceScalar',scalar_divide),
        ('native.source_scalar_U',('SourceSymbolPair',),'SourceScalar',source_scalar_U),
        ('native.source_scalar_Y',('SourceSymbolPair',),'SourceScalar',source_scalar_Y),
        ('native.monomial',('ModeChart','EvenGrassmannMask','SourceScalar','SourceSymbolPair'),'EvenGrassmannKernel',monomial),
        ('native.kernel_add',('EvenGrassmannKernel','EvenGrassmannKernel'),'EvenGrassmannKernel',kernel_add),
        ('native.action_from_polynomial',('EvenGrassmannKernel','SourceScalar'),'KernelActionProgram',action_from_polynomial),
        ('native.quartic_source_program',('ModeChart','SourceScalar','SourceScalar','SourceSymbolPair'),'KernelActionProgram',quartic_source_program),
    ]
    return [dict(id=key,input_types=inputs,output_type=output,execute=func,version=version()) for key,inputs,output,func in specs]
