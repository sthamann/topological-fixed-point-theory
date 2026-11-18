"""Original E8 CX regular-product zero modes, with ORIGINAL lattice cocycles.

Native source computation on the entire existing conformal-weight-one block
(240 root states plus8 Cartan states). No physical H/time/P/state is chosen.
The result tests the proposed DIRECT coordinate identification only.
"""
from fractions import Fraction
from pathlib import Path
import json
from tfpt_explorer.current_block_geometry import _actual_current_source
from tfpt_explorer.source_program import _field


def _get(i,a):
    if type(i) is not int or not 0<=i<5 or type(a) is not int or not 0<=a<4:
        raise ValueError('use existing carrier i=0..4 and family a=0..3')
    s=_actual_current_source();g=s['chevalley']
    ci,cp=next(iter(_field({'kind':'C','a':a})[1].items()))
    xi,xp=next(iter(_field({'kind':'X','i':i,'a':a})[1].items()))
    c,x=tuple(g.roots[ci]),tuple(g.roots[xi])
    v=tuple(u+w for u,w in zip(c,x))
    dot=lambda p,q:sum(Fraction(u*w,4) for u,w in zip(p,q))
    assert dot(c,x)==0 and dot(v,v)==4
    kappa=Fraction(cp*xp*g.eps(ci,xi))
    return g,ci,xi,c,x,v,kappa,dot


def cx_zero_mode(i,a):
    """Sparse o(:C_a X_ia:), including original OPE and input phases.

    Acting on an original root state |p> at L0=1, the h=2 zero mode selects
    (v,p)=-2 and returns kappa*epsilon(v,p)|p+v>. Bimultiplicativity uses
    only ORIGINAL root cocycles eps(C,p)eps(X,p), never a new extension.
    The8 Cartan states have no output at this grade by charge/energy.
    """
    g,ci,xi,c,x,v,kappa,dot=_get(i,a)
    action={};opposite={}
    nc,nx=g.opp[ci],g.opp[xi]
    # Actual vertex adjoint, computed from the SAME original 2-cocycle.
    adj=g.eps(nc,ci)*g.eps(nc,xi)*g.eps(nx,ci)*g.eps(nx,xi)
    assert adj==1
    for pi,praw in enumerate(g.roots):
        p=tuple(praw)
        if dot(v,p)==-2:
            q=tuple(u+w for u,w in zip(p,v))
            assert q in g.ridx
            action[pi]=(g.ridx[q],kappa*g.eps(ci,pi)*g.eps(xi,pi))
        if dot(v,p)==2:
            q=tuple(u-w for u,w in zip(p,v))
            assert q in g.ridx
            opposite[pi]=(g.ridx[q],kappa*adj*g.eps(nc,pi)*g.eps(nx,pi))
    assert len(action)==len(opposite)
    for p,(q,k) in action.items():
        assert opposite[q]==(p,k), 'actual lattice adjoint must be matrix adjoint'
    return {'source':g,'charge_doubled':v,'kappa':kappa,'adjoint_sign':adj,
            'C':c,'X':x,'action':action,'adjoint_action':opposite,
            'root_dimension':len(g.roots),'cartan_dimension':8}


def coordinate_commutator_packet(i,a):
    """A native nonzero source matrix element for proposed Re/Im coordinates."""
    z=cx_zero_mode(i,a);g=z['source'];T=z['action'];Td=z['adjoint_action']
    p=g.ridx[tuple(-u for u in z['C'])]
    target,coef=T[p]
    assert tuple(g.roots[target])==z['X'] and p not in Td
    assert Td[target]==(p,coef) and coef*coef==1
    # Full sparse diagonal of [T,Tdagger], including every original root leg.
    comm={}
    for j in range(len(g.roots)):
        tt=Td.get(j)
        td=T.get(j)
        v=Fraction(0)
        if tt:
            middle,c=tt;back,d=T[middle];assert back==j;v+=c*d
        if td:
            middle,c=td;back,d=Td[middle];assert back==j;v-=c*d
        if v:comm[j]=v
    assert comm[p]==-1 and all(v in (-1,1) for v in comm.values())
    return {'i':i,'a':a,'whole_original_grade_dimension':len(g.roots)+8,
            'CX_charge_doubled':list(z['charge_doubled']),
            'original_regular_OPE_coefficient':str(z['kappa']),
            'vertex_adjoint_sign':z['adjoint_sign'],
            'T_rank_on_grade_one':len(T),'commutator_nonzero_rank_on_grade_one':len(comm),
            'witness_incoming_doubled_charge':list(g.roots[p]),
            'witness_T_outgoing_doubled_charge':list(g.roots[target]),
            'witness_T_coefficient':str(coef),
            'witness_Tdagger_annihilates_incoming':True,
            'witness_T_Tdagger_commutator':'-1',
            'coordinate_definition':'X_Re=(T+Tdagger)/2, X_Im=(T-Tdagger)/(2i)',
            'witness_coordinate_commutator':'-I/2',
            'T_sparse':[{'from':p,'to':q,'coefficient':str(c)} for p,(q,c) in sorted(T.items())],
            'Tdagger_sparse':[{'from':p,'to':q,'coefficient':str(c)} for p,(q,c) in sorted(Td.items())],
            'physical_Higgs_identified':False}


def main():
    rows=[coordinate_commutator_packet(i,a) for a in range(4) for i in range(5)]
    frames=[]
    for a in range(4):
        weak=[r for r in rows if r['a']==a and r['i'] in (3,4)]
        frames.append({'coordinate_family_frame':a,'weak_components':[3,4],
                       'transport_rule':'conjugate original fields, Fuchs connection, line projector and time together',
                       'weak_coordinate_witnesses':[r['witness_coordinate_commutator'] for r in weak],
                       'new_physical_line_chosen':False})
    out={'scope':'DIRECT CX conformal zero-mode coordinate test in original E8 source; not a TFPT no-go',
         'whole_original_grade_one_dimension':248,'original_CX_pairs':20,
         'actual_all_twenty_cocycles_and_adjoints_checked':True,
         'original_four_family_frames':frames,'packets':rows,
         'conclusion':'raw CX zero-mode real/imaginary fields fail commuting-coordinate identity before physical P',
         'physical_UY_value':None,'physical_HB':None,'physical_time_selected':False}
    path=Path(__file__).with_name('native-cx-zero-mode-result.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='packets'},indent=2))


if __name__=='__main__':main()
