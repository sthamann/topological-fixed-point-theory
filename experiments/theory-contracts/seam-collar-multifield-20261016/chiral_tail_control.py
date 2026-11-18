"""Exact tail certificates for the supplied quarter CAR source, not P1 selection."""
from fractions import Fraction as F
import argparse,json
from pathlib import Path


def certificate(cutoff=4):
    if type(cutoff) is not int or cutoff<0:raise ValueError('nonnegative integer cutoff required')
    # l=4 log4, alpha=1/4; q=exp(-l)=1/256, q^alpha=1/4.
    q=F(1,256);first=(F(1,4),F(1,64))
    controls={}
    for p in (1,2):
        total=sum(x**p for x in first)/(1-q**p)
        tail=q**(p*(cutoff+1))*total
        zcut=F(1)
        for n in range(cutoff+1):
            for x in first:zcut*=1+(x*q**n)**p
        controls[str(p)]={'trace_R_power':str(total),'omitted_trace_bound':str(tail),
                          'Z_cutoff_exact':str(zcut),
                          'log_Z_tail_upper_bound':str(tail)}
    return {'source':'declared charged quarter chiral source, positive particle/hole excitation representation',
            'alpha':'1/4','ell':'4 log4','q':str(q),'cutoff':cutoff,
            'heat_operator_norm_cutoff_error':str(first[0]*q**(cutoff+1)),
            'partition_tail_bounds':controls,
            'normalized_choi_vector_squared_error_upper_bound':controls['2']['omitted_trace_bound'],
            'vector_space':'common HS/Choi space, or boundary graph CAR with declared filled input-face reference; not raw Slater vectors in empty doubled Fock space',
            'analytic_basis':'PROOF.txt section8: explicit countable occupation basis and geometric trace bounds',
            'P1_P2_physical_source_selection':False,'absolute_partition_scalar_derived':False}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cutoff',type=int,default=4);p.add_argument('--out',type=Path)
    a=p.parse_args();out=certificate(a.cutoff)
    if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
