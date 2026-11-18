"""Second path: literal source field words and SymPy's GF(2) matrix rank."""
import json
from sympy.polys.matrices import DomainMatrix
from sympy import GF

from .engine import save_json
from .research import Task, Transport


def main():
    group=Transport(); reports=[]
    for name in ['reflection-0','reflection-1','rotation-0','dihedral-0']:
        kind,s=name.rsplit('-',1);task=Task(group,kind,int(s))
        rank=DomainMatrix.from_list([[(r >> j)&1 for j in range(task.n)] for r in task.rows],GF(2)).rank()
        assert rank==task.solution['rank']
        columns=0
        for coordinate in range(task.n):
            v=1 << coordinate;values=task.values(v);bits=[]
            for label,word in task.relations:
                for n in group.source.units:
                    out=n;phase=1
                    # Every variable has order dividing 4, so its inverse is
                    # its third power. No normal-form inverse/multiplication
                    # implementation is used along this verification path.
                    for variable,sign in reversed(word):
                        for _ in range(1 if sign==1 else 3):
                            p,out=group.field(values[variable],out);phase*=p
                    assert out==n
                    bits.append(int(phase==-1))
            assert bits==task.residual(v)
            columns+=1
        reports.append({'task':name,'sympy_GF2_rank':rank,'literal_source_field_columns':columns,
                        'field_generators_per_relation':8,'relators':len(task.relations)})
    result={'status':'passed','source':group.source.identity,'independent_paths':reports,
            'scope':'Every linear coefficient independently checked in the original cocycle field action, plus separate GF2 rank implementation'}
    save_json('results/independent-check.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':main()
