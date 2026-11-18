"""Learn symbolic character maps by choosing discriminating experiments.

This is finite exact program synthesis, not neural training. The hypothesis
grammar is given; the target coefficients are not read from the oracle matrix.
The oracle returns actual original word evaluations only. Whole-map equality is
verified after learning, separately from the experiment selection loop.
"""
import argparse
import json
from pathlib import Path
import random

from .decomposed import packet, check_lemma
from .engine import save_json
from .mathcore import digest
from .research import Transport, Task, act, xor, linear
from .frontier import forgetful


ATOMS=['I','A','K','AK']


def hypotheses(group):
    # Collapse different expressions denoting the very same source operator.
    maps={}
    for mask in range(16):
        names=[name for i,name in enumerate(ATOMS) if mask >> i&1]
        columns=tuple(xor(group.atoms[name][j] for name in names) for j in range(8))
        if columns not in maps or len(names)<len(maps[columns]): maps[columns]=names
    return [{'columns':list(columns),'atoms':names} for columns,names in maps.items()]


def learn(task,word,mode,seed):
    rng=random.Random(seed); active,prompt=packet(task,word)
    grammar=hypotheses(task.group); terms={}; queries=[]
    for variable in active:
        candidates=list(grammar); asked=set(); base_cursor=0
        while len(candidates)>1:
            if mode=='active':
                # Maximize worst-case elimination; candidates are equally weighted.
                def score(x):
                    groups={}
                    for h in candidates:
                        y=act(h['columns'],x);groups[y]=groups.get(y,0)+1
                    return (len(candidates)-max(groups.values()),len(groups),-x)
                x=max(range(1,256),key=score)
            elif mode=='basis':
                x=1 << base_cursor;base_cursor+=1
            else:
                x=rng.choice([v for v in range(1,256) if v not in asked])
            assert x not in asked
            asked.add(x)
            assignment=x << (8*variable)
            actual,k,e=task.group.word(word,task.values(assignment))
            assert k==e==0
            before=len(candidates)
            candidates=[h for h in candidates if act(h['columns'],x)==actual]
            if not candidates:raise ArithmeticError('Hypothesis grammar does not contain this original map')
            queries.append({'variable':task.names[variable],'input_character':x,'observed_output':actual,
                            'hypotheses_before':before,'hypotheses_after':len(candidates)})
        terms[task.names[variable]]=candidates[0]['atoms']
    candidate={'terms':terms}; report=check_lemma(task,word,active,candidate)
    assert report['status']=='verified'
    return {'candidate':candidate,'queries':queries,'query_count':len(queries),'report':report,
            'grammar_semantic_size':len(grammar),'mode':mode}


def run(mode,outdir,seed=20261024):
    group=Transport(); library={}; records=[]; fresh_queries=0
    for name in ['reflection-0','reflection-1','rotation-0','dihedral-0']:
        kind,s=name.rsplit('-',1); task=Task(group,kind,int(s)); rows=[]; lemma_refs=[]
        for label,word in task.relations:
            active,prompt=packet(task,word)
            key=digest({'prompt':prompt,'source':group.source.identity,'mode':mode,'seed':seed})
            reused=key in library
            if not reused:
                library[key]=learn(task,word,mode,seed+len(library));fresh_queries+=library[key]['query_count']
            result=library[key]
            assert check_lemma(task,word,active,result['candidate'])['status']=='verified'
            row=[[] for _ in task.names]
            for i in active:row[i]=result['candidate']['terms'][task.names[i]]
            rows.append(row);lemma_refs.append({'label':label,'key':key,'reused':reused})
        program={'rows':rows};report=task.check(program);assert report['status']=='verified'
        consumed=forgetful(task,task.compile(program))
        records.append({'task':name,'program':program,'report':report,'lemmas':lemma_refs,
                        'forgetful_classification_from_learned_program':consumed})
    result={'method':mode,'seed':seed,'status':'verified','tasks':records,'learned_lemmas':library,
            'distinct_lemmas':len(library),'fresh_oracle_queries':fresh_queries,
            'independent_basis_checks_for_distinct_lemmas':sum(x['report']['input_basis_vectors'] for x in library.values()),
            'query_model':'Single-variable evaluation of one original group word; full oracle rows not exposed',
            'input_grammar':'XOR sums of I,A,K,AK; finite and source justified',
            'neural_model_used':False,'physical_completion':False,'source':group.source.sources,
            'source_id':group.source.identity,'code_sha256':digest(Path(__file__).read_text())}
    save_json(Path(outdir)/f'{mode}-{seed}.json',result)
    return result


def restore(path):
    """New process: restore maps, verify source, compose, use actual consumer."""
    data=json.loads(Path(path).read_text()); group=Transport()
    assert data['source_id']==group.source.identity
    consumed=[]
    for record in data['tasks']:
        kind,s=record['task'].rsplit('-',1); task=Task(group,kind,int(s))
        rows=[]
        for (label,word),ref in zip(task.relations,record['lemmas'],strict=True):
            assert ref['label']==label
            active,prompt=packet(task,word); lemma=data['learned_lemmas'][ref['key']]
            assert check_lemma(task,word,active,lemma['candidate'])['status']=='verified'
            row=[[] for _ in task.names]
            for i in active:row[i]=lemma['candidate']['terms'][task.names[i]]
            rows.append(row)
        program={'rows':rows};assert program==record['program']
        assert task.check(program)==record['report']
        consumer=forgetful(task,task.compile(program))
        assert consumer==record['forgetful_classification_from_learned_program']
        consumed.append({'task':task.id,'underlying_Gamma_classes':consumer['underlying_Gamma_classes']})
    return {'restored_from':str(path),'source_verified':True,'new_learning_queries':0,
            'fresh_certificate_verification':True,'actual_consumer_outputs':consumed}


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/active-search')
    p.add_argument('--random-repetitions',type=int,default=20)
    p.add_argument('--restore')
    args=p.parse_args();summary=[]
    if args.restore:
        result=restore(args.restore);save_json(Path(args.out)/'restore.json',result);print(json.dumps(result,indent=2));return
    for mode in ['basis','active']:
        r=run(mode,args.out);summary.append({k:r[k] for k in ['method','distinct_lemmas','fresh_oracle_queries','status']})
    random_counts=[]
    for i in range(args.random_repetitions):
        r=run('random',args.out,20261024+i);random_counts.append(r['fresh_oracle_queries'])
    summary.append({'method':'random','repetitions':len(random_counts),'min_queries':min(random_counts),
                    'max_queries':max(random_counts),'mean_queries':sum(random_counts)/len(random_counts),'status':'all verified'})
    save_json(Path(args.out)/'summary.json',summary);print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
