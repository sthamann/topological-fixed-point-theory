"""Verified lemma reuse: construct one transport word at a time.

Each accepted lemma is an exact linear MAP identity, stronger than a sampled
constraint check. No target elimination operation is supplied to the proposer.
"""
import argparse
import json
from pathlib import Path
import time
import urllib.request

from .engine import parse_response, save_json
from .mathcore import digest, parity
from .research import Task, Transport


def packet(task, word, library=False):
    active = sorted(set(i for i,s in word))
    names = [task.names[i] for i in active]
    labels = {i:(f'M{i}' if i<4 else task.names[i].upper()) for i in active}
    spec = [{'operator_name':labels[i],'definition':f'D_({task.names[i]}) * j^{task.projections[i][0]} * t0^{task.projections[i][1]}'} for i in active]
    text = ('Normalize ONE word of lattice-field automorphisms. Character parameters c0,c1,c2,c3,r,t are vectors, '
        'whereas M0,M1,M2,M3,R,T are group operators; do not confuse these types. Characters are column vectors in V=F2^8. '
        'D_x D_y=D_(x+y), j^4=t0^2=1, t0 j t0=j^-1, j D_x j^-1=D_(A*x), '
        't0 D_x t0^-1=D_(K*x), A^2=K^2=I, AK=KA. The character D_x is on the LEFT of j^k t0^e. '
        'The supplied word has identity projected action, so it is D_C. Derive its FULL character C '
        'as a linear combination of the listed character variables. Do not set variables to zero or impose extra constraints. '
        'Output {"terms":{variable_name:[coefficient_atoms],...}} with exactly the keys '+json.dumps(names)+'. '
        'Each coefficient is a list of atoms I,A,K,AK summed by XOR. '
        '[] is the zero map. The result must be an identity for ALL input characters, not just a correct zero set. '
        'Operator definitions: '+json.dumps(spec)+'. Word in composition order LEFT to RIGHT: '+
        ' * '.join(labels[i]+('^-1' if s==-1 else '') for i,s in word)+'.')
    if library:
        text += (' Reusable proved NORMAL FORM LEMMA: encode D_c j^k t0^e as (c,k,e). '
            'Multiplication: (c,k,e)*(d,l,f)=(c + A^(k mod2) K^e d, (k+(-1)^e*l) mod4, (e+f) mod2). '
            'Inverse: (c,k,e)^-1=(A^(k mod2) K^e c, ((-1)^(e+1)*k) mod4, e). '
            'Use these elementary identities on the requested word. They are not the target word result.')
    return active,text


def check_lemma(task, word, active, candidate):
    if not isinstance(candidate,dict) or set(candidate)!={'terms'}:
        return {'status':'invalid','error':'Expected terms only'}
    terms = candidate['terms']
    if not isinstance(terms,dict) or set(terms)!={task.names[i] for i in active}:
        return {'status':'invalid','error':'Wrong variable names'}
    expanded = [[] for _ in task.names]
    for i in active:
        expanded[i] = terms[task.names[i]]
    try:
        rows = task.compile({'rows':[expanded]})
    except (ValueError,TypeError) as e:
        return {'status':'invalid','error':str(e)}
    for i in active:
        for bit in range(8):
            vector = 1 << (i*8+bit)
            c,k,e = task.group.word(word,task.values(vector))
            assert k==e==0
            proposed = sum(parity(r & vector) << j for j,r in enumerate(rows))
            if c != proposed:
                return {'status':'counterexample','assignment':{task.names[j]:(1 << bit if j==i else 0) for j in active},
                        'original_character':c,'your_character':proposed,
                        'meaning':'Your formula gives a different character on this basis assignment. Recheck inverse and conjugation actions.'}
    return {'status':'verified','certificate':'Exact equality of linear character maps on every input basis vector',
            'input_basis_vectors':8*len(active)}


def propose(messages, names, model, temperature=0, seed=20261024, think=False):
    atom={'type':'array','maxItems':8,'items':{'type':'string','enum':['I','A','K','AK']}}
    terms={'type':'object','properties':{name:atom for name in names},'required':names,'additionalProperties':False}
    schema={'type':'object','properties':{'terms':terms},'required':['terms'],'additionalProperties':False}
    body = {'model':model,'messages':messages,'stream':False,'think':think,'format':schema,
            'options':{'temperature':temperature,'presence_penalty':0,'repeat_penalty':1,
                       'seed':seed,'num_predict':2048 if think else 512,'num_ctx':8192}}
    req = urllib.request.Request('http://127.0.0.1:11434/api/chat',data=json.dumps(body).encode(),headers={'Content-Type':'application/json'})
    started=time.monotonic(); result=json.load(urllib.request.urlopen(req,timeout=240))
    return {'raw':result['message']['content'],'seconds':time.monotonic()-started,
            'generated_tokens':result.get('eval_count'),'prompt_tokens':result.get('prompt_eval_count'),
            'eval_seconds':result.get('eval_duration',0)/1e9,'done_reason':result.get('done_reason')}


def run(task, outdir, model='qwen3.6:35b', max_attempts=3, feedback=True, temperature=0, library=False, think=False):
    modelinfo=json.load(urllib.request.urlopen('http://127.0.0.1:11434/api/tags'))
    modelinfo=next(x for x in modelinfo['models'] if x['name']==model)
    outdir=Path(outdir); lemmas=outdir/'lemmas'; assembled=[]; records=[]
    for label,word in task.relations:
        active,prompt=packet(task,word,library)
        identity={'prompt':prompt,'model':modelinfo,'source':task.group.source.identity,
                  'checker':digest(Path(__file__).read_text()),'max_attempts':max_attempts,'feedback':feedback,
                  'think':think,'library':library,'temperature':temperature,'presence_penalty':0,
                  'num_predict':2048 if think else 512,'num_ctx':8192}
        key=digest(identity); path=lemmas/f'{key}.json'
        cached=path.exists()
        if cached:
            record=json.loads(path.read_text()); assert record['identity']==identity
            for attempt in record['attempts']:
                if attempt['candidate'] is not None:
                    assert check_lemma(task,word,active,attempt['candidate'])==attempt['report']
        else:
            messages=[{'role':'system','content':'Construct an exact symbolic linear map from the given group relations. Return JSON only.'},
                      {'role':'user','content':prompt}]
            attempts=[]
            for index in range(max_attempts):
                response=propose(messages,[task.names[i] for i in active],model,temperature,think=think)
                try:
                    candidate=parse_response(response['raw']); report=check_lemma(task,word,active,candidate)
                except (ValueError,TypeError) as e:
                    candidate=None; report={'status':'invalid','error':str(e)}
                attempts.append({**response,'candidate':candidate,'report':report})
                record={'identity':identity,'attempts':attempts,'status':report['status']}; save_json(path,record)
                print(json.dumps({'task':task.id,'lemma':label,'attempt':index+1,'report':report}),flush=True)
                if report['status']=='verified' or not feedback:
                    break
                if len(attempts)>=2 and attempts[-1]['candidate']==attempts[-2]['candidate']:
                    record['stopped_reason']='identical_candidate_repeated'; save_json(path,record)
                    break
                messages.extend([{'role':'assistant','content':response['raw']},
                    {'role':'user','content':json.dumps(report)+' Repair the exact character map; JSON only.'}])
        records.append({'label':label,'path':str(path),'reused':cached,'attempts':len(record['attempts']),
                        'status':record['status'],'first_status':record['attempts'][0]['report']['status']})
        if record['status']=='verified':
            row=[[] for _ in task.names]
            for i in active: row[i]=record['attempts'][-1]['candidate']['terms'][task.names[i]]
            assembled.append(row)
        save_json(outdir/f'{task.id}.json',{'task':task.id,'source':task.group.source.identity,'lemmas':records,
            'complete':False,'accepted_block_rows':len(assembled),'host_target_solver_exposed':False})
    complete=len(assembled)==len(task.relations)
    program={'rows':assembled}
    report=task.check(program) if complete else {'status':'incomplete'}
    result={'task':task.id,'source':task.group.source.identity,'lemmas':records,'complete':complete,'program':program,
            'report':report,'host_target_solver_exposed':False,'model':modelinfo}
    save_json(outdir/f'{task.id}.json',result)
    print(json.dumps({'task':task.id,'complete':complete,'report':report}),flush=True)
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tasks',default='reflection-0,reflection-1,rotation-0,dihedral-0')
    p.add_argument('--out',default='results/decomposed')
    p.add_argument('--model',default='qwen3.6:35b')
    p.add_argument('--attempts',type=int,default=3)
    p.add_argument('--no-feedback',action='store_true')
    p.add_argument('--temperature',type=float,default=0)
    p.add_argument('--library',action='store_true')
    p.add_argument('--think',action='store_true')
    args=p.parse_args(); group=Transport()
    for name in args.tasks.split(','):
        kind,s=name.rsplit('-',1)
        run(Task(group,kind,int(s)),args.out,args.model,args.attempts,not args.no_feedback,args.temperature,args.library,args.think)


if __name__=='__main__': main()
