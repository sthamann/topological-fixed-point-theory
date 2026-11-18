"""Model-derived transport constraints, checked against original group words.

The model constructs a symbolic presentation. The oracle cannot be called as a
solver by the model; it only accepts or produces a separating counterexample.
All mathematical statements concern explicitly marked finite transport objects.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random
import time
import urllib.request

from .engine import save_json, parse_response
from .mathcore import Source, exact_gf2, digest, parity


def xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def act(columns, vector):
    return xor(c for i, c in enumerate(columns) if vector >> i & 1)


def linear(rows, n):
    return exact_gf2(rows, [0] * len(rows), n)


class Transport:
    def __init__(self, source=None):
        self.source = source or Source()
        self.A = self.dual(self.source.matrices['J'])
        self.K = self.dual(self.source.matrices['T'])
        self.atoms = {'I': [1 << i for i in range(8)], 'A': self.A, 'K': self.K,
                      'AK': [act(self.A, x) for x in self.K]}
        for i in range(256):
            assert act(self.A, act(self.A, i)) == i
            assert act(self.K, act(self.K, i)) == i
            assert act(self.A, act(self.K, i)) == act(self.K, act(self.A, i))

    @staticmethod
    def dual(matrix):
        dual = matrix.inv().T
        return [sum((int(dual[i, j]) % 2) << i for i in range(8)) for j in range(8)]

    def multiply(self, x, y):
        c, k, e = x
        d, l, f = y
        if e:
            d = act(self.K, d)
        if k % 2:
            d = act(self.A, d)
        return c ^ d, (k + (-l if e else l)) % 4, e ^ f

    def inverse(self, x):
        c, k, e = x
        k_inv = (k if e else -k) % 4
        d = act(self.K, c) if e else c
        d = act(self.A, d) if k_inv % 2 else d
        return d, k_inv, e

    def word(self, word, values):
        out = (0, 0, 0)
        for index, sign in word:
            value = values[index]
            out = self.multiply(out, value if sign == 1 else self.inverse(value))
        return out

    def field(self, x, n):
        c, k, e = x
        phase = 1
        if e:
            v, n = self.source.apply('T', 0, n)
            phase *= v
        for _ in range(k):
            v, n = self.source.apply('J', 0, n)
            phase *= v
        phase *= (-1 if parity(c & sum((v % 2) << i for i, v in enumerate(n))) else 1)
        return phase, n

    def audit(self):
        # Exact relations on generators determine equality of lattice extension
        # automorphisms. Random non-generator vectors are additional controls.
        for n in self.source.units:
            for word in [('J',)*4, ('T',)*2, ('T','J','T','J')]:
                assert self.source.word(word, 0, 0, n) == (1, n)
        for name, generator in [('A', (0,1,0)), ('K',(0,0,1))]:
            for i in range(8):
                conjugate = self.multiply(self.multiply(generator, (1 << i,0,0)), self.inverse(generator))
                assert conjugate == (self.atoms[name][i],0,0)
                for n in self.source.units:
                    s, v = self.field(self.inverse(generator), n)
                    z, v = self.field((1 << i,0,0), v); s *= z
                    z, v = self.field(generator, v); s *= z
                    assert (s,v) == self.field(conjugate,n)
        rng = random.Random(20261024)
        for _ in range(120):
            x, y = [(rng.randrange(256),rng.randrange(4),rng.randrange(2)) for _ in range(2)]
            n = tuple(rng.randrange(-2,3) for _ in range(8))
            sy, v = self.field(y,n); sx, v = self.field(x,v)
            assert (sx*sy,v) == self.field(self.multiply(x,y),n)
        return {'generator_relations': 'exact', 'dual_character_actions': 'exact on all basis generators',
                'original_cocycle_field_product_controls':120, 'source_id': self.source.identity}


class Task:
    def __init__(self, group, kind='reflection', s=0):
        self.group, self.kind, self.s = group, kind, s
        self.names = ['c0','c1','c2','c3']
        self.projections = [(1,0)]*4
        self.relations = [('meridian product', [(i,1) for i in range(4)])]
        if kind in ('rotation','dihedral'):
            r = len(self.names); self.names.append('r'); self.projections.append((1,0))
            self.relations.append(('rotation fourth power',[(r,1)]*4))
            for i in range(4):
                self.relations.append((f'rotation equivariance {i}',[(r,1),(i,1),(r,-1),((i+1)%4,-1)]))
        if kind in ('reflection','dihedral'):
            t = len(self.names); self.names.append('t'); self.projections.append((s,1))
            self.relations.append(('reflection square',[(t,1)]*2))
            for i in range(4):
                self.relations.append((f'reflection equivariance {i}',[(t,1),(i,1),(t,-1),((s-i)%4,1)]))
        if kind == 'dihedral':
            self.relations.append(('dihedral compatibility',[(t,1),(r,1),(t,-1),(r,1)]))
        self.n = 8*len(self.names)
        assert not any(self.residual(0))
        columns = [self.residual(1 << i) for i in range(self.n)]
        self.rows = [sum(columns[j][i] << j for j in range(self.n)) for i in range(8*len(self.relations))]
        self.solution = linear(self.rows,self.n)
        self.id = f'{kind}-{s}'

    def values(self, vector):
        return [(vector >> (8*i) & 255, k, e) for i, (k,e) in enumerate(self.projections)]

    def residual(self, vector):
        result = []
        for label, word in self.relations:
            c,k,e = self.group.word(word,self.values(vector))
            assert k == e == 0, (label,'base action mismatch')
            result.extend((c >> i) & 1 for i in range(8))
        return result

    def compile(self, program):
        if not isinstance(program,dict) or set(program) != {'rows'}:
            raise ValueError('Expected exactly {"rows": [...]}')
        blocks = program['rows']
        if not isinstance(blocks,list) or not 1 <= len(blocks) <= 16:
            raise ValueError('Use 1..16 block rows')
        rows = []
        for blockrow in blocks:
            if not isinstance(blockrow,list) or len(blockrow) != len(self.names):
                raise ValueError(f'Each block row has {len(self.names)} blocks')
            columns = []
            for block in blockrow:
                if not isinstance(block,list) or len(block)>8 or any(a not in self.group.atoms for a in block):
                    raise ValueError('Each block lists atoms I,A,K,AK; [] is zero; sum is XOR')
                columns.extend(xor(self.group.atoms[a][i] for a in block) for i in range(8))
            rows.extend(sum(((v >> i)&1) << j for j,v in enumerate(columns)) for i in range(8))
        return rows

    def check(self, program):
        try:
            rows = self.compile(program)
        except (ValueError,TypeError) as e:
            return {'status':'invalid','error':str(e)}
        proposed = linear(rows,self.n)
        for vector in proposed['null_basis']:
            residual = self.residual(vector)
            if any(residual):
                j = next(i for i,b in enumerate(residual) if b)
                return {'status':'underconstrained','assignment':dict(zip(self.names,[x[0] for x in self.values(vector)])),
                        'failing_original_relation':self.relations[j//8][0],
                        'nonzero_character':sum(residual[(j//8)*8+i] << i for i in range(8)),
                        'meaning':'Your equations accept this assignment but the original group relation fails.'}
        for vector in self.solution['null_basis']:
            if any(parity(r & vector) for r in rows):
                return {'status':'overconstrained','assignment':dict(zip(self.names,[x[0] for x in self.values(vector)])),
                        'meaning':'All original group relations hold for this assignment, but your equations reject it.'}
        return {'status':'verified','rank':proposed['rank'],'dimension':self.n-proposed['rank'],
                'proof':'Equality of kernels by both basis containments over F2', 'host_target_solver_exposed':False}

    def packet(self):
        action = {'A':self.group.A,'K':self.group.K}
        text = ('Construct ALL necessary and sufficient homogeneous F2-linear equations for the character variables. '
                'You must derive the equations, not solve a supplied matrix. E8 character space V=F2^8. '
                'There are fixed ORIGINAL lattice-field lifts j and t0 with j^4=t0^2=1 and t0*j*t0=j^-1. '
                'D_c are sign characters, D_c D_d=D_(c+d), j D_c j^-1=D_(A c), t0 D_c t0^-1=D_(K c). '
                'A^2=K^2=I and AK=KA. Coefficients of character columns are on the LEFT. '
                'Every meridian is M_i=D_ci j. Ordered marked meridians satisfy M0*M1*M2*M3=1. ')
        if self.kind in ('rotation','dihedral'):
            text += ('An unknown rotation lift R=D_r j must satisfy R^4=1 and '
                     'R M_i R^-1=M_(i+1 mod4) for ALL i. ')
        if self.kind in ('reflection','dihedral'):
            text += (f'An unknown reflection lift T=D_t j^{self.s} t0 must satisfy T^2=1 and '
                     f'T M_i T^-1=M_(({self.s}-i) mod4)^-1 for ALL i. ')
        if self.kind == 'dihedral':
            text += 'Also impose T R T^-1=R^-1. '
        text += ('Give block equations in variables '+','.join(self.names)+'. Each block row means sum(block_operator * variable)=0 in V. '
                 'A block is a list of atoms I,A,K,AK summed by XOR, [] means zero. '
                 'Every block row must contain exactly '+str(len(self.names))+' blocks in that variable order. '
                 'Return only the JSON object with rows. '
                 'At most 16 block rows. Do not invent or omit a condition. Actual action columns as 8-bit integer masks: '+json.dumps(action))
        return text

    def classification(self):
        # Simultaneous conjugation by D_g, in the fixed projected marking.
        gauge = []
        for i in range(8):
            g = (1 << i,0,0)
            values = [self.group.multiply(self.group.multiply(g,(0,k,e)),g)[0] for k,e in self.projections]
            vector = sum(c << (8*j) for j,c in enumerate(values))
            assert not any(self.residual(vector))
            gauge.append(vector)
        rank_g = linear(gauge,self.n)['rank']
        dimension = self.n-self.solution['rank']
        return {'variables':self.n,'constraint_rank':self.solution['rank'],'solution_dimension':dimension,
                'solution_count':1 << dimension,'gauge_rank':rank_g,'gauge_orbit_size':1 << rank_g,
                'character_gauge_stabilizer_order':1 << (8-rank_g),'quotient_dimension':dimension-rank_g,
                'character_gauge_classes':1 << (dimension-rank_g),
                'equivalence':'simultaneous conjugation by D_g only, fixed projected operators; no mapping class quotient'}


def model_call(messages, think, variables, seed=20261024):
    block = {'type':'array','items':{'type':'string','enum':['I','A','K','AK']},'maxItems':8}
    row = {'type':'array','items':block,'minItems':variables,'maxItems':variables}
    schema = {'type':'object','properties':{'rows':{'type':'array','items':row,'minItems':1,'maxItems':16}},
              'required':['rows'],'additionalProperties':False}
    body = {'model':'qwen3.6:35b','messages':messages,'stream':False,'think':think,'format':schema,
            'options':{'temperature':0,'seed':seed,'num_predict':8192}}
    start = time.monotonic()
    req = urllib.request.Request('http://127.0.0.1:11434/api/chat',data=json.dumps(body).encode(),headers={'Content-Type':'application/json'})
    reply = json.load(urllib.request.urlopen(req,timeout=240))
    return {'raw':reply['message']['content'],'thinking':reply['message'].get('thinking',''),
            'seconds':time.monotonic()-start,'prompt_tokens':reply.get('prompt_eval_count'),
            'generated_tokens':reply.get('eval_count'),'done_reason':reply.get('done_reason')}


def run(task, strategy, outdir, seed=20261024):
    messages = [{'role':'system','content':'You construct exact mathematical presentations from definitions. Return the requested JSON.'},
                {'role':'user','content':task.packet()}]
    identity = {'task':task.id,'packet':task.packet(),'source':task.group.source.identity,'strategy':strategy,
                'seed':seed,'script':digest(Path(__file__).read_text()),'checker_access':'only counterexample or accept',
                'model':json.load(urllib.request.urlopen('http://127.0.0.1:11434/api/tags'))['models']}
    # Keep only the tested model metadata, not other inventory.
    identity['model'] = next(x for x in identity['model'] if x['name']=='qwen3.6:35b')
    path = Path(outdir)/f'{task.id}-{strategy}-{seed}.json'
    if path.exists():
        old = json.loads(path.read_text())
        if old['identity'] != identity:
            raise ValueError('Changed run identity: use a new output directory')
        for attempt in old['attempts']:
            if attempt.get('program') is not None:
                assert task.check(attempt['program']) == attempt['report']
        return old
    attempts = []
    think = strategy != 'direct'
    for index in range(1 if strategy in ('direct','reasoning') else 3):
        common = Path(outdir)/f'{task.id}-reasoning-{seed}.json'
        if index == 0 and strategy in ('self_review','counterexample') and common.exists():
            shared = json.loads(common.read_text())
            for key in ('packet','source','seed','script','model'):
                assert shared['identity'][key] == identity[key]
            response = {key:value for key,value in shared['attempts'][0].items() if key not in ('program','report')}
            response['origin'] = 'shared_initial_proposal_from_reasoning'
        else:
            response = model_call(messages,think,len(task.names),seed)
        try:
            program = parse_response(response['raw']); report = task.check(program)
        except (ValueError,TypeError) as e:
            program = None; report = {'status':'invalid','error':str(e)}
        attempts.append({**response,'program':program,'report':report})
        result = {'identity':identity,'attempts':attempts,'status':report['status']}
        save_json(path,result)
        print(json.dumps({'task':task.id,'strategy':strategy,'attempt':index+1,'report':report,'seconds':response['seconds']}),flush=True)
        # Self review always gets its review opportunity, without oracle results.
        if report['status']=='verified' and strategy != 'self_review':
            break
        messages.append({'role':'assistant','content':response['raw']})
        if strategy == 'counterexample':
            messages.append({'role':'user','content':'Exact original group-word checker result: '+json.dumps(report)+'. Repair the COMPLETE equation system. Output JSON only.'})
        elif strategy == 'self_review':
            messages.append({'role':'user','content':'Independently review your preceding derivation for missing conditions, inverse and conjugation mistakes, and unnecessary restrictions. Re-derive from the definitions and return the complete corrected JSON.'})
        else:
            break
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--strategy',choices=['direct','reasoning','self_review','counterexample','audit'],default='audit')
    parser.add_argument('--tasks',default='reflection-0,reflection-1,rotation-0,dihedral-0')
    parser.add_argument('--out',default='results/research-v3')
    parser.add_argument('--seed',type=int,default=20261024)
    args = parser.parse_args()
    group = Transport()
    if args.strategy == 'audit':
        audit = group.audit()
        classifications = {}
        rng = random.Random(42)
        for name in args.tasks.split(','):
            kind,s = name.rsplit('-',1); task = Task(group,kind,int(s))
            for _ in range(200):
                v = rng.randrange(1 << task.n)
                assert task.residual(v) == [parity(r & v) for r in task.rows]
            classifications[name] = task.classification()
        report = {'audit':audit,'classifications':classifications,'status':'exact finite classification in fixed marked category',
                  'physical_completion':False}
        save_json(Path(args.out)/'classification.json',report); print(json.dumps(report,indent=2))
    else:
        for name in args.tasks.split(','):
            kind,s = name.rsplit('-',1)
            run(Task(group,kind,int(s)),args.strategy,args.out,args.seed)


if __name__ == '__main__':
    main()
