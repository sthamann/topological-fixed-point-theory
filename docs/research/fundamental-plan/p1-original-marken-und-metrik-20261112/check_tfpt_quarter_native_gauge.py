"""Compare the 32 original quarter lifts with the fixed SM gauge axes."""
from contextlib import redirect_stdout
from collections import Counter
from itertools import combinations,product
from pathlib import Path
import io
import sys

repo=Path(sys.argv[1])
src=repo/"docs/research/fundamental-plan/check_tfpt_charge_clock_20261112.py"
# Reuse the native lifts and their actual signed vector-weight D8 frame.
prefix=src.read_text().split("\nws=[]",1)[0]
ns={}
with redirect_stdout(io.StringIO()):
    exec(compile(prefix,str(src),"exec"),ns)
quarter_lifts=ns["g_perm"]
assert len(quarter_lifts)==32
cycles=Counter()
def cycle_lengths(permutation):
    seen=set()
    lengths=[]
    for i in range(len(permutation)):
        if i in seen:
            continue
        x=i
        length=0
        while x not in seen:
            seen.add(x)
            length+=1
            x=permutation[x]
        lengths.append(length)
    return tuple(sorted(lengths))

all_y=[]
for weak in combinations(range(8),2):
    other=[j for j in range(8) if j not in weak]
    for colour in combinations(other,3):
        for signs in product((-1,1),repeat=5):
            values=[0]*8
            for j,s in zip(weak+colour,signs):
                values[j]=s*(3 if j in weak else 2)
            all_y.append(tuple(values))
assert len(set(all_y))==17920
outcomes=Counter()
for g in sorted(quarter_lifts):
    cycles[cycle_lengths(tuple(i for i,_ in g))]+=1
    # A fixed hypercharge U(1) permits only a sign on its generator.
    count=0
    for y in all_y:
        gy=tuple(s*y[i] for i,s in g)
        count+=gy==y or gy==tuple(-x for x in y)
    outcomes[count]+=1
print('Unsigned coordinate cycle types of all compatible G in the D8 chart:',dict(cycles))
print('Count of possible fixed hypercharge subgroups per G:',dict(outcomes))
assert cycles=={(4,4):32} and outcomes=={0:32}
assert all(not any(ns["apply"](g,r)==r for r in ns["roots"]) for g in ns["compatible"])
# G-fixed Cartan vectors in this signed frame have constant absolute
# coordinate values along each 4-cycle. The orbit average lies there.
for g in quarter_lifts:
    def rotate(y):
        return tuple(sign*y[index] for index,sign in g)
    y=all_y[0]
    yy=[y]
    for _ in range(3):
        yy.append(rotate(yy[-1]))
    avg=tuple(sum(z) for z in zip(*yy))
    assert rotate(avg)==avg
print('All 240 roots moved by every compatible G. A G-fixed Cartan vector has support 0, 4, or 8; no orbit mean gives a hypercharge profile with 5 nonzero axes.')
