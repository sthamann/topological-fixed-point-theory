"""Character counts on unoriented and oriented opposite-flag edges."""
import sys
from pathlib import Path
from itertools import combinations
from collections import Counter
import sympy as s
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"fortsetzung"))
from apartment import FLAGS,linear_maps,LINES

G=list(linear_maps())
def act(g,f):
    p,line=FLAGS[f]
    return FLAGS.index((g[p],tuple(sorted(g[x] for x in line))))
edges=[(i,j) for i,j in combinations(range(21),2) if FLAGS[i][0] not in FLAGS[j][1]
       and FLAGS[j][0] not in FLAGS[i][1]]
inc=s.zeros(21,len(edges))
for k,(i,j) in enumerate(edges):
    inc[i,k]=-1
    inc[j,k]=1

def sign(p):
    return (-1)**sum(x>y for x,y in combinations(p,2))

def char(g):
    perm=[act(g,i) for i in range(21)]
    vertices=sum(x==y for x,y in enumerate(perm))
    edge=sum(int((perm[i],perm[j])==(i,j))-int((perm[j],perm[i])==(i,j)) for i,j in edges)
    return 1-vertices+edge

def main():
    assert inc.rank()==20 and len(edges)==84
    chars={g:char(g) for g in G}
    # Stabilizer of u induces S4 on lines not through u, equivalently
    # on the four planes f=1 of a fixed-point affine representation dual.
    u=7
    lines=[line for line in LINES if u not in line]
    def sign_at_u(g):
        return sign([lines.index(tuple(sorted(g[x] for x in line))) for line in lines])
    groups=[
        ("whole incidence",G,lambda g:1),
        ("one mark",[g for g in G if g[u]==u],lambda g:1),
        ("one mark, sign",[g for g in G if g[u]==u],sign_at_u),
        ("unordered apartment",[g for g in G if {g[x] for x in (1,2,4)}=={1,2,4}],lambda g:1),
        ("apartment sign",[g for g in G if {g[x] for x in (1,2,4)}=={1,2,4}],lambda g:sign([g[x] for x in (1,2,4)])),
        ("one flag",[g for g in G if act(g,0)==0],lambda g:1)]
    for name,H,chi in groups:
        mult=s.Rational(sum(chi(g)*chars[g] for g in H),len(H))
        print(name,len(H),"cycle character multiplicity",mult)
    print("cycle rank",len(edges)-inc.rank(),"; character values with frequencies",Counter(chars.values()))
    assert sum(chars.values())==0

if __name__=="__main__":
    main()
