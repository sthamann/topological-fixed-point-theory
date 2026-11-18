"""Exact E8 winding and circle spectral flow for the old marked seam legs."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product

D=[]
for i,j in combinations(range(8),2):
    for s,t in product((-1,1),repeat=2):
        D.append(tuple(s*(k==i)+t*(k==j) for k in range(8)))
S=[tuple(F(x,2) for x in signs) for signs in product((-1,1),repeat=8) if sum(signs)%4==0]
p1=(F(1,2),)*8
p2=tuple(F(-1 if k in (0,5) else 1,2) for k in range(8))
v=(0,0,0,0,0,0,0,-1)
w=(0,2,4,-2,0,4,2,0)
dot=lambda x,y:sum(a*b for a,b in zip(x,y))
gamma=tuple(a-b for a,b in zip(p1,p2))
twov=tuple(2*x for x in v)
assert len(D)==112 and len(S)==128 and p1 in S and p2 in S and gamma in D
assert all((dot(w,p)-dot(twov,p))%2==0 for p in D+S)
assert all(dot(twov,p)%2==(p in S) for p in D+S)
assert all(tuple(x+y for x,y in zip(p,v))[-1]==-F(1,2) for p in (p1,p2))
assert all(sum(x*x for x in tuple(y+z for y,z in zip(p,v)))==2 for p in (p1,p2))
assert all(sum(x*x for x in tuple(y+z for y,z in zip(p,twov)))==4 for p in (p1,p2))
assert dot(twov,gamma)==0
def sf(m):
    # Exact crossings of n+1/2+s*m, 0<s<1, of the complete circle spectrum.
    sign=(m>0)-(m<0)
    cross=tuple(-F(2*n+1,2*m) for n in range(-abs(m)-1,abs(m)+1)
                if m and 0<-F(2*n+1,2*m)<1)
    return len(cross)*sign, cross
root_flows=Counter(sf(int(dot(twov,p)))[0] for p in D+S)
assert root_flows=={-2:14,-1:64,0:84,1:64,2:14}
assert sf(-1)==(-1,(F(1,2),)) and sf(1)==(1,(F(1,2),))
print('Old g8^4 and exp(2pi i <v,p>) agree on all 240 E8 roots and hence on L.')
print('Marked original odd legs p1,p2 have the same flux -1 and one crossing each at s=1/2.')
print('Their common middle current has flux 0; reversed charges in the same angular orientation have flux +1.')
print('Root flux counts:',dict(sorted(root_flows.items())), 'sum',sum(m*n for m,n in root_flows.items()))
print('Actual v shift preserves both original ground weights h=1; 2v as a momentum shift gives h=2, so use 2v only as the closed torus cocharacter for the circle index.')
assert (dot(p1,gamma),dot(p2,gamma))==(1,-1)
assert tuple(x-y for x,y in zip(p1,gamma))==p2
assert tuple(x+y for x,y in zip(p2,gamma))==p1
assert tuple(x+y for x,y in zip(p1,gamma)) not in D+S
assert tuple(x-y for x,y in zip(p2,gamma)) not in D+S
print('The original middle root connects the two legs as one irreducible SU(2) doublet. If both are simultaneous fibres, the circle flux is -2 and no projector of rank one commutes with the full middle-root action.')
for split in (F(1,4),F(3,4)):
    cross=sf(-1)[1]
    parts=tuple(-sum(lo<x<hi for x in cross) for lo,hi in ((0,split),(split,1)))
    assert sum(parts)==-1
    print('Serially glued charge lines, split',split,': signed crossings',parts,'sum -1.')
sign_p2,sign_p1=((-1)**(int(dot(twov,p))%2) for p in (p2,p1))
assert sign_p2==sign_p1==-1
print('Original ordered point word -i/4: one inner phase +i/4; two inner phases -i/4.')
