"""Interval exclusion for all 216 integer splits of the lepton word lengths."""
from itertools import combinations, product
from mpmath import iv
import math

iv.dps=35
I=iv.mpf
one=I(1)
phi=1/(6*iv.pi)+3/(256*iv.pi**4)
lam=iv.sqrt(phi*(1-phi))
ud=I(55)/117*(I(34)/47)*(I(3)/26)*(1-phi)**2/phi
tab_e=tuple(map(I,("0.466","1.085","0.899")))
tab_d=tuple(map(I,("0.936","0.943","0.480")))
half_unit=I("0.0005")
ed_lo=math.prod(a-half_unit for a in tab_e)/math.prod(a+half_unit for a in tab_d)
ed_hi=math.prod(a+half_unit for a in tab_e)/math.prod(a-half_unit for a in tab_d)
x_min=9*(ed_lo-1)/(9-ed_lo)
x_max=9*(ed_hi-1)/(9-ed_hi)
# Enlarge the true outward enclosure before a binary splitting of endpoints.
x_bounds=(math.nextafter(float(x_min.a),-math.inf),math.nextafter(float(x_max.b),math.inf))
x_whole=I(list(x_bounds))
lengths=(8,5,3)
pairs=tuple(combinations(range(3),2))
r_e=12*phi**2/7
r_mu=8*phi/7
masses=(r_e*r_mu,r_mu,one)
mass_det=math.prod(masses)
F_goal=(sum(y*y for y in masses)**3)/(mass_det**2)
C_goal=(sum((masses[i]*masses[j])**2 for i,j in pairs)**3)/(mass_det**4)
D_scale=2*lam**sum(lengths)/243
F_diag=(4*lam**16+4*lam**10+64*lam**6)/9
powers=[lam**(2*n) for n in range(2*sum(lengths)+1)]

def interval_for_mixing(p,x):
    A=144+81*p
    C=16*p
    k=ud*A*(1+x/9)/144
    coefficient=256*ud*(1+x/9)/144-16*x/81
    m=144*coefficient+p*(81*coefficient-16)
    assert m.a>0
    w=(m+iv.sqrt(m*m+1024*k*C))/512
    t2=x*A/(256*w+C)
    return A,64*w*t2/81,p*(1+16*t2/81)

def weights(left):
    right=tuple(lengths[i]-left[i] for i in range(3))
    X=powers[left[0]+right[1]]+powers[left[1]+right[0]]
    Y=powers[left[1]+right[2]]+powers[left[2]+right[1]]
    W=tuple(powers[left[i]+left[j]+right[k]+right[l]] for i,j in pairs for k,l in pairs)
    return X,Y,W

def bounds(p,x,X,Y,W):
    A,u,v=interval_for_mixing(p,x)
    F=F_diag+X*u+Y*v
    D=D_scale*A*(1+x)
    f=F**3/(D**2*F_goal)
    if not (f.a<=1<=f.b):
        return f,None
    terms=((I(4)/9+u)**2,4*v/9,u*v,4*v/9,I(256)/81,64*u/9,u*v,64*u/9,(I(16)/9+v)**2)
    C=sum(a*b for a,b in zip(W,terms))
    return f,C**3/(D**4*C_goal)

max_depth=0
interval_count=0
closest_margin=math.inf
for left in product(*(range(n+1) for n in lengths)):
    X,Y,W=weights(left)
    # For p>=1, D<=225*D_scale*(1+x_max)*p and F>=Y*p,
    # so p strictly above this cutoff cannot fit the Frobenius invariant.
    cutoff=F_goal*(225*D_scale*(1+x_whole))**2/Y**3
    p_upper=max(2,math.ceil(float(cutoff.b))+2)
    stack=[(0.,1.,*x_bounds,0)]
    lower=1.
    while lower<p_upper:
        high=min(2*lower,p_upper)
        stack.append((lower,high,*x_bounds,0))
        lower=high
    while stack:
        lo,hi,xlo,xhi,depth=stack.pop()
        interval_count+=1
        max_depth=max(max_depth,depth)
        f,c=bounds(I([lo,hi]),I([xlo,xhi]),X,Y,W)
        if c is None:
            continue
        if not (c.a<=1<=c.b):
            closest_margin=min(closest_margin,min(abs(float(c.a)-1),abs(float(c.b)-1)))
            continue
        if depth==45:
            raise RuntimeError((left,lo,hi,xlo,xhi,str(f),str(c)))
        if (hi-lo)/(1+lo) > (xhi-xlo)/x_bounds[0]:
            mid=(lo+hi)/2
            stack.extend(((lo,mid,xlo,xhi,depth+1),(mid,hi,xlo,xhi,depth+1)))
        else:
            mid=(xlo+xhi)/2
            stack.extend(((lo,hi,xlo,mid,depth+1),(lo,hi,mid,xhi,depth+1)))

assert closest_margin>0
print("All 216 integer word-length splits excluded by outward-rounded interval inequalities.")
print("Intervals:",interval_count,"max depth:",max_depth,"closest recorded exterior-power invariant margin:",closest_margin)
