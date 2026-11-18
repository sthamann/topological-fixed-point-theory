from pathlib import Path
import sys, io, runpy
from contextlib import redirect_stdout
repo=Path(sys.argv[1])
src=repo/"experiments/theory-contracts/original-paper-source-reconciliation-20261102/contribution_20261106Q/source"
sys.path[:0]=[str(repo),str(src)]
from native_cx_zero_mode import cx_zero_mode
prior=repo/"experiments/theory-contracts/recursive-mirrored-source-20261103/check_seam_ordered_yukawa_20261111.py"
with redirect_stdout(io.StringIO()):
 d=runpy.run_path(str(prior))
def actions(ii,f):
 ops=[cx_zero_mode(i,f) for i in ii]
 return ops[0]["source"], tuple({q:sum(q in row[key] for row in ops) for q in range(240)} for key in ("action","adjoint_action"))
def legs(g,j,v):
 return [v[g.ridx[tuple(2*x for x in d["matter_vertex"](j,k)[0])]] for k in range(3)]

sp=d["sp"]
groups=((0,1,2),(3,4))
unit=lambda i:[int(i==j) for j in range(3)]
color={"e":(0,3),"u":(2,1),"d":(2,1)}
weak={"e":(1,0),"u":(1,2),"d":(1,0)}
for count_spec,indices in zip((color,weak),groups):
 for f in range(4):
    g,(forward,back)=actions(indices,f)
    for tag,rows in d["slots"].items():
        for p,q,_ in rows:
            assert legs(g,p,forward)==legs(g,q,forward)==[0,0,0]
            expected=([0,0,0],[0,0,0]) if f==0 else tuple(
                [m*z for z in unit(f-1)] for m in count_spec[tag])
            assert (legs(g,p,back),legs(g,q,back))==expected
    # Gauge roots are in the original carrier. These exact sums commute
    # with the grade-one SU3 and SU2 transitions, including the Cartan outputs.
    for i,j in ((0,1),(1,2),(3,4)):
        root=tuple(2*(int(k==i)-int(k==j)) for k in range(8))
        nroot=tuple(-z for z in root)
        for mark in (forward,back):
            assert mark[g.ridx[root]]==mark[g.ridx[nroot]]==0
            for idx,p in enumerate(g.roots):
                q=tuple(x+y for x,y in zip(p,root))
                if q in g.ridx:
                    assert mark[idx]==mark[g.ridx[q]]
print("Original CX color and weak return counts, e L/R, u L/R, d L/R:",color,weak)
print("The summed positive return forms commute with the actual carrier root generators.")

# Add any positive diagonal baseline shared by e and d on the same leg,
# and allow different nonnegative color and weak weights on each leg.
Gl,Gr,cl,cr,wl=sp.symbols("Gl Gr cl cr wl",positive=True)
le,ld=Gl+wl,Gl+wl+2*cl
re,rd=Gr+3*cr,Gr+cr
assert sp.expand(ld-le)==2*cl
assert sp.expand(3*rd-re)==2*Gr
# L_d/L_e >=1 and R_d/R_e >=1/3 imply a lower bound on the
# electron/down lightest singular-value ratio in the nested quarter chart.
lam=d["lY"]
mp=d["mp"]
table_max=lam*(mp.mpf(d["TAB"]["e"][1])+d["half_printing"])/(
    mp.mpf(d["TAB"]["d"][1])-d["half_printing"])
bound=lam/mp.sqrt(3)
assert table_max<bound
print("Even any positive linear CX-end norm gives e/d >=lambda/sqrt(3) =",
      mp.nstr(bound,12), "above the old rounded table ceiling",mp.nstr(table_max,12))

# Same positive family-dependent distribution of a dimensionless color heat
# time: joint propagation sees n_L+n_R=3 for both e,d; separate norms do not.
assert sum(color["e"])==sum(color["d"])==3
x,y=sp.symbols("x y",nonnegative=True)
# For independent X,Y in [0,1] with the same distribution, this identity
# gives 2(f3-f1 f2) as E[(X-Y)^2(X+Y)] >= 0.
assert sp.expand((x-y)**2*(x+y)-(x**3+y**3-x*y**2-y*x**2))==0
# An unsuppressed atom of weight p forces f1>=p and f2>=f3.
# Thus even separately averaged first-family light blocks have e/d >=lambda sqrt(p).
ratio_bound=table_max/lam
assert ratio_bound<mp.mpf(1)/2
minimum_time=-2*mp.log(ratio_bound)
assert minimum_time>mp.mpf("1.39")
print("Common color-time averaging is trivial on the joint e/d return (equal total count 3).")
print("Separate diagonal norms from positive color-time mixtures have residual |det E/det D| <=1.")
print("If one of four equal-weight times is unsuppressed, the first-family factor cannot go below 1/2.")
print("With bounded times, the light-ratio test alone requires a maximum dimensionless heat time >",
      mp.nstr(minimum_time,12))
