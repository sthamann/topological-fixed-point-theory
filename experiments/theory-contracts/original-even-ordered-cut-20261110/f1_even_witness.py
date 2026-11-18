"""Solve for a first-local-mode input whose early order vanishes, if one exists."""
from contextlib import redirect_stdout
from pathlib import Path
import io,json,runpy,sys
with redirect_stdout(io.StringIO()):
    a=runpy.run_path(str(Path(__file__).with_name("f1_even_fused_words.py")))
sp,K,q,r,B,M,V,old,src=[a[x] for x in ("sp","K","q","r","B","M","V","old","src")]
fused=a["fused"]
def num(x):
    return K.from_sympy(sp.expand(x))
def express(states,vector):
    return src["combine"]((w,coef*num(u)) for state,u in zip(states,vector,strict=True) for w,coef in state.items())
d=q+r
weak=sp.eye(8)[:,3]-sp.eye(8)[:,4]
Y=sp.Matrix([-2,-2,-2,3,3,0,0,0])/6
assert all(x==0 for x in (M*d)[:3]) and (M*d).dot(Y)==sp.Rational(1,2) and (M*d).dot(weak)==-1
cases={}
for alpha in (sp.Rational(1,4),sp.Rational(3,4)):
    drift=sp.Matrix([j+alpha-4 for j in range(8)])
    co=[num(x) for x in old["cl"](V.T*drift/8)]
    prepared=[old["series_response"](f,co,4) for f in src["polys"]]
    inputs=[old["input_word_map"](v) for v in prepared]
    raw_q=[old["vertex_zero"](v) for v in prepared]
    Gq=src["gram"](raw_q)
    rows={}
    for mode in (-2,-1,0,1,2):
        early=[fused(v,mode,False) for v in inputs]
        late=[fused(v,mode,True) for v in inputs]
        A,B=src["gram"](early),src["gram"](late)
        ker=A.to_DM().convert_to(K).nullspace().to_Matrix()
        stack=src["rank"](src["gram"](early+late))
        survivors=(ker.conjugate()*B*ker.T).applyfunc(sp.expand)
        item={"early_rank":src["rank"](A),"late_rank":src["rank"](B),"all_endpoints_rank":stack,
              "rank_of_late_on_exact_early_nullspace":src["rank"](survivors)}
        if mode in (-1,2):
            if ker.rows:
                for row in ker.tolist():
                    c=sp.Matrix(row)
                    y=express(early,c)
                    z=express(late,c)
                    assert src["inner"](y,y)==0
                    n=src["inner"](z,z)
                    if n:
                        nr=sp.expand((c.H*Gq*c)[0])
                        item["one_exact_early_null_input"]=list(map(str,row))
                        item["late_raw_norm"]=str(n)
                        item["late_per_q_input_norm"]=str(sp.simplify(n/nr))
                        break
        if mode==2:
            assert not any(early) and all(set(v)<={()} for v in late)
            coeffs=[K.to_sympy(v.get((),K.zero)) for v in late]
            item["bare_integer_vertex_mode"]=-2
            item["normone_end_coefficients"]=list(map(str,coeffs))
            item["final_energy"]=str(16+drift.dot(d))
        rows[str(mode)]=item
    cases[str(alpha)]=rows
print(json.dumps({"end_charge":list(map(str,d)),"cases":cases},indent=2))
