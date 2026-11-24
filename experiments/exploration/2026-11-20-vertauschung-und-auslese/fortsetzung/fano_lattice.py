"""Finite exact audits of Fano panels, parity completion, and the E8 lattice."""
from collections import Counter
from itertools import combinations, permutations, product
from math import comb, factorial
from fractions import Fraction


def support(indices):
    return sum(1 << i for i in indices)


def span(rows):
    code = {0}
    for r in rows:
        code |= {x ^ r for x in code}
    return sorted(code)


def panels():
    # Coordinates 1..7 represent nonzero points of F_2^3; bit zero is reserved
    # for the extra parity coordinate. The ith nonzero point occupies bit i.
    return sorted({support((a,b,a^b)) for a,b in combinations(range(1,8),2)})


def compatible_families(rows, selected=()):
    if len(selected)==7:
        yield selected
        return
    for i, row in enumerate(rows):
        if len(rows)-i < 7-len(selected):
            break
        rest = [r for r in rows[i+1:] if (r & row).bit_count()==1]
        yield from compatible_families(rest,selected+(row,))


def even_extension(code):
    return sorted(x | (x.bit_count()%2) for x in code)


def construction_roots(code):
    # z represents the normalized real lattice vector z/sqrt(2).
    # At norm two, sum z_i^2=4 leaves just one entry +/-2 or four entries +/-1.
    roots=set()
    for i,sgn in product(range(8),(-1,1)):
        z=[0]*8
        z[i]=2*sgn
        roots.add(tuple(z))
    for c in code:
        if c.bit_count()!=4:
            continue
        block=[i for i in range(8) if c>>i&1]
        for signs in product((-1,1),repeat=4):
            z=[0]*8
            for i,s in zip(block,signs):
                z[i]=s
            roots.add(tuple(z))
    return sorted(roots)


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def simple_system(roots):
    # Powers of five determine a regular vector on all the integer roots.
    order=lambda r:sum(x*5**i for i,x in enumerate(r))
    positive={r for r in roots if order(r)>0}
    assert len(positive)*2==len(roots)
    decomposable={tuple(x+y for x,y in zip(a,b)) for a in positive for b in positive}
    simple=sorted(positive-decomposable)
    n=len(simple)
    gram=[[Fraction(dot(a,b),2) for b in simple] for a in simple]
    assert all(gram[i][j] in (0,-1) for i in range(n) for j in range(i))
    graph=[{j for j in range(n) if gram[i][j]==-1} for i in range(n)]
    center,=[i for i,g in enumerate(graph) if len(g)==3]
    arms=[]
    for first in graph[center]:
        previous, current, length=center, first, 1
        while len(graph[current])==2:
            following,=graph[current]-{previous}
            previous,current=current,following
            length+=1
        arms.append(length)
    return simple, sorted(arms)


def primitive_and_symmetry_tests(roots, simple, selected, perp):
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    basis = Matrix(simple)
    def elementary_divisors(subroots):
        c = Matrix(subroots)*basis.inv()
        assert all(x.is_Integer for x in c)
        return tuple(abs(int(x)) for x in smith_normal_form(c,domain=ZZ).diagonal() if x)
    a_div=elementary_divisors(selected)
    d_div=elementary_divisors(perp)
    joint_div=elementary_divisors(selected+perp)
    assert a_div==(1,1,1) and d_div==(1,1,1,1,1) and joint_div==(1,1,1,1,1,1,1,4)
    print("Smith divisors for A3, complement and direct sum in E8:",a_div,d_div,joint_div)
    coords = [[int(v) for v in row] for row in (Matrix(selected+perp)*basis.inv()).tolist()]
    def mod4_pair(row,character):
        return sum(x*y for x,y in zip(row,character))%4
    primitive_character=next(
        chi for chi in product(range(4),repeat=8)
        if any(x%2 for x in chi) and all(mod4_pair(row,chi)==0 for row in coords))
    root_coords = (Matrix(roots)*basis.inv()).tolist()
    cosets = Counter(mod4_pair([int(v) for v in row],primitive_character) for row in root_coords)
    assert [cosets[k] for k in range(4)]==[52,64,60,64]
    print("Z4 glue class 0,1,2,3 root counts:",[cosets[k] for k in range(4)])


    # Count ALL embedded A3 root subsystems, and all invariant ones for this S3.
    # For any embedded A3, the chosen chamber determines a unique simple triple
    # (up to swapping endpoints). The six positive roots complete that triple.
    def positive_root(z):
        return z if tuple(reversed(z))>tuple([0]*8) else tuple(-v for v in z)
    positive=sorted({positive_root(z) for z in roots})
    rid={z:i for i,z in enumerate(positive)}
    def add(*vs):
        return tuple(map(sum,zip(*vs)))
    systems=set()
    for b in positive:
        adjacent=[a for a in positive if dot(a,b)==-2]
        for a,c in combinations(adjacent,2):
            if dot(a,c)==0:
                systems.add(frozenset(rid[z] for z in (a,b,c,add(a,b),add(b,c),add(a,b,c))))

    # GL_2(F2) on the first two bits, preserving the direct-sum complement
    # F2*4. It has orbits on coordinates {0},{4},{1,2,3},{5,6,7}.
    generators=[]
    for p in ([0,2,1,3],[0,2,3,1]):
        generators.append(tuple((k&4)+p[k%4] for k in range(8)))
    def line_action(perm,signs):
        return [rid[positive_root(tuple(s*z[k] for s,k in zip(signs,perm)))] for z in positive]
    actions = [line_action(perm,[1]*8) for perm in generators]
    singleton_orbits = ((0,),(4,),(1,2,3),(5,6,7))
    sign_actions = [line_action(range(8),[-1 if k in block else 1 for k in range(8)])
                    for block in singleton_orbits]
    def fixed(family,transforms):
        return [sys for sys in family if all(frozenset(tr[i] for i in sys)==sys for tr in transforms)]
    fixed_s3=fixed(systems,actions)
    fixed_signed=fixed(fixed_s3,sign_actions)
    gl3 = []
    for u,v,w in permutations(range(1,8),3):
        if u^v==w:
            continue
        p = [0,u,v,u^v,w,u^w,v^w,u^v^w]
        gl3.append(line_action(p,[1]*8))
    assert len(gl3)==168
    fixed_gl3=fixed(fixed_s3,gl3)
    print("automorphism group of Fano incidence alone, invariant A3 count:",len(gl3),len(fixed_gl3))
    assert (len(systems),len(fixed_s3),len(fixed_signed))==(7560,18,0)
    print("all embedded A3, S3-invariant, and invariant including S3-commuting signs:",
          len(systems),len(fixed_s3),len(fixed_signed))


def main():
    lines=panels()
    raw=span(lines)
    ext=even_extension(raw)
    assert Counter(c.bit_count() for c in raw)=={0:1,3:7,4:7,7:1}
    assert Counter(c.bit_count() for c in ext)=={0:1,4:14,8:1}
    assert all((u&v).bit_count()%2==0 for u,v in product(ext,repeat=2))
    basis=[]
    for w in raw:
        if w not in span(basis):
            basis.append(w)
    assert len(basis)==4
    linear_choices=[]
    for vals in product(range(2),repeat=len(basis)):
        extension=span([b|v for b,v in zip(basis,vals)])
        linear_choices.append(all(x.bit_count()%4==0 for x in extension))
    assert sum(linear_choices)==1

    candidates=[support(triple) for triple in combinations(range(1,8),3)]
    designs=list(compatible_families(candidates))
    assert len(designs)==30 and all(len(span(d))==16 for d in designs)
    print("punctured and extended weights:",sorted(Counter(x.bit_count() for x in raw).items()),
          sorted(Counter(x.bit_count() for x in ext).items()))
    print("doubly even linear one-bit extensions:",sum(linear_choices),"/",len(linear_choices))
    print("panel controls with seven distinct 3-subsets on seven points:",len(designs),"/",comb(35,7))

    roots=construction_roots(ext)
    assert len(roots)==240
    simple,arms=simple_system(roots)
    assert len(simple)==8 and arms==[1,2,4]
    # The unscaled preimage of this four-dimensional binary code has index 16
    # in Z^8. Rescaling each coordinate by 1/sqrt(2) has volume factor 2^-4.
    print("Construction-A covolume, roots, simple rank, E8 arm lengths:",1,len(roots),len(simple),arms)

    # A conjugation-invariant S3 complement in GL3(F2), fixed point 4 and
    # preserved projective line {1,2,3}. No coordinate triple supports a
    # weight-four codeword, so any 3-coordinate section realizes A1^3, not A3.
    sections=Counter(sum(all(z[j]==0 for j in range(8) if j not in ids) for z in roots)
                     for ids in combinations(range(8),3))
    assert sections=={6:56}
    print("3-coordinate root counts over all 56 sections:",dict(sections))

    # Intrinsic coordinate signs of the Construction-A lattice fix all binary
    # codewords. An invariant real subspace for *every* coordinate sign is a
    # coordinate subspace (diagonal characters are pairwise distinct).
    # Thus Fano incidence + the lattice alone cannot pick a canonical A3
    # invariant under these symmetries. The following explicit extra pair tests
    # what the orthogonal root subsystem would be after choosing an A3.
    tetrad = sorted(i for i in range(8) if (ext[1]>>i)&1)
    assert len(tetrad)==4
    i,j = tetrad[:2]
    selected=[z for z in roots if z[i]==z[j] and all(z[k]==0 for k in range(8) if k not in tetrad)]
    perp=[z for z in roots if all(dot(z,a)==0 for a in selected)]
    assert len(selected)==12 and len(perp)==40
    print("after choosing a pair in a tetrad, roots in A3 and orthogonal D5:",len(selected),len(perp))
    primitive_and_symmetry_tests(roots,simple,selected,perp)


if __name__=="__main__":
    main()
