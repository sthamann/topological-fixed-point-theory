"""Independent simplicial, graph and history models with comparable readouts."""
from collections import Counter, defaultdict
from itertools import combinations
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def complex_data(name, n):
    if name == "interval":
        return [(i,i+1) for i in range(n)], []
    if name == "circle":
        return [(i,(i+1)%n) for i in range(n)], []
    if name == "figure8":
        return ([(i,(i+1)%n) for i in range(n)] +
                [(0,n), *((i,i+1) for i in range(n,2*n-2)), (2*n-2,0)]), []
    if name == "disk":
        return [], [(n,i,(i+1)%n) for i in range(n)]
    if name == "sphere":
        return [], [(pole,i,(i+1)%n) for pole in (n,n+1) for i in range(n)]
    if name == "projective_plane":
        # Six vertex minimal triangulation with two incident faces per edge.
        return [], [(0,1,2),(0,1,3),(0,2,4),(0,3,5),(0,4,5),
                    (1,2,5),(1,3,4),(1,4,5),(2,3,4),(2,3,5)]
    if name in {"cylinder","moebius","torus","klein"}:
        closed = name in {"torus","klein"}
        height = n if closed else 2
        faces = []

        def v(i,j):
            if i == n:
                i = 0
                if name == "moebius":
                    j = height-1-j
                if name == "klein":
                    j = -j
            return i*height+j%height

        for i in range(n):
            for j in range(height if closed else height-1):
                a,b,c,d = v(i,j),v(i+1,j),v(i+1,j+1),v(i,j+1)
                faces.extend(((a,b,c),(a,c,d)))
        return [], faces
    raise ValueError(f"unknown topology {name}")


TOPOLOGIES = ("interval","circle","figure8","disk","sphere","cylinder",
              "moebius","torus","klein","projective_plane")
DYNAMICS = ("local_generator","glued_boundary","oriented_local_system","hidden_history","face_gluing")


def gf2_rank(matrix):
    basis = {}
    for row in matrix.tolist():
        value = sum((int(x)%2)<<i for i,x in enumerate(row))
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def boundaries(name, n):
    explicit, raw_faces = complex_data(name,n)
    faces = sorted({tuple(sorted(face)) for face in raw_faces})
    edges = sorted({tuple(sorted(edge)) for edge in explicit} |
                   {edge for face in faces for edge in combinations(face,2)})
    vertices = sorted({v for edge in edges for v in edge})
    lookup = {v:i for i,v in enumerate(vertices)}
    d1 = sp.zeros(len(vertices),len(edges))
    d2 = sp.zeros(len(edges),len(faces))
    edge_index = {e:i for i,e in enumerate(edges)}
    for column,(a,b) in enumerate(edges):
        d1[lookup[a],column],d1[lookup[b],column] = -1,1
    for column, face in enumerate(faces):
        for i in range(3):
            edge = face[:i]+face[i+1:]
            d2[edge_index[edge],column] = (-1)**i
    if d1*d2 != sp.zeros(len(vertices),len(faces)):
        raise ArithmeticError("invalid boundary composition")
    return d1,d2,edges,faces


def rational_rank(matrix):
    return DomainMatrix.from_Matrix(matrix).convert_to(sp.QQ).rank()


def signature(d1,d2,edges,faces):
    sizes = (d1.rows,d1.cols,d2.cols)

    def betti(a,b):
        return [sizes[0]-a,sizes[1]-a-b,sizes[2]-b]

    incidence = Counter(edge for f in faces for edge in combinations(f,2))
    rational = betti(rational_rank(d1),rational_rank(d2))
    mod2 = betti(gf2_rank(d1),gf2_rank(d2))
    edge_index = {edge:i for i,edge in enumerate(edges)}
    face_dual = defaultdict(list)
    for edge,count in incidence.items():
        if count != 2:
            continue
        row = d2[edge_index[edge],:]
        i,j = [k for k in range(row.cols) if row[k]]
        sign = -int(row[i]*row[j])
        face_dual[i].append((j,sign))
        face_dual[j].append((i,sign))
    orientable = None
    if faces and max(incidence.values(),default=0) <= 2:
        orientable = True
        assigned = {}
        for start in range(len(faces)):
            if start in assigned:
                continue
            pending = [(start,1)]
            while pending:
                i,sign = pending.pop()
                if i in assigned:
                    orientable &= assigned[i] == sign
                    continue
                assigned[i] = sign
                pending.extend((j,sign*twist) for j,twist in face_dual[i])
    # Boundary components are connected components of the boundary edge graph.
    boundary_neighbors = defaultdict(set)
    for (a,b),count in incidence.items():
        if count == 1:
            boundary_neighbors[a].add(b)
            boundary_neighbors[b].add(a)
    left = set(boundary_neighbors)
    boundary_count = 0
    while left:
        boundary_count += 1
        stack = [min(left)]
        while stack:
            v = stack.pop()
            if v in left:
                left.remove(v)
                stack.extend(boundary_neighbors[v] & left)
    return dict(cells=sizes,euler=sizes[0]-sizes[1]+sizes[2],
                betti_Q=rational,betti_F2=mod2,orientable=orientable,
                boundary_components=boundary_count if faces else None)


def mod2_cocycle(d1,d2):
    # Row-reduce the face constraints over F2, then seek a vector outside
    # the coboundary space. A valid sign twist has trivial face holonomy.
    rows = [[int(x)%2 for x in row] for row in d2.T.tolist()]
    pivots = []
    for j in range(d1.cols):
        k = next((i for i in range(len(pivots),len(rows)) if rows[i][j]),None)
        if k is None:
            continue
        r = len(pivots)
        rows[r],rows[k] = rows[k],rows[r]
        for i,row in enumerate(rows):
            if i != r and row[j]:
                rows[i] = [a^b for a,b in zip(row,rows[r])]
        pivots.append(j)
    rank = gf2_rank(d1)
    for free in sorted(set(range(d1.cols))-set(pivots)):
        cochain = [int(j == free) for j in range(d1.cols)]
        for i,pivot in enumerate(pivots):
            cochain[pivot] = rows[i][free]
        if gf2_rank(d1.col_join(sp.Matrix([cochain]))) > rank:
            return cochain
    return [0]*d1.cols


def response(matrix, observe):
    initial = sp.eye(matrix.rows)[:,0]
    vector = initial
    moments = []
    for _ in range(4):
        moments.append((initial.T*vector)[0])
        vector = matrix*vector
    n = min(observe,matrix.rows)
    q = matrix[:n,:n]
    q2 = (matrix**2)[:n,:n]
    defect = q2-q*q
    return {"return_moments":list(map(str,moments)), "trace_defect":str(defect.trace()),
            "omitted_path_rank":rational_rank(defect),
            "operator_dimension":matrix.rows, "fixed_modes":matrix.rows-rational_rank(matrix-sp.eye(matrix.rows))}


def run(topology,n,dynamics,metric="unit"):
    d1,d2,edges,faces = boundaries(topology,n)
    invariants = signature(d1,d2,edges,faces)
    # Conductances are hypotheses, explicitly selected rather than fitted.
    weights = [sp.Integer(1) if metric == "unit" else sp.Rational(1,2) if i%2 else sp.Integer(1)
               for i in range(len(edges))]
    w = sp.diag(*weights)
    observe = d1.rows
    twist = [0]*d1.cols
    if dynamics == "local_generator":
        # First-order positive transfer supplied on vertices before any readout.
        laplacian = d1*w*d1.T
        matrix = sp.eye(laplacian.rows)-laplacian/(1+max(sum(abs(x) for x in row) for row in laplacian.tolist()))
    elif dynamics in {"glued_boundary","oriented_local_system"}:
        # No primitive clock. Eliminate edge amplitudes from a positive block form.
        if dynamics == "oriented_local_system":
            twist = mod2_cocycle(d1,d2)
        incidence = sp.Matrix(d1)
        for i,value in enumerate(twist):
            if value:
                for j in range(incidence.rows):
                    if incidence[j,i] == 1:
                        incidence[j,i] = -1
        matrix = (sp.eye(incidence.rows)+incidence*w*incidence.T).inv(method="DM")
    elif dynamics == "face_gluing":
        # Edge histories couple through vertices and filled triangles.
        # The two incidence operators share the same cell complex.
        laplacian = d1.T*d1+d2*d2.T
        metric_form = sp.diag(*(1/value for value in weights))
        # Return a self-adjoint edge transfer in the unit metric chart.
        # The conductance perturbation is additive to the Hodge operator.
        matrix = (laplacian+metric_form).inv(method="DM")
        observe = min(d1.rows, d1.cols)
    elif dynamics == "hidden_history":
        # The readout keeps vertices only. The transition acts on vertices and edges.
        boundary = d1.applyfunc(abs)*w
        adjacency = sp.zeros(d1.rows).row_join(boundary).col_join(
            boundary.T.row_join(sp.zeros(d1.cols)))
        degree = sp.diag(*(sum(row) for row in adjacency.tolist()))
        laplacian = degree-adjacency
        matrix = (sp.eye(laplacian.rows)+laplacian).inv(method="DM")
    else:
        raise ValueError("unknown dynamics scenario")
    return dict(status="EXACT_SCOPED",topology=topology,resolution=n,dynamics=dynamics,
                metric=metric,invariants=invariants,twist_support=sum(twist),
                response=response(matrix,observe), physical_gate_passes=[],
                scope="finite simplicial topology and exact rational transfers with declared conductances")


def discriminate(records):
    # Compare topology-independent scalar readouts before exposing exact topology labels.
    keys = ("return_moments","trace_defect","omitted_path_rank","fixed_modes")
    pairs = list(combinations(records,2))
    return sorted((dict(observable=k,pairs_separated=sum(a["response"][k]!=b["response"][k] for a,b in pairs),
                        comparable_pairs=len(pairs)) for k in keys),
                  key=lambda item:(-item["pairs_separated"],item["observable"]))
