"""Test whether native one-window marginals determine the R032 joint source."""
from itertools import combinations
import json

from flint import arb, arb_mat, ctx

from r113 import ball, conditional_word, squared_singular_values, WINDOWS
from source import NATIVE, source_refs
from verify import digest


def independent_direction(rows, mixed):
    """Use an invertible minor to vary the mixed mark with all fitted forms fixed."""
    count = len(rows)
    columns = len(mixed)
    best = None
    for support in combinations(range(columns), count + 1):
        matrix = arb_mat([[row[i] for i in support] for row in (*rows, mixed)])
        separation = abs(matrix.det()).lower()
        if separation > 0 and (best is None or separation > best[0]):
            best = separation, support
    if best is None:
        raise ArithmeticError("mixed native mark could not be separated")
    _, support = best
    # The signed maximal minors of a rectangular matrix form an exact null vector.
    direction = [arb(0)] * columns
    for local, column in enumerate(support):
        minor = arb_mat([[row[i] for i in support if i != column] for row in rows]).det()
        direction[column] = (-1)**local * minor
    amplitude = max(abs(x).upper() for x in direction)
    return [x / amplitude for x in direction], support


def recover_empty_gram(probes, radii_squared, step):
    """Invert unnormalized empty probes with independent early and late norm insertions."""
    nodes = [(-step*u).exp() for u in radii_squared]
    basis = arb_mat([[node**i for node in nodes] for i in range(len(nodes))])
    inverse = basis.inv()
    return inverse * probes * inverse.transpose()


def replay(root):
    domain = root / NATIVE
    fixture = json.loads(domain.read_bytes())
    source_refs(root, fixture["source_originals"])
    radii = ("1/4", "3/10", "7/20", "2/5", "11/20")
    n = len(radii)
    edges = list(combinations(range(n), 2))
    with ctx.workprec(256):
        fields = list(map(ball, radii))
        u = [x*x for x in fields]
        momentum, width = ball("1/4"), ball("2/5")
        spectra = [(name != "Ye", squared_singular_values(fixture["couplings"][name]))
                   for name in ("Yu", "Yd", "Ye")]
        ref, _ = conditional_word(spectra, arb(0), arb(0), momentum, width)
        kernels = [[conditional_word(spectra, x, y, momentum, width) for y in fields] for x in fields]
        k0 = [[value/ref for value, _ in row] for row in kernels]
        ky = [[value*mark/ref for value, mark in row] for row in kernels]
        if not ref > 0 or any(not k0[i][j] > 0 or not k0[i][j].overlaps(k0[j][i])
                              for i in range(n) for j in range(n)):
            raise ArithmeticError("invalid native unmarked kernel or fixed empty reference")

        def coefficient(form, i, j):
            return form[i][j] + form[j][i] - form[i][i] - form[j][j]

        weighted_marginals = [[arb(0) for _ in edges] for _ in fields]
        for edge, (i, j) in enumerate(edges):
            weighted_marginals[i][edge] = k0[i][j] - k0[i][i]
            weighted_marginals[j][edge] = k0[j][i] - k0[j][j]
        y_row = [coefficient(ky, i, j) for i, j in edges]
        mixed_kernel = [[u[i]*ky[i][j] for j in range(n)] for i in range(n)]
        mixed = [coefficient(mixed_kernel, i, j) for i, j in edges]
        direction, support = independent_direction([*weighted_marginals, y_row], mixed)
        if any(not sum(x*y for x, y in zip(row, direction)).contains(0) for row in (*weighted_marginals, y_row)):
            raise ArithmeticError("native constraint closure did not replay")
        delta = [[arb(0) for _ in range(n)] for _ in range(n)]
        for (i, j), value in zip(edges, direction):
            delta[i][j] = delta[j][i] = value
            delta[i][i] -= value
            delta[j][j] -= value

        unit = arb(1) / (n*(n+1))
        step = unit / 16
        if any(not step * sum(abs(value) for value in row) < unit / 2 for row in delta):
            raise ArithmeticError("symmetric perturbation exceeded the PSD margin")
        bases = [[[unit*(1 + (i == j)) + sign*step*delta[i][j] for j in range(n)] for i in range(n)]
                 for sign in (1, -1)]

        def read(form, mass):
            return sum(mass[i][j]*form[i][j] for i in range(n) for j in range(n))

        u_kernel = [[u[i]*k0[i][j] for j in range(n)] for i in range(n)]
        marks = [[read(form, mass) for form in (k0, u_kernel, ky, mixed_kernel)] for mass in bases]
        gap = 2 * step * sum(x*y for x, y in zip(direction, mixed))
        if (gap.contains(0) or not (marks[0][3]-marks[1][3]).overlaps(gap) or
                any(not a.overlaps(b) for a, b in zip(marks[0][:3], marks[1][:3]))):
            raise ArithmeticError("fitted marks or withheld mixed gap did not replay")
        # Each entry of the empty Gram is K0_ij times the corresponding mass.
        # Simulate independent quadratic insertions on the two window norms.
        # The inverse itself only receives the unmarked samples.
        probe_step = ball("8")
        for mass, values in zip(bases, marks):
            gram = [[mass[i][j]*k0[i][j] for j in range(n)] for i in range(n)]
            probes = arb_mat([[sum(gram[i][j]*(-probe_step*(left*u[i]+right*u[j])).exp()
                                   for i in range(n) for j in range(n))
                               for right in range(n)] for left in range(n)])
            recovered = recover_empty_gram(probes, u, probe_step)
            if any(not recovered[i,j].overlaps(gram[i][j]) for i in range(n) for j in range(n)):
                raise ArithmeticError("unmarked quadratic insertion inverse did not recover the joint Gram")
            readout = sum(recovered[i,j]*u[i]*(ky[i][j]/k0[i][j]) for i in range(n) for j in range(n))
            if not readout.overlaps(values[3]):
                raise ArithmeticError("empty Gram did not recover the mixed mark")
        return {"status":"EXACT_CONDITIONAL_R032_REFLECTED_MEASURE_AMBIGUITY", "source":str(WINDOWS),
                "native_domain_sha256":digest(domain), "radii":radii, "histories":n*n,
                "fit_constraint_rank":len(weighted_marginals)+1,
                "rank_with_mixed_mark":len(weighted_marginals)+2,
                "separating_edges":[edges[i] for i in support],
                "same_bare_and_weighted_radial_marginals":True, "both_reflected_grams_positive":True,
                "reference":"same free charged R032 reference at both ends",
                "marks_over_reference":[dict(zip(("Z0","ZU","ZY","ZUY"), [x.str(15) for x in v])) for v in marks],
                "nonzero_ZUY_difference":gap.str(15),
                "connected_difference":(gap/marks[0][0]).str(15),
                "whole_empty_gram_recovers_mixed_mark":True,
                "independent_empty_quadratic_probes":n*(n+1)//2,
                "empty_probe_step":str(probe_step), "held_reference_across_probes":True,
                "scope":"conditional two-window charged84 sector, common colour-centre projected trace",
                "physical_R113_source_derived":False}
