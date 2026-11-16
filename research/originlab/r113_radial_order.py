"""Rigorous association checks on the original conditional Full96 radial chart."""

from itertools import combinations, product
from pathlib import Path
import json

from flint import arb, ctx

from reconstruction import require
from verify import digest


RADIAL_DATA = Path("docs/research/fundamental-plan/evidence/R113-FULL96-RADIAL-BALLS-20261116.json")


def controls(root):
    path = root/RADIAL_DATA
    data = json.loads(path.read_text())
    require(data["input_role"]=="conditional_R032_original_full96_including_majorana",
            "radial test requires the reviewed conditional chart")
    size = len(data["radii"])
    indices = range(size)
    pairs = tuple(combinations(indices,2))
    with ctx.workprec(160):
        u = [arb(r)**2 for r in data["radii"]]
        kernel = [[arb(v) for v in row] for row in data["reflected_kernel_over_conditional_reference"]]
        late = [[arb(v) for v in row] for row in data["late_mark_over_conditional_reference"]]
        require(size==5 and all(len(m)==size and all(len(row)==size for row in m) for m in (kernel,late)),
                "expected the complete five-radius chart")
        require(all(u[i]<u[k] for i,k in pairs) and all(v>0 for row in kernel for v in row),
                "positive kernel in increasing U order required")
        minor_ratios = [kernel[i][j]*kernel[k][l]/(kernel[i][l]*kernel[k][j])-1
                        for (i,k),(j,l) in product(pairs,repeat=2)]
        y = [[value/weight for value,weight in zip(row,base)]
             for row,base in zip(late,kernel)]
        early_drops = [y[i][j]-y[i+1][j] for i,j in product(range(size-1),indices)]
        late_drops = [y[i][j]-y[i][j+1] for i,j in product(indices,range(size-1))]
        require(all(v>0 for v in (*minor_ratios,*early_drops,*late_drops)),
                "strict TP2 and strict decrease of the normalized Y mark required")
        prior_cross_ratio_floor = arb(243)/250
        require(all(prior_cross_ratio_floor*(1+v)>1 for v in minor_ratios),
                "the fermion kernel must compensate the certified amount of negative association")
        # The two-site angular integral on S^3 is 2 I_1(z)/z.
        # Its positive series in (r*s)^2 certifies TP2 for every positive
        # link stiffness. Unit stiffness here is only a numerical replay.
        radii = [arb(r) for r in data["radii"]]
        angular = []
        for r in radii:
            row = []
            for s in radii:
                argument = r*s
                row.append(2*argument.bessel_i(1)/argument)
            angular.append(row)
        angular_minors = [angular[i][j]*angular[k][l]-angular[i][l]*angular[k][j]
                          for (i,k),(j,l) in product(pairs,repeat=2)]
        require(all(v>0 for v in angular_minors), "conditional O(4) radial link must have positive minors")

        def connected(weights):
            z,zu,zy,zuy = (sum((weight*entry(i,j) for (i,j),weight in weights.items()),arb(0))
                           for entry in (lambda i,j: kernel[i][j],
                                         lambda i,j: u[i]*kernel[i][j],
                                         lambda i,j: late[i][j],
                                         lambda i,j: u[i]*late[i][j]))
            return zuy/z-zu*zy/z**2

        independent = connected({(i,j):arb(1)/25 for i,j in product(indices,repeat=2)})
        a,b = 0,size-1
        supported = tuple(product((a,b),repeat=2))
        paired = connected({(i,j):arb(99 if i==j else 1)/200 for i,j in supported})
        crossed = connected({(i,j):arb(1 if i==j else 99)/200 for i,j in supported})
        require(independent<0 and paired<0 and crossed>0,
                "uniform, positively and negatively associated controls must separate the signs")
        return {"status":"CERTIFIED_CONDITIONAL_RADIAL_ORDER","data_sha256":digest(path),
                "positive_two_by_two_minors":len(minor_ratios),
                "smallest_cross_ratio_excess_lower":str(min(v.lower() for v in minor_ratios)),
                "sufficient_positive_boson_cross_ratio_floor":"243/250",
                "local_O4_radial_link":{"strict_positive_minors_at_control_stiffness_one":len(angular_minors),
                                        "formula":"2 I_1(a*r*s)/(a*r*s), all positive a",
                                        "premises":"open scalar nearest-neighbor links and radial site factors"},
                "smallest_mark_drop_lower":{"early":str(min(v.lower() for v in early_drops)),
                                            "late":str(min(v.lower() for v in late_drops))},
                "connected_UY_under_strictly_positive_TP2_boson_weights":"strictly negative",
                "controls":{"uniform_independent":str(independent),"positive_association_same_bare_marginals":str(paired),
                            "negative_association_same_bare_marginals":str(crossed)},
                "scope":"R032 spatial two-window chart with frozen sea for Y; TP2 is a premise on the boson weights",
                "original_native_time_or_measure_identified":False}
