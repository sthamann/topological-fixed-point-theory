"""Certified negative mixing for reflected radial histories beyond TP2."""

from itertools import combinations, product
import json

from flint import arb, ctx
import sympy as sp

from r113_radial_order import RADIAL_DATA
from reconstruction import require
from verify import digest


RADIAL_PROFILE_SHA256 = "7607ca22beea3e5e7a19f29138e6d73935ec8af5c44a88e221e4f97b4c62f3a4"
RADII = tuple(map(sp.Rational,("1/4","3/10","7/20","2/5","11/20")))
LATE_BASELINES = (47,48,49,51)


def moment_bound(early, late, early_square, late_square, cross, connected=None):
    """Bound the mixing on the reviewed profile from five relative moments."""
    ve,vl = early_square-early**2,late_square-late**2
    covariance = cross-early*late
    require(all(v.is_Rational for v in (early,late,early_square,late_square,cross))
            and RADII[0]**2<=early<=RADII[-1]**2 and RADII[0]**2<=late<=RADII[-1]**2
            and ve>=0 and vl>=0 and covariance**2<=ve*vl,
            "expected consistent real radial moments on the reviewed interval")
    baseline = next(b for r,b in zip(RADII[1:],LATE_BASELINES) if early<=r**2)
    upper = -(sp.Rational(77,10))*ve-baseline*covariance+7*sp.sqrt(ve*(vl+(late-early)**2))
    return {"early_variance":str(ve),"late_variance":str(vl),"radial_covariance":str(covariance),
            "baseline":baseline,"conditional_connected_UY_upper_bound":str(upper),
            "strict_negative_bound":bool(upper<0),
            "supplied_response_exceeds_bound":bool(connected>upper)
                if connected is not None and connected.is_real else None,
            "physical_source_identity_proved":False}


def controls(root):
    path = root / RADIAL_DATA
    require(digest(path)==RADIAL_PROFILE_SHA256,"review the changed Full96 radial profile before applying its bounds")
    data = json.loads(path.read_text())
    with ctx.workprec(160):
        u = [arb(r)**2 for r in data["radii"]]
        late = data["late_mark_over_conditional_reference"]
        empty = data["reflected_kernel_over_conditional_reference"]
        y = [[arb(l) / arb(k) for l, k in zip(lr, kr)] for lr, kr in zip(late, empty)]
        size = len(u)
        indices = range(size)
        require(size == 5 and all(len(row) == size for row in y) and len(y) == size,
                "five common early and late values required")
        require(all(u[i] < u[i+1] for i in range(size-1)), "radial values must be strictly increasing")

        early_slope = arb(77)/10
        error_slope = arb(7)
        require(all(y[i][j]-y[i+1][j] > early_slope*(u[i+1]-u[i])
                    for i,j in product(range(size-1),indices)), "early slope enclosure failed")

        # At a mean inside cell k, the anchored backward secant in the
        # second variable is a convex combination of its grid endpoint
        # secants in rows k and k+1. Zero-length secants contribute zero
        # weight. This covers the whole continuous mean cell.
        mean_cells = []
        for k, baseline in enumerate(LATE_BASELINES):
            endpoints = (k,k+1)
            secants = [(y[i][h]-y[i][j])/(u[j]-u[h])
                       for i,h,j in product(endpoints,endpoints,indices) if j!=h]
            require(all(abs(v-baseline) < error_slope for v in secants),
                    "anchored late slope enclosure failed")
            mean_cells.append({"left_radius":data["radii"][k],"right_radius":data["radii"][k+1],
                               "baseline":baseline,"tested_endpoint_secants":len(secants),
                               "largest_residual_upper":str(max(abs(v-baseline).upper() for v in secants))})

        def evaluate(weights):
            scale = sum(sum(row) for row in weights)
            mu = sum(weights[i][j]*u[i] for i,j in product(indices,repeat=2))/scale
            variance = sum(weights[i][j]*(u[i]-mu)**2 for i,j in product(indices,repeat=2))/scale
            cross = sum(weights[i][j]*(u[i]-mu)*(u[j]-mu)
                        for i,j in product(indices,repeat=2))/scale
            connected = sum(weights[i][j]*(u[i]-mu)*y[i][j]
                            for i,j in product(indices,repeat=2))/scale
            return variance,cross,connected

        # G = 11^t + aa^t is positive semidefinite and strictly positive
        # entrywise. The nonmonotone a makes an ordered minor negative.
        a = (0,2,1,3,2)
        gram = [[1+a[i]*a[j] for j in indices] for i in indices]
        ordered_pairs = tuple(combinations(indices,2))
        negative_minors = sum(gram[i][j]*gram[k][l] < gram[i][l]*gram[k][j]
                              for (i,k),(j,l) in product(ordered_pairs,repeat=2))
        variance,cross,connected = evaluate(gram)
        require(negative_minors and cross>0 and connected < -(early_slope-error_slope)*variance,
                "non-TP2 positive Gram control failed")
        reversed_endpoints = [[int((i,j) in ((0,size-1),(size-1,0))) for j in indices] for i in indices]
        _,negative_cross,positive_connected = evaluate(reversed_endpoints)
        require(negative_cross<0 and positive_connected>0,"negative covariance control must change the sign")
        return {"status":"CERTIFIED_CONDITIONAL_REFLECTED_RADIAL_BOUND",
                "data_sha256":digest(path),"early_secant_slope_lower":"77/10",
                "late_anchored_slope_error_upper":"7","mean_cells":mean_cells,
                "bound":"C_UY <= -(7/10) Var(U_early) - b Cov(U_early,U_late), b in {47,48,49,51}",
                "nonnegative_radial_covariance_implies_strict_negative_mixing_for_nonzero_variance":True,
                "reflected_gram_control":{"negative_ordered_minors":negative_minors,"variance":str(variance),
                                          "radial_covariance":str(cross),"connected_UY":str(connected)},
                "anticorrelated_control":{"radial_covariance":str(negative_cross),
                                           "connected_UY":str(positive_connected)},
                "scope":"nonnegative symmetric selected five-radius weights and the reviewed normalized Y profile",
                "physical_source_identity_proved":False}
