"""Compare the flag-return channel with independently imported TFPT word projectors."""
import argparse
import importlib.util
from pathlib import Path

import sympy as sp
from flag_process import I, N, change


def load_bridge(root):
    checker = root / "experiments/theory-contracts/compiler-clifford-bridge/checker.py"
    spec = importlib.util.spec_from_file_location("tfpt_clifford_bridge", checker)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    return bridge


def main(root):
    bridge = load_bridge(root)
    src, _ = bridge.inherited(root)
    q = bridge.finite_data(src)["q"]
    generators = bridge.generators()
    native = sp.zeros(16)
    for v in src.W16:
        word = bridge.monomial(v, generators)
        column = word.vec() / 2
        factor = sp.Rational(1, 3) if q[v] == 0 else sp.Rational(2, 3)
        if not any(v):
            factor = 1
        native += factor * column * column.H

    a, b = change(0), change(1)
    galleries = [I, a, b, a*b, b*a, a*b*a]
    masks = [m.applyfunc(lambda x: int(x != 0)) for m in galleries]
    assert sum(masks, sp.zeros(N)) == sp.ones(N)
    weights = [sum(mask.multiply_elementwise(m), sp.S.Zero) / N**2
               for mask, m in zip(masks, galleries)]
    assert weights == [sp.Rational(1, 21)] * 6
    conditional = [p / sum(weights) for p in weights]
    x = sp.Matrix([[0, 1], [1, 0]])
    z = sp.diag(1, -1)
    y = sp.I * x * z
    gamma = [sp.kronecker_product(u, v)
             for u, v in [(sp.eye(2), x), (sp.eye(2), z), (x, y), (z, y), (y, y)]]
    channels = [sp.eye(16)] + [(sp.eye(16) + sp.kronecker_product(g.conjugate(), g))/2
                                for g in gamma]
    flag = sum((p*k for p, k in zip(conditional, channels)), sp.zeros(16))
    differences = sum(x != 0 for x in flag - native)
    assert differences == 0
    print("Independently imported q* word projectors vs flag-return channel:", 256-differences, "/ 256")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[4])
    main(parser.parse_args().repo)
