"""Exact controls for the source recognition theorem, not a P1 source search.

The CAR matrix units are imported from the existing original-source checker.
All-mode uniqueness is proved in DERIVATION.txt, not extrapolated from samples.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from math import factorial
from pathlib import Path
import importlib.util
import json


ROOT = Path(__file__).resolve().parents[3]
OLD = ROOT / 'experiments/theory-contracts/original-admissible-source-closure-20261025/check_w8_stress_cyclic.py'
spec = importlib.util.spec_from_file_location('existing_original_car_units', OLD)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
LOW, HIGH = -48, 48
SEA = (1 << -LOW) - 1


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


@lru_cache(None)
def current_on_mask(mask, channel, n):
    """j_(channel,n) = sum_l a*_(8l+channel) a_(8(l+n)+channel)."""
    out = {}
    for target in range(LOW + ((channel - LOW) % 8), HIGH, 8):
        source = target + 8*n
        if LOW <= source < HIGH:
            moved = old.fermion_move(mask, source-LOW, target-LOW)
            if moved is not None:
                bits, sign = moved
                out[bits] = out.get(bits, 0) + sign
    return tuple(sorted((bits, c) for bits, c in out.items() if c))


def act(v, channel, n):
    out = {}
    for mask, coefficient in v.items():
        for bits, sign in current_on_mask(mask, channel, n):
            out[bits] = out.get(bits, 0) + coefficient*sign
    return {bits: c for bits, c in out.items() if c}


@lru_cache(None)
def car_word(word):
    if not word:
        return ((SEA, 1),)
    return tuple(sorted(act(dict(car_word(word[1:])), *word[0]).items()))


@lru_cache(None)
def vacuum_recursion(word):
    if not word:
        return 1
    channel, n = word[0]
    if n < 0:
        return 0
    return sum(n*vacuum_recursion(word[1:j]+word[j+1:])
               for j in range(1, len(word))
               if word[j] == (channel, -n))


def energy(mask, alpha):
    return sum((F(r)+alpha)*(((mask >> (r-LOW)) & 1)-int(r < 0))
               for r in range(LOW, HIGH))


@lru_cache(None)
def thermal_recursion(word, q):
    if not word:
        return F(1)
    channel, n = word[0]
    coefficient = F(n)/(1-q**n)
    return sum((coefficient*thermal_recursion(word[1:j]+word[j+1:], q)
                for j in range(1, len(word))
                if word[j] == (channel, -n)), F(0))


@lru_cache(None)
def normal_order(word):
    """Independent CCR sorting, then exact infinite Gibbs trace below."""
    for j in range(len(word)-1):
        a, b = word[j], word[j+1]
        # All creators, sorted by mode, precede all annihilators.
        key = lambda x: (int(x[1] > 0), x[0], abs(x[1]))
        if key(a) > key(b):
            out = dict(normal_order(word[:j]+(b, a)+word[j+2:]))
            if a[1] > 0 and b == (a[0], -a[1]):
                for monomial, c in normal_order(word[:j]+word[j+2:]):
                    out[monomial] = out.get(monomial, 0)+a[1]*c
            return tuple((m, c) for m, c in out.items() if c)
    return ((word, 1),)


def infinite_gibbs_trace(word, q):
    total = F(0)
    for monomial, coefficient in normal_order(word):
        creators = Counter((a, -n) for a, n in monomial if n < 0)
        annihilators = Counter((a, n) for a, n in monomial if n > 0)
        if creators != annihilators:
            continue
        value = F(coefficient)
        for (_, n), k in creators.items():
            # Exact sum (1-q^n) sum_l q^(nl) l!/(l-k)!, no cutoff.
            value *= n**k*factorial(k)*(q**n/(1-q**n))**k
        total += value
    return total


def main():
    words = {()}
    for a, b, n, m in product(range(8), range(8), (1, 2), (1, 2)):
        words.update(permutations(((a, n), (a, -n), (b, m), (b, -m))))
    for colors in ((0, 0, 0), (0, 0, 1), (0, 1, 2)):
        word = tuple(item for a in colors for item in ((a, 1), (a, -1)))
        words.update(permutations(word))
    for word in sorted(words):
        need(dict(car_word(word)).get(SEA, 0) == vacuum_recursion(word),
             'original folded CAR moment versus recognition recursion')

    inputs = [{SEA: 1}]+[act({SEA: 1}, a, -1) for a in range(8)]
    commutators = 0
    labels = tuple(product(range(8), (-2, -1, 1, 2)))
    for v in inputs:
        for (a, n), (b, m) in product(labels, repeat=2):
            left = act(act(v, b, m), a, n)
            for mask, c in act(act(v, a, n), b, m).items():
                left[mask] = left.get(mask, 0)-c
            left = {mask: c for mask, c in left.items() if c}
            right = {mask: n*c for mask, c in v.items()} if a == b and n == -m else {}
            need(left == right, 'actual folded current CCR on original CAR inputs')
            commutators += 1

    energy_checks = 0
    for alpha in (F(1, 4), F(3, 4)):
        for v in inputs:
            for a, n in labels:
                for mask, c in v.items():
                    for bits, _ in current_on_mask(mask, a, n):
                        need(energy(bits, alpha)-energy(mask, alpha) == -8*n,
                             'same original energy, not a replacement time')
                        energy_checks += 1

    thermal_words = {()}
    for a, b, n, m in product(range(2), range(2), (1, 2), (1, 2)):
        thermal_words.update(permutations(((a, n), (a, -n), (b, m), (b, -m))))
    thermal_words.update(permutations(((0, 1), (0, -1))*3))
    q = F(1, 3)  # q = exp(-8 beta), a diagnostic exact temperature.
    for word in sorted(thermal_words):
        need(thermal_recursion(word, q) == infinite_gibbs_trace(word, q),
             'KMS recursion versus independently normal-ordered infinite trace')

    # An elementary condition countercontrol: stationary is weaker than ground.
    # rho_A=|1><1|, rho_B=(|0><0|+|2><2|)/2 have same <N>, different <a*2 a2>.
    rho_a, rho_b = {1: F(1)}, {0: F(1, 2), 2: F(1, 2)}
    moment = lambda rho, k: sum((p*F(factorial(n), factorial(n-k))
                                for n, p in rho.items() if n >= k), F(0))
    need(moment(rho_a, 1) == moment(rho_b, 1) == 1, 'equal oscillator covariance')
    need(moment(rho_a, 2) == 0 and moment(rho_b, 2) == 1, 'distinct fourth moments')

    result = {
        'status': 'PASS_SCOPED',
        'original_folded_CAR_vacuum_word_comparisons': len(words),
        'original_folded_CAR_CCR_vector_checks': commutators,
        'original_energy_transition_checks_both_quarters': energy_checks,
        'thermal_word_comparisons_to_exact_infinite_trace': len(thermal_words),
        'CAR_window': '[−48,48): six levels per actual folded channel',
        'checked_word_bound': 'at most six mode-energy units; no sampled boundary leakage',
        'stationary_same_covariance_fourth_moments': ['0', '1'],
        'folded_neutral_quarter_clock': 'identity: H=8N implies exp(-i*pi*H/2)=1',
        'all_mode_proof': 'analytic representation/word recursion in DERIVATION.txt',
        'P1_local_source_image_derived': False,
        'complete_TFPT_solution': False,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
