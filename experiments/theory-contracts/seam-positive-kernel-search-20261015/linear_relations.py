"""Exact, declared linear relations among moments of the existing word algebra."""
from fractions import Fraction
import sympy as sp
from word_algebra import WordAlgebra


def _q(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise ValueError('Relation coefficients must be integers or rational strings')
    try:
        result = Fraction(value)
        float(result)  # numerical search backend must be able to represent it
        return result
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        raise ValueError('Invalid finite rational relation coefficient') from exc


def normalize_linear_relations(m, degree, relations=None):
    if relations is None:
        return []
    if not isinstance(relations, list):
        raise ValueError('linear_relations must be a list')
    algebra = WordAlgebra(m)
    normalized = []
    for row in relations:
        if not isinstance(row, dict) or not {'terms', 're'} <= set(row) or set(row) - {'terms', 're', 'im'}:
            raise ValueError('Relation requires terms and re, optionally im')
        if not isinstance(row['terms'], list):
            raise ValueError('Relation terms must be a list')
        combined = {}
        for term in row['terms']:
            if not isinstance(term, dict) or not {'word', 're'} <= set(term) or set(term) - {'word', 're', 'im'}:
                raise ValueError('Each term requires word and re, optionally im')
            if not isinstance(term['word'], (list, tuple)):
                raise ValueError('Relation word must be a token list')
            sign, word = algebra.reduce_word(term['word'])
            if len(word) > 2 * degree:
                raise ValueError('Relation word exceeds available moment degree')
            re, im = combined.get(word, (Fraction(0), Fraction(0)))
            combined[word] = (re + sign*_q(term['re']), im + sign*_q(term.get('im', 0)))
        terms = [{'word': list(word), 're': str(re), 'im': str(im)}
                 for word, (re, im) in sorted(combined.items(), key=lambda item: (len(item[0]), item[0]))
                 if re or im]
        normalized.append({'terms': terms, 're': str(_q(row['re'])), 'im': str(_q(row.get('im', 0)))})
    return normalized


def check_linear_relations(relations, moments, m, degree):
    """Reject missing moments or a nonzero exact residual, without float tolerances."""
    rows = normalize_linear_relations(m, degree, relations)
    for index, row in enumerate(rows):
        value = sp.Integer(0)
        for term in row['terms']:
            word = tuple(term['word'])
            if word not in moments:
                raise ValueError(f'Linear relation {index}: missing exact word {word}')
            coefficient = sp.Rational(term['re']) + sp.I*sp.Rational(term['im'])
            value += coefficient*moments[word]
        target = sp.Rational(row['re']) + sp.I*sp.Rational(row['im'])
        if sp.expand(value-target) != 0:
            raise ValueError(f'Linear relation {index} fails exact replay')
    return len(rows)
