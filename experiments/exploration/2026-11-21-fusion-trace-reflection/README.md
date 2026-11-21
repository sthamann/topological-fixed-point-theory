# TFPT exploration — 2026-11-21

Independent search across graph operator algebras, fusion rules, reflection-positive Markov processes, and categorical anomalies.

**New outcome:** The apparent conflict between the finite 5-dimensional braid-invariant tracial state and the infinite PF/KMS history state disappears when the latter is **conditioned on the same terminal fusion charge**. Exactly one compatible stationary AF trace exists under the chosen primitive graph (unique positive PF ray). On level 5 / terminal tau it reduces to I5/5 and a preceding 2/5 + 3/5 split.

A reflected two-step Doob transition has spectrum {1,phi^-4}, giving the finite positive generator H=4 log(phi)(I-Pi) without fitting the rate. This is a 0+1D process, not original P1 or physical relativity.

**Critical no-gos:** The CK gauge/KMS flow is identity on the finite AF/braid observables, so the two times have not been identified. Fib itself and its same-oriented tensor powers cannot be condensed by a nontrivial connected ribbon étale algebra into the E8 vacuum (Booker-Davydov theorem). A family of primitive n-state graphs with polynomials x^n-x-1 has growth eigenvalues approaching 1, defeating minimum-entropy selection without a separately selected binary alphabet.

## Reproduce

```sh
python check_new.py
```

Expected exact test ending:

```text
AF level-n tr_n = phi^-n(Tr X0+phi Tr Xtau) compatible n=1..16
Conditioned n=5 total tau: M5 trace I5/5, preceding 4-event sectors 2/5,3/5
P eigenvalues 1,-phi^-2; P-adjoint P=P^2 eigen 1,phi^-4
Baseline 3/16 labeled 2-state graphs primitive, 2/16 invertible; both Fib
ALL EXACT CHECKS PASSED
```

Detailed reports: IDEAS.md and VERDICT.md. Original P1, c3, P2 physical fermions, 3+1, gravity **unresolved**. No ledger or Atlas updates.
