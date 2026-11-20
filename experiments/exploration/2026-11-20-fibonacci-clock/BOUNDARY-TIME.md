# Boundary source loophole: nontrivial modular clock from the *same* 2+3 splitting

## Difference between two observable algebras

The quantum braid image is irreducible in each four-anyon total-charge sector, giving dimensions 2 and 3. Any full-braid invariant faithful density on their direct sum has the form

\[
\rho=\frac{p_2}{2}I_2\oplus\frac{p_3}{3}I_3,\qquad
p_2,p_3>0,\quad p_2+p_3=1.
\]

For the **charge-conserving algebra** \(\mathcal A_0=M_2(\mathbb C)\oplus M_3(\mathbb C)\), this rho is central, and its Tomita modular flow is the identity. It does **not** generate a time direction.

If a physically constructed *boundary source* adds permitted inter-sector transitions to make the larger algebra \(\mathcal A_{\rm ext}=M_5(\mathbb C)\), rho is no longer central, unless \(p_2/2=p_3/3\). The modular evolution of a charge-changing matrix unit E between the blocks is

\[
\sigma_t^\rho(E)=e^{it\log(\rho_2/\rho_3)}E,\qquad
\omega=\log\left(\frac{\rho_3}{\rho_2}\right)
=\log\left(\frac{p_3/3}{p_2/2}\right).
\]

The modular log-density is **exactly a scalar plus the relative U(1) sector charge**:

\[
\log\rho=a I_5+\beta Y,\qquad
Y=\operatorname{diag}(1/2,1/2,-1/3,-1/3,-1/3),\qquad
\beta=\frac65\log\left(\frac{p_2/2}{p_3/3}\right).
\]

Thus relative modular time and relative charge share the **same generator** *if inter-sector physical observables exist*. This is a general fact of a two-block algebra, not a unique signature of Fibonacci.

## Exactly tested weight choices

| Boundary state choice | Sector probabilities p2,p3 | Modular gap omega |
|---|---|---|
| Uniform among five fusion channels | 2/5, 3/5 | 0 |
| Equal sector probabilities | 1/2, 1/2 | log(2/3) |
| Quantum dimensions squared (proposed canonical topological weight) | 1/(1+phi²), phi²/(1+phi²) | log(1+sqrt(5)/3) ≈ +0.556958542 |
| Dimension weighted by quantum dimension | 2/(2+3phi), 3phi/(2+3phi) | log(phi) ≈ +0.481211825 |

All densities are braid invariant. Therefore **the modular rate is NOT selected by braiding and fusion alone**. Assigning the quantum-dimension-squared state is a new state-selection postulate whose physical justification is open. Even under that assumption, the modular flow only rotates charge-exchange operators; it is NOT full spacetime dynamics.

Test with SymPy:

    python modular_boundary_clock.py

## Physical interpretation, explicitly conditional

A local boundary charge reservoir could allow both 2 and 3 blocks as physical source states while conserving total charge including the reservoir. This is a plausible *operational boundary-to-state/time interface* which the original closed finite braiding system lacks. It also makes the formerly unused relative U(1) a possible modular source action. However it has no original weighted TFPT P1 kernel, no OS seam, no absolute c3, no physical clock calibration and no 3+1 Lorentz continuum.

**Status unchanged: Stufe 1, not F1 or physical time.**
