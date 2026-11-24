# CODE und tatsächliche Ausgabe

Rechenstand vom 20. November 2026. Python 3.14.3, SymPy 1.14.0. Alle unten aufgeführten Aufrufe wurden erfolgreich beendet. Alle Kernvergleiche beruhen auf exakter Arithmetik. Dezimalangaben veranschaulichen die daneben hergeleiteten Werte. Für den unabhängig importierten Vergleich bezeichnet TFPT_REPO das Stammverzeichnis des Repositorys. Der geprüfte Originalcheckout war a1273187b82488d34d2a474470152a41d42d4055.

## Hauptbefunde

Die Dateien coxeter_sichtbarkeit.py und dihedral.py kontrollieren den Auswahlsatz. markenwirkung.py kontrolliert die Gleichverteilung der drei effektiven Beschreibungen. originalvergleich.py vergleicht die Flaggenmatrix mit unabhängig geladenen Originalwörtern. Die universellen Ausschlüsse stützen sich auf die Beweise in den begleitenden Markdown Dateien.

## apartment.py

```sh
python3 apartment.py
```

```text
flags, ordered opposite pairs, apartments: 21 168 28
opposite pairs per apartment / apartments per u / pairs per u: 6 4 24
longest successful return, each conditional pair, each conditional u: 1/21 1/168 1/7
fixed initial flag: four u outside initial line, conditional probabilities 1/4 each
distinct embedded A3/D5 from this process: 7 with 12/40 roots
averaged spectrum check / constant old points / six zero-sum old directions: 0, 0, 1/2
covariant under all 168 linear maps: OK
```

## continuous_time.py

```sh
python3 continuous_time.py
```

```text
Native half step normalized grade-two Choi weight: -sqrt(6)/24 + sqrt(3)/48 + 1/16
Half-step-only lower bound on two-step error: 28/9 - 16*sqrt(3)/9 0.031909675433107035
All intermediate times lower bound: -4/9 + 2*3**(3/4)/9 0.06211267932328392
All-times CPTP / Euclidean Gram family, grade-two return: 5/9 ; two-step excess: 1/9
Static five-choice reversible CP dilation, two-step response: 2/9,11/18
Fresh-choice suspension, scalar reflected 2x2 determinant: -1/36
Two local balanced/nonbalanced Hermitian jumps: decay ratio >=3; native log(3)/log(3/2)<3
```

## coxeter_sichtbarkeit.py

```sh
python3 coxeter_sichtbarkeit.py
```

```text
A2 q=2: raw averaged weights [1,4,4,4,4,4]; active 1/5; idle 1/21
A2 q=3: conjugation invariant; individual active weights 11,11,11,9,9
Uniform dihedral control (3<=m<=18, 2<=q<=8): [(3, 2)]
Unequal B2/I2(4) q_s=2,q_t=3: central longest has unique weight 36
Rank three: A3 excluded by positive P(-2); B3,H3 excluded by order
Rank >=4: 5!-1 > 3^4; subsequent order lower bound grows by >=6 per rank
Reducible A1 x A1: active weights (2,3,6)
A1: one active label, criterion vacuous
Raw empty/Lueders channel sector values: 1,5/21,13/21
```

## dihedral.py

```sh
python3 dihedral.py
```

```text
I2(2), parameters (2,2): total 9, distinct active class means [Fraction(2, 1), Fraction(4, 1)], equalizable False
I2(3), parameters (2,2): total 21, distinct active class means [Fraction(4, 1)], equalizable True
I2(3), parameters (3,3): total 52, distinct active class means [Fraction(9, 1), Fraction(11, 1)], equalizable False
I2(4), parameters (2,4): total 135, distinct active class means [Fraction(8, 1), Fraction(10, 1), Fraction(17, 1), Fraction(64, 1)], equalizable False
I2(6), parameters (2,2): total 189, distinct active class means [Fraction(4, 1), Fraction(14, 1), Fraction(16, 1), Fraction(64, 1)], equalizable False
parameter-system control: 1 / 175 [(3, 2, 2)] (includes formal types without asserted existence)
class weights 1/3 each suffice even with nonuniform six chart weights
active endpoint mass each: 4/21 ; conditional on being active: 1/5
```

## fano_lattice.py

```sh
python3 fano_lattice.py
```

```text
punctured and extended weights: [(0, 1), (3, 7), (4, 7), (7, 1)] [(0, 1), (4, 14), (8, 1)]
doubly even linear one-bit extensions: 1 / 16
panel controls with seven distinct 3-subsets on seven points: 30 / 6724520
Construction-A covolume, roots, simple rank, E8 arm lengths: 1 240 8 [1, 2, 4]
3-coordinate root counts over all 56 sections: {6: 56}
after choosing a pair in a tetrad, roots in A3 and orthogonal D5: 12 40
Smith divisors for A3, complement and direct sum in E8: (1, 1, 1) (1, 1, 1, 1, 1) (1, 1, 1, 1, 1, 1, 1, 4)
Z4 glue class 0,1,2,3 root counts: [52, 64, 60, 64]
automorphism group of Fano incidence alone, invariant A3 count: 168 0
all embedded A3, S3-invariant, and invariant including S3-commuting signs: 7560 18 0
```

## flag_process.py

```sh
python3 flag_process.py
```

```text
Flags 21 ; exact Hecke and braid relations; Burau parameter 1/2
Characteristic polynomial of mean of square transfers:
(x - 1)*(4*x - 1)**8*(64*x**2 - 80*x + 23)**6/4503599627370496
Normalized double exchange overlap, independent of r: 2/3
Reflected squared distance a²,b²: 1/1024
Reflected squared distance ab,ba: 1/512
Same scalar observation 1/32 1/32 ; after a 1/256 0
Gram of six positive basis words: [1, 1/2, 1/2, 1/4, 1/4, 1/8] ; OS dimension: 6
Pure word a²b^-² on invariant rank two: has nontrivial standard block
Raw relative position probabilities: [1/21, 2/21, 2/21, 4/21, 4/21, 8/21]
Return rate: 2/7 ; recorded conditional weights: [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
Five visible orthogonal modes by parity: 3 odd and 2 even
Return-conditioned M4 channel: 256/256; without return selection: different
Elementary conditioned readout rank vs two distinct composed: 8 4
Conditional successful scattering equals both Hecke generators; minimum hold 1/2
Local coin probabilities: [1/9, 4/9, 4/9] ; channel eigenvalues on X,Y,Z: 2/3,0,1/3
amplitude braid defect coefficient: -4/27
Born transfer braid defect coefficient: 40/729
Normalized affine panel transform, braid solutions: [-1, 0, 1]
```

## historical_instrument.py

```sh
python3 historical_instrument.py
```

```text
spectrum counted CP transfer: {3/16 - sqrt(57)/16: 1, 3/16 + sqrt(57)/16: 1, -1/2: 1, 3/8: 1}
stationary density ratio after Perron normalization: (2*sqrt(2)/7 + sqrt(57)/7)**2
noncentral modular angular gap in dimensionless units: 2 log((sqrt(57)+2 sqrt(2))/7)
equal P(aabb)=P(bbaa)= (2337 - 259*sqrt(57))/5472
P(next a | aabb), P(next a | bbaa) = 1283/1796 - 665*sqrt(57)/16164 513/1796 + 665*sqrt(57)/16164
difference: 385/898 - 665*sqrt(57)/8082 approx -0.19248204645571004
hex tensor state: stationary spectral ratio = (4*sqrt(3)/15 + sqrt(273)/15)**2
hex geometry P(next a | aabb), P(next a | bbaa): 0.06481134741718261 0.9351886525828174
hex geometric signed difference: 11970/3019 - 149310*sqrt(273)/510211
```

## hodge_readout.py

```sh
python3 hodge_readout.py
```

```text
two given double exchanges: unique minimum E(i) = 2
with third symmetric slope: unique hexagonal minimum E = 4
eight shortest hidden transports: unique minimum sum 4e = 272
Gaussian amplitudes of eight hidden histories at square: ['1/3', '1/3', '1/3', '1/3', 'sqrt(10)/10', 'sqrt(10)/10', 'sqrt(10)/10', 'sqrt(10)/10']
one hyperbolic history, trace 6 ; max amplitude 2/trace = 1/3 ; axis (x-1)^2+y^2=2
mirror Gaussian Gram diagonal p = y/sqrt(y**2 + 1) ; off diagonal 1; det = -1/(y**2 + 1)
square and hex determinants: -1/2 -4/7
parabolic subgroup kernel: reflected aa/bb Gram [[1/2, 0], [0, 1/2]] ; values aaaa, aabb: 1/2 0
finite positive-word pairs checked: 3969
```

## kernel_audit.py

```sh
python3 kernel_audit.py
```

```text
Ordinary full-group Gram lower bounds: 1/2 1/4
Two kernels agree on the entire positive half; hidden vacuum response: 0 1/16
```

## kontrollen/modular_winding.py

```sh
python3 kontrollen/modular_winding.py
```

```text
Distinct positive braids at degrees 0..14: [1, 2, 4, 7, 12, 20, 33, 54, 88, 143, 232, 376, 609, 986, 1596]
Critical fugacity: -1/2 + sqrt(5)/2 ; Gibbs endpoint probabilities: [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
Subcritical fugacity 1/2, in e,a,b,(ab),(ba),(aba) order: [184/585, 116/585, 116/585, 64/585, 64/585, 41/585]
Squared critical norm of partial hidden a^2 b^(-2) insertion: 3/2 - sqrt(5)/2
Alternative representation (a^2 a)(b^2 a)^(-1) relative norm squared: -1/2 + sqrt(5)/2
Primitive isometry coefficient test at <ba|X*X|ab>: [[0, 1], [0, 0]]
Stationary connected range-projection variance: -2 + sqrt(5)
```

## kontrollen/path_os.py

```sh
python3 kontrollen/path_os.py
```

```text
Nonlazy n, paths to g, paths to identity: [(0, 0, 1), (1, 0, 0), (2, 0, 4), (3, 0, 0), (4, 1, 28), (5, 0, 0), (6, 20, 244), (7, 0, 0), (8, 307, 2412)]
At operational time tau, first hidden-loop coefficient in exp(tau)*p_tau(g): 1/6144
Exact time-reflected principal minors: [5/16, 835/16384, 88242633/17179869184, 6919264750817/18014398509481984]
Lazy discrete heat returns and hidden endpoints: [(2, 5/16, 0), (4, 167/1024, 1/4096), (6, 6605/65536, 65/65536)]
Unconditioned winding moments orders 2,4,6: [tau, 3*tau**2 + tau, 15*tau**3 + 15*tau**2 + tau]
Length-four strings with trivial coarse observables but nonidentity: 8 / 256 ; distinct coarse-null lift endpoints: 9
Pulled back curvature vs intrinsic conical curvature, in pi units: 8 0
Added symmetric hidden-jump menu: new derivative at g = eta/ 4 ; projected derivative identically zero.
Repulsion m, pi^2 times winding rate(l=1), relative radius, product: 1 3 5 15
Repulsion m, pi^2 times winding rate(l=1), relative radius, product: 3 1 11 11
Repulsion m, pi^2 times winding rate(l=1), relative radius, product: 5 3/5 17 51/5
```

## kontrollen/reference_clock.py

```sh
python3 kontrollen/reference_clock.py
```

```text
Central braid coordinates (w,M) = (6, -1, 0, 0, -1)
Cyclic phase order and first surviving positive central exponent: [(1, 1), (2, 1), (3, 1), (4, 2), (5, 5), (6, 1), (7, 7), (8, 4), (9, 3), (10, 5), (11, 11), (12, 2)]
vacuum moments of one-mark channel, n=0..8: [1, 1, 1, 1, 1, 1, 1, 1, 1]
primitive moments of one-mark channel, n=0..8: [1, 0, 0, 0, 0, 0, 0, 0, 0]
central nonempty moments of one-mark channel, n=0..8: [1, 0, 0, 0, 0, 0, 0, 0, 0]
hidden moments of one-mark channel, n=0..8: [1, 0, 1/4, 0, 15/128, 0, 19/256, 0, 895/16384]
exp(tau) times hidden autocorrelation through order eight: 179*tau**8/132120576 + 19*tau**6/184320 + 5*tau**4/1024 + tau**2/8 + 1
Checked endpoint pointer covariance on 2589 histories and all four elementary exchanges; active grading 3 odd, 2 even.
Haar charged correlation is exactly exp(-tau) for all nonzero windings.
```

## lorentz_control.py

```sh
python3 lorentz_control.py
```

```text
Lie-generator minor 4*b**2*(a - d)*(a + d)/(a*d - b**2)**2
q=2 raw Born effect rank 3; native imaginary direction on M2(R); Lie rank 3+3
hex full quadratic raw Born effect rank 3; native imaginary direction on M2(R); Lie rank 3+3
Hermitian cone determinant: t**2 - x**2 - y**2 - z**2
```

## markenwirkung.py

```sh
python3 markenwirkung.py
```

```text
Parity-preserving label permutations preserving multiplication: 2 / 12
Triangle relabelings preserving all flag products: 2 / 6
Squared normalized HS norms of raw odd galleries: [1/2, 1/2, 1/8]
Three-chart covariance: 6 x 3 x 6 blocks; fixed named pair braids in 1 / 3 charts
Equal raw visibility in every active class: chart weights 1/(q + 1) ; unique thick order q=2; raw counts [1, 4, 4, 4, 4, 4]
Uniform chart raw excess of each transposition over each 3-cycle: q*(q - 2)*(q - 1)/3
Permutations of reverse Weyl types giving all six successes: 1 / 720
Uniform recovered channel on the 1/6 probability grid: 1 / 462
```

## originalvergleich.py

```sh
python3 originalvergleich.py --repo "$TFPT_REPO"
```

```text
Independently imported q* word projectors vs flag-return channel: 256 / 256
```

## paired_frame.py

```sh
python3 paired_frame.py
```

```text
D8 vector roots / half-spin weights / parity of minus count: 112 128 1
D3(A3), D5 roots: 12 40
glue degree 0,1,2,3 roots: [52, 64, 60, 64]
generator norms in D3*, D5*, sum: 3/4 5/4 2
D3+D5 among frame triples / S3-invariant triples: 56 [(2, 4, 6), (3, 5, 7)]
S3 diagonal sign lifts, all and up to coordinate-sign gauge: 256 16
preserve active difference lines, all and up to pair-preserving gauge: 32 8
also swap reference lines for odd permutations, all / gauge classes: 16 4
also fix appended positive parity coordinate, all / gauge classes: 8 2
also preserve positivity of the other six point axes: (sign masks for a,b) [(128, 0)]
alternative, all 7 old point vectors stay positive; lifts / swap lifts: [(0, 0), (1, 0)] [(1, 0)]
hexagon H0 and H1 (transposition): +1, -1
```

## positive_history.py

```sh
python3 positive_history.py
```

```text
X,Y eigenspectra {1/2: 1, -1: 1} ; (XY)^3 = [[1/8, 0], [0, 1/8]]
sum of off-identity coefficient majorants <= 1
words aa bb Gram before common damping r^4 [[17/32, 5/16], [5/16, 17/32]] normalized contrast 7/17
words aabb bbaa Gram before common damping r^8 [[41/256, 17/128], [17/128, 41/256]] normalized contrast 7/41
positive histories through length 0 OS rank 1
positive histories through length 1 OS rank 3
positive histories through length 2 OS rank 4
positive histories through length 3 OS rank 4
exchange-symmetric pure boundary vector: minimal OS rank 2; mixed-history squared distance before r^8: 7/128
hex H1 / tensor-square / traceless pair overlaps: 1/4 1/16 1/4
tensor-square active span: q + 1/q = 14; contrasts aa/bb, mixed order 90/97 1080/1177
```

## signed_clifford.py

```sh
python3 signed_clifford.py
```

```text
all 16 signed lifts of one 3-cycle, fixed native quarter:
  (group, Pauli kernel, J exactly conserved) -> number [((48, 8, True), 8), ((96, 16, False), 8)]
kernel commutant gamma monomials (zero-based): [(), (0, 1)]
minimal extension invariant densities: (I+s*(I tensor Y))/4, -1 <= s <= 1
all stabilizer subgroup obstructions + attaining flag fibres: [4, 4, 2, 2, 2]
14 flags on one J sector: 8 cube vertices + 6 face centres
pointwise invariant determinant winding is divisible by 2 (order 48) or 4 (order 96)
positive angular winding one: {(0, 0, 0, 1)} minimum winding with one zero mode: 3
rank-one current coefficient times 2 pi: 1/4
Theta²=-I, Theta fixes 5/5 gamma, Theta J Theta^-1=-J: unique symmetric state I4/4
any relative orbit weights: Lüders map per sector = (rho+tr(rho)I2)/3
Kramers-even pointwise loop ranks: 0,2,4. One moving-flag loop current: 1/(8 pi)
self-reflected pair loop has rank two and current 1/(4 pi); explicit full lifted cover: 2 x 14
for 10 signed Kramers tests, Pf(u_P(z)B)/Pf(B)=z; det u_P=z²
M4 bare fastest spectrum: 1, (1/5)x5, (3/5)x10
uniform empty + five shortest tests: native 1, (1/3)x5, (2/3)x10 on 16/16 basis words
fastest law on 0.1 grid: 1 / 1001 (unique analytically on simplex)
full 16x16 native vs six-word comparison: 256/256; 12*Phi = diag(12,[4]*5,[8]*10)
binary CP instrument effect dimensions, total / fixed all-input probability: 36 / 20
same full transfer and half-probability record, conditional gamma1: [Fraction(1, 1), Fraction(-1, 3)]
```

## word_coordinates.py

```sh
python3 word_coordinates.py
```

```text
full Gram lower bound: 1/3
rank of reflected Gram for words of length <=2: 4
f(a^2) = f(b^2) = 5/512
OS squared distance of a^2 and b^2: 7/65536
f(a^3), f(b^2 a), difference: 7/8192 0 7/8192
Hecke central transfer (AB)^3: Matrix([[1/8, 0], [0, 1/8]])
reflection invariance, braid, full-group mixed-word samples and Gram bounds: OK
```

## Unveränderter ursprünglicher TFPT Prüfer

```sh
python3 "$TFPT_REPO/experiments/theory-contracts/marked-seam-cp-lift-20260921/checker.py" --repo "$TFPT_REPO" --out /tmp/marked-seam-comparison.json
```

```text
{"verdict": "PARTIAL", "count": 65, "weights": ["7/12", "1/12", "0"]}
```
