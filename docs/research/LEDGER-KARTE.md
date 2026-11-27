# Ledger-Karte: was zu welchem Teil der Theorie bereits feststeht

Erzeugt von `research/originlab/ledger_map.py --write` aus `verification/status_ledger.csv`. Nicht von Hand ändern. Die Klasse ist aus dem Statustext des Ledgers abgelesen und keine neue Bewertung; maßgeblich bleiben das Ledger und die Originalquellen.

**784 aktive Behauptungen:** 256 exakt, 410 exakt, Teil bedingt oder offen, 45 bedingt, 39 offen, 3 Kill-Test, 2 gesetzt, 29 Regel oder Synthese.

Zweck: Vor jeder Aussage „neu“ oder „offen“ hier nachsehen. Suche im Volltext mit `ledger_map.py --find 'Muster'`.

## Übersicht

| Gebiet | Blöcke | exakt | exakt, Teil bedingt oder offen | bedingt | offen | Kill-Test | gesetzt | Regel oder Synthese |
|---|---|---|---|---|---|---|---|---|
| Eingaben, Anker und Selbstkonsistenz | CAR 11, ARCH 11, ANCHOR 4, OVERDET 4, TFPT 4, DOILY 4, SYNTH 3, P2 3, AX 2, BOOT 2, META 2, CARRIER 2, SEL 1, FOUR 1, CENTER 1, FIRSTPRINCIPLES 1, P1 1, K5 1, PG32 1 | 27 | 25 | 1 | 1 | 0 | 2 | 3 |
| Naht, Rand und Transfer | SEAM 120, FTRANSFER 8, FTR 6, DYN 6, GNET 5, QFT 4, DHR 3, TRANSFER 2, BOUNDARY 1, CMNORM 1, MARKLOCAL 1, RATE 1, DECOUPLING 1, DYNAMICS 1, RESIDUAL 1, RES 1, SHEET 1 | 51 | 91 | 1 | 12 | 0 | 0 | 8 |
| E8, Gitter, Uhren und Codes | E8 42, FORM 23, ARF 5, GALOIS 4, ARITH 4, COX 3, CURVE 3, CFIN 3, COXETER 2, H2 1, CAS 1, MCKAY 1, QEC 1, DUAL 1, GOLD 1, GAUSSIAN 1, NORMALFORM 1 | 63 | 28 | 1 | 1 | 0 | 0 | 4 |
| Standardmodell, Eichstruktur und 4D | GATE 46, PS 17, QFT4D 13, CHIRAL4D 5, SM 2, DIMENSION 2, QCD 1, HIGGS 1, TFPT4D 1, GAUGE 1, PARENT 1 | 14 | 54 | 8 | 5 | 3 | 0 | 6 |
| Feinstruktur und Elektromagnetismus | ALPHA 9, EM 5 | 3 | 11 | 0 | 0 | 0 | 0 | 0 |
| Flavor, Massen und CP | FLAV 63, DIAMOND 13, CP 5, SEED 4, FLAVOR 2, LEP 1, THETA13 1 | 41 | 41 | 4 | 2 | 0 | 0 | 1 |
| Twistor und Himmelskugel | CELEST 23, WOIT 9 | 18 | 10 | 2 | 2 | 0 | 0 | 0 |
| Gravitation, Horizont und Geometrie | QGEO 76, GRAV 28, HOR 22, QGAMB 13, HYP 10, FLATAWAY 7, QG 3, VGEO 2, CAUSAL 1, SCALE 1, GRAVITY 1, BULK 1 | 17 | 128 | 6 | 11 | 0 | 0 | 3 |
| Kosmologie, Vorhersagen und Daten | FR 27, COSMO 6, PRED 3, DATA 2, DM 2, CCC 2, OBSERVATORY 1, LAMBDA 1, KILL 1, AS 1, UVBRANCH 1, METROLOGY 1, OBS 1 | 10 | 15 | 19 | 4 | 0 | 0 | 1 |
| Methodik, Pruefung und Sonstiges | REDTEAM 6, CONTRACT 4, CORRECTIONS 3, RH 2, REVIEW 2, AUDIT 2, TOPO 2, REG 1, WITNESS 1, TRANSLATE 1, EXTREV 1, MESSAGE 1 | 12 | 7 | 3 | 1 | 0 | 0 | 3 |

## Eingaben, Anker und Selbstkonsistenz

### ARCH (11)

- **ARCH.E8DEG.01** [exakt] Organizing simplification: the compiler atoms ARE the E8 fundamental Casimir invariant degrees deg(E8)={2,8,12,14,18,20,24,30} (exponents+1) [L]. CLEAN/non-trivial readings: 2=|Z2|, 8=rank … · `v66_e8_casimir_degrees.py`
- **ARCH.CFT.01** [exakt] Boundary-CFT bridge (WZW central charges = compiler atoms) + SU(3)<->SU(4) reconciliation [I]/[L]; NOT a closure. (A) level-1 WZW central charges c=dim(g)/(1+h^v) are the atoms: … · `v61_cft_bridge.py`
- **ARCH.CORE.01** [exakt] Compiler core - the whole integer skeleton from one pair (g_car,N_fam)=(5,3): rank E8=g+n=8 (D5(+)A3 rank split), |Z2|=g-n=2 (carrier-family), |mu4|=(g+n)/2=4; (N_fam,|mu4|,g_car)=(3,4,5) … · `v53_compiler_core.py`
- **ARCH.MICRO.01** [exakt] Compiler micro-lemmas (review synthesis, all exact identities compressing scattered factors): (1) SPINE QUOTIENT LADDER - the recurring O(1) factors 2/3,3/4,4/5,5/4,4/3 are adjacent … · `v74_compiler_micro_lemmas.py`
- **ARCH.REVIEW.01** [exakt] Four validated review identities (passed arithmetic AND the anti-numerology/admissible-invariant discipline) + one disciplined rejection: (1) HYPERCHARGE LUCAS-BINET: carrier charges … · `v79_review_identities.py`
- **ARCH.SPINE.01** [exakt] SPINE TETRAHEDRON - the compiler microkernel as ONE finite object (review proposal 2026-06-10, validated and typed): the spine T={2,3,4,5}={|Z2|,N_fam,|mu4|,g_car} equals … · `v91_spine_tetrahedron.py`
- **ARCH.QUAD.01** [exakt, Teil bedingt oder offen] EXTERNAL-REVIEW VALIDATION + the DEGREE-2 pattern named: (1) SEED NORMAL FORM [I] - phi0 = (|mu4|/N_fam) c3 + Omega_adm c3^|mu4| exactly: the seed is a linear seam term plus a topological … · `v106_review_validation.py`
- **ARCH.TRIAD.01** [exakt] SECOND EXTERNAL-REVIEW VALIDATION (2026-06-11b): the anchor ratio triad, the 121 audit lemma and the normal-form micro-identities, all exact: (1) MICRO-IDENTITIES [I]: Omega_adm = 2 p1 p2 = … · `v119_review_validation_2.py`
- **ARCH.RRCAR.01** [exakt, Teil bedingt oder offen] THE RIEMANN-ROCH CARRIER READING (external-review proposal 2026-06-14): (g_car, N_fam)=(5,3) are the two canonical invariants of ONE object, the four-point seam divisor D=mu4 on P^1. [E]: … · `v189_riemann_roch_carrier.py`
- **ARCH.RRCAR.02** [exakt, Teil bedingt oder offen] RIEMANN-ROCH CARRIER HARDENED FROM g_car TO D5 ITSELF (external-review point 7, 2026-06-14): v189 read g_car=5=h^0(P^1,O(mu4)) and matched rank(D5); this goes one step further -- the 5-dim … · `v197_rr_carrier_clifford_d5.py`
- **ARCH.K3.01** [exakt, Teil bedingt oder offen] KUMMER/K3 UNIFICATION OF SEAM + CARRIER-16 + E8 (2026-06-17): one smooth Kummer/K3 surface carries the seam (pillowcase = the elliptic-involution quotient of the square torus tau=i, j=1728, … · `v260_k3_kummer_unification.py`

### CAR (11)

- **CAR.PAS.01** [exakt] g_car = 5 unique Pascal solution (proof + check) · `v2_carrier_pascal.py`
- **CAR.SM.01** [exakt] N_fam=3; Omega_adm=48; b1=41/10; dim S+ = 16 · `v2_carrier_pascal.py`
- **CAR.UNIQ.01** [exakt] g_car=5 is the UNIQUE rank with N_fam in Z+ and g+N_fam=8; split (b,s)=(3,2) unique from b+s=5,b^2+s^2=13; Tr Y=Tr Y^3=0 (anomaly-free) · `v14_carrier_uniqueness.py`
- **CAR.PASCAL.01** [exakt] Red-team Target B residual REDUCED [I]/[L]: the carrier 'K=2 Pascal truncation' flagged by the red team as a free postulate is NOT free. K=(g-1)/2 is the Pascal-row MIDPOINT: … · `v83_e8net_holomorphic_uniqueness.py`
- **CAR.LADDER.01** [exakt, Teil bedingt oder offen] PASCAL LADDER THEOREM (Target B residual sharpened to zero slack): the closure 2^{g-1} = sum_{k<=K} C(g,k) has for every K exactly one solution g = 2K+1 (K=1..6 scanned over g=1..60; … · `v108_pascal_ladder.py`
- **CAR.PAIR.01** [exakt, Teil bedingt oder offen] SHEET-PAIRING LEMMA (the first exact representation-theoretic anchor of the QBL programme; connects red-team Target B to the sheet Z2): character-exact weight-multiset results from the D5 … · `v109_sheet_pairing.py`
- **CAR.PAIR.02** [exakt, Teil bedingt oder offen] CALDERON-SHEET SELECTION THEOREM (which Calderon involutions certify a scalar kernel; the v109 follow-up): on the carrier H = S+ (+) S- the certified bilinear data of an involution eps … · `v110_calderon_sheet.py`
- **CAR.QTRANS.01** [exakt, Teil bedingt oder offen] QUADRATIC-TRANSPORT THEOREM (the transport half of the QBL analytic core becomes a theorem): grade seam transport by Clifford degree in the 10 field operators on the 5-slot Fock space (S = … · `v111_quadratic_transport.py`
- **CAR.COUNT.01** [exakt, Teil bedingt oder offen] SELF-COUNTING CHANNEL (the Pascal closure is an IDENTITY of the certified channel, not an extra requirement; closes the state half of the QBL analytic core at the structural level): (1) … · `v112_selfcounting_channel.py`
- **CAR.QFREE.01** [exakt, Teil bedingt oder offen] ONE KERNEL IS THE WHOLE NET (the QBL input merges with the R2/holomorphy premise; anchors the carrier choice to the Gate-A objects): (1) MAJORANA BOOKKEEPING [I] - level-1 central charges … · `v113_quasifree_kernel.py`
- **CAR.SM.ANOM.01** [exakt, Teil bedingt oder offen] THE CARRIER HALF-SPINOR IS ONE ANOMALY-FREE SM GENERATION + THE SM RG (2026-06-21, TOE attack point 3): connects the discrete compiler to constructive 4D QFT, closing the … · `v310_carrier_sm_anomaly.py`

### ANCHOR (4)

- **ANCHOR.GEN.01** [exakt] anchor-first: a=(1,1,2) elementary symmetric e_k=(4,5,2)=(|mu4|,g_car,|Z2|), c3=1/(2 e1 pi)=1/(8pi); power sums p_n=2+2^n give |R(E8)|=p1p2p3=240, dim E8=248, rank=p4-p3=8, binary ladder … · `v23_anchor_generator.py`
- **ANCHOR.VGEO.01** [exakt, Teil bedingt oder offen] v_geo = CLOSED METROLOGY STRUCTURE (R_+ unit torsor; one dimensionful calibration necessary by dimensional analysis -- the No-Unit theorem v153/v725; honest reframe, not a derived mass and … · `v78_vgeo_floor.py`
- **ANCHOR.VGEO.02** [exakt, Teil bedingt oder offen] NO-UNIT THEOREM (v_geo re-typed from open gap to CLOSED METROLOGY STRUCTURE -- R_+ unit torsor; one dimensionful calibration necessary by dimensional analysis, v153/v725; external review … · `v153_no_unit_theorem.py`
- **ANCHOR.AFFINE.01** [exakt] ANCHOR.AFFINE.01 + FLAVOR.BIDIR.01: the anchor affine normal form and the flavor bidirectional checksum -- two exact corpus compressions, ONE module from two probes (18/18 + 15/15 checks, … · `v832_anchor_flavor_checksum.py`

### DOILY (4)

- **DOILY.ANCHOR.CLASSIFY.01** [exakt] v880 (part 2) -- the EXCLUSIVE anchor classification of PG(3,2) lines: EVERY of the 15 doily (isotropic) lines carries the local weight multiset {2 - q*} = {2,1,1} = the anchor a = (1,1,2) … · `v880_finite_anchor_theorems.py`
- **DOILY.HODGE.FACTOR.01** [exakt] v880 (part 3) -- the doily IS the multiplication table Lambda^2 E x Lambda^2 E -> Lambda^4 E of the five-slot carrier: via the v845 parity lift the three iota-supports of every doily line … · `v880_finite_anchor_theorems.py`
- **DOILY.PLUCKER.GR25.01** [exakt, Teil bedingt oder offen] THE GR(2,5) PLUCKER READING + THE S6 SPECHT CHANNELS + THE HADAMARD ANCHOR + THE VACUUM SELECTION (closed 2026-08-09, round CV; promoted round CVIII in v888; probe … · `v888_finite_pfaffian_signature.py`
- **DOILY.PFAFFIAN.WICK.01** [exakt, Teil bedingt oder offen] THE FERMIONIC PFAFFIAN SIGNATURE OF THE FINITE COMPILER (closed 2026-08-09, round CVII; promoted round CVIII in v888; probe k6_pfaffian_selfhosting_probe.py 42/42 K6PFAFFIAN-CLOSED, exact … · `v888_finite_pfaffian_signature.py`

### OVERDET (4)

- **OVERDET.FLOOR.01** [bedingt] THE UNCONDITIONAL IMPROBABILITY FLOOR + THE L1<->L2 BRIDGE (2026-06-27, v432): the defensible over-determination number that survives even if a skeptic REJECTS the declared formula grammar … · `v432_overdet_floor.py`
- **OVERDET.FLOOR.02** [exakt, Teil bedingt oder offen] HARDENING THE UNCONDITIONAL FLOOR -- THE ASSUMPTION-MINIMAL COUNTING FLOOR (2026-06-29, v436): OVERDET.FLOOR.01 (v432) multiplies the counting census by SUBJECTIVE chance assignments … · `v436_overdet_floor2.py`
- **OVERDET.WITNESS.MAP.01** [exakt, Teil bedingt oder offen] The multiply-vs-compress FRAMEWORK for over-determination (2026-06-25; classification refined by v428/OVERDET.WITNESS.RECLASS.01): a reviewer's methodological correction made … · `v427_overdet_witness_map.py`
- **OVERDET.WITNESS.RECLASS.01** [exakt, Teil bedingt oder offen] The honest correction of v427 (2026-06-25): applying v427's own multiply-vs-compress test to v427 itself. [E] THE ONE GENERATOR (v236/TOPO.BRIESKORN.01): the (2,3,5) Brieskorn singularity … · `v428_overdet_reclass.py`

### TFPT (4)

- **TFPT.SELFLOOP.01** [Regel oder Synthese] THE BIRD'S-EYE RETHINK (2026-06-22, refining TFPT.IRREDUCIBLE.01/v352). Stepping all the way back: TFPT is not a LINEAR theory (axioms -> theorems) at all -- it is a CLOSED SELF-CONSISTENT … · `v353_selfloop_capstone.py`
- **TFPT.IRREDUCIBLE.01** [Regel oder Synthese] THE FRAMEWORK REDUCTION (2026-06-22) -- the deepest honest answer to 'the inputs are not axioms, they are back-determined by the theory'. BOTH P1 (c3=1/(8pi)) and P2 (g_car=5) are … · `v352_framework_irreducible.py`
- **TFPT.POSITIVE_DESCENT.MASTER.01** [exakt] Abstract positive-descent theorem and G_net finite hypothesis checks; the prime packet instance is maintained by Prime Front. Original module ID retained; original complete mixed record is … · `v817_positive_descent_master.py`
- **TFPT.TOE.COMPLETE.01** [offen] THE MASTER TOE-COMPLETENESS CONTRACT [O] (registered 2026-08-29, review wave 4). Eight gates, all named; NOTHING claimed closed. T1 structure postulate -> P1/P2/compiler (AX.P1.01, … · `v1012_toe_gate_witnesses.py`

### P2 (3)

- **P2.PARTITION.01** [exakt, Teil bedingt oder offen] g_car = 5 AS A COROLLARY OF THE FOUR MARKS (2026-07-20; the partition sharpening of ARCH.CORE.01/v53): a multiset of exactly 3 POSITIVE INTEGERS summing to 4 is UNIQUE -- {1,1,2}, the … · `v491_p2_partition_corollary.py`
- **P2.TYPING.01** [exakt] THE P2 WEIGHT-TYPING POSTULATE HARDENED (2026-07-21; types the v491/P2.PARTITION.01 check-9 residual): the anchor weights a = (1,1,2) are the Birkhoff-Grothendieck SPLITTING EXPONENTS of … · `v499_p2_weights_deligne_bg.py`
- **P2.QUADRATIC.EXHAUSTION.01** [exakt, Teil bedingt oder offen] v880 (part 4) -- QUADRATIC ONE-STEP EXHAUSTION FORCES g_car = 5 ([E] theorem / [O] premise split): the demand Lambda^2 E ^ Lambda^2 E = Lambda^{g-1} E (one quadratic composition step … · `v880_finite_anchor_theorems.py`

### SYNTH (3)

- **SYNTH.INVENTORY.01** [exakt] THE RESIDUAL INVENTORY (one constant, one anchor, one clock to find; the COMPLETE machine-pinned gap list): (A) ONE-CONSTANT INVENTORY [I] - 2/3 = |Z2|/N_fam appears exactly in seven … · `v105_residual_inventory.py`
- **SYNTH.INVENTORY.02** [exakt] RESIDUAL INVENTORY UPDATED (post v110-v122; supersedes the CONTENTS of SYNTH.INVENTORY.01 while keeping its cardinality): (1) thirteen new ledger rows machine-pinned (CAR.PAIR.02, … · `v123_inventory_update.py`
- **SYNTH.INVENTORY.03** [exakt] TERMINAL RESIDUAL INVENTORY (post premise-(A) closure; supersedes the CONTENTS of SYNTH.INVENTORY.02): after the (A)-chain (v160-v165) premise (A) is NO LONGER an independent structural … · `v167_inventory_post_closure.py`

### AX (2)

- **AX.P1.01** [gesetzt] c3 = 1/(8 pi) boundary seam kernel -- a DECLARED input but NOT a free knob: it reduces to the parabolic anchor a=(1,1,2) via c3 = 1/(2 e1(a) pi) = 1/(8pi) (e1(a)=4=|mu4|; … · `-`
- **AX.P2.01** [gesetzt] g_car = 5 carrier interface (3+2) -- a DECLARED input but an OVER-DETERMINED bootstrap fixed point, not a free choice: forced three independent ways by the E8 closure (bootstrap … · `v1005_relative_pencil_bridges.py`

### BOOT (2)

- **BOOT.CLASS.01** [exakt] D5(+)A3 is the UNIQUE familyful cyclic-glue decomposition of E8 (16-dim half-spinor + 3-family A3); glue Z4=|mu4| · `v15_bootstrap_classification.py`
- **BOOT.01** [exakt] bootstrap mu^2-5mu+4=0; g_car=5 three ways; 8=rankE8=h(D5)=phi(30) · `v6_bootstrap.py`

### CARRIER (2)

- **CARRIER.PETERSEN.RADIAL.01** [exakt] CARRIER.PETERSEN.RADIAL.01 (+ TRANSPORT.SIXTHROOT.01): the K5 edge machine, the Petersen carrier geometry with the radial dynamics whose sixth power IS the deployed transport spectrum, + … · `v808_petersen_sixthroot.py`
- **CARRIER.CUTCODE.01** [exakt] CARRIER.CUTCODE.01 (+ CARRIER.MOMENT.INCIDENCE.01): the Arf split as the cut classification of K5 + the one-line master moment with the frozen 210 decision, ONE module from two probes … · `v810_cutcode_moment.py`

### META (2)

- **META.DISCIPLINE.01** [exakt, Teil bedingt oder offen] The discipline audit (one module, 17 checks, ~200 s; discovery probe discipline_audit_probe.py 17/17, verdict DISCIPLINE-MEASURED). The suite parses ITSELF and freezes its … · `v649_discipline_audit.py`
- **META.LOOKELSEWHERE.01** [exakt] The look-elsewhere audit -- the bingo budget of the typed [C] numeric coincidences, quantified (one module, 14 checks, ~2 s; discovery probe look_elsewhere_probe.py 14/14, verdict … · `v664_look_elsewhere.py`

### CENTER (1)

- **CENTER.NORM.01** [exakt, Teil bedingt oder offen] CENTER BUDGET NORMS (2026-06-16): the Sheet-Diamond center budget (7,11,13) is the three LOCAL NORMS of the theory. C=R+Q diag(1,0,0), row sums (7,11,13): 7=N_Z[omega](3+2 omega) … · `v230_center_budget_norms.py`

### FIRSTPRINCIPLES (1)

- **FIRSTPRINCIPLES.BOUNDARY.01** [exakt, Teil bedingt oder offen] THE HONEST FIRST-PRINCIPLES BOUNDARY, mapped precisely (2026-06-22): answers (a) can M_Planck/v_geo be computed from first principles? and (b) can F_transfer (Koide, m_p/m_e, eta_B, axion) … · `v339_firstprinciples_boundary.py`

### FOUR (1)

- **FOUR.ROUTES.01** [Regel oder Synthese] THE FOUR BLACK-HOLE-BIRTH SOLUTION ROUTES INVESTIGATED (2026-06-22): direct answer to 'which closes ALL problems and is 100% provable' -- NONE individually; they CONVERGE on the same two … · `v343_four_routes_analysis.py`

### K5 (1)

- **K5.SIXSTEP.TRANSPORT.01** [exakt] K5.SIXSTEP.TRANSPORT.01: the six-step composition and the 10-dim transport deploy -- the (2/3)^6 rate is canonical on three spatial routes, the deployed six is the clock exponent, and the … · `v814_k5_sixstep_transport.py`

### P1 (1)

- **P1.INDEX.KMS.01** [exakt, Teil bedingt oder offen] P1.INDEX.KMS.01: the closed normalization identity c3 = 1/(I_Jones * beta_angle) = 1/(4 * 2pi) = 1/(8pi), its convention audit, and the sharpened P1 contract, one probe (13/13 checks incl. … · `v813_p1_index_kms.py`

### PG32 (1)

- **PG32.FLAGS.01** [exakt] PG32.FLAGS.01: the PG(3,2) flag normal form for the 105 + 35 vacuum completion, with the incidence kill test typed as a split, ONE module from one probe (26/26 checks, verdict … · `v834_pg32_flag_completion.py`

### SEL (1)

- **SEL.THM.01** [exakt, Teil bedingt oder offen] Boundary Carrier Selection Theorem (Theorem A): admissible boundary-carrier data SELECT (D5+A3)+mu4=E8 uniquely. Boundary mu4->H1=Z^3->A3,N_fam=3; minimal 16-spinor->D5 (only n=5); … · `v47_selection_theorem.py`

## Naht, Rand und Transfer

### SEAM (120)

- **SEAM.EQUIV.CONTINUUM.02** [exakt, Teil bedingt oder offen] THE CONTINUUM REALISATION SHARPENED BY THE BOOTSTRAP (2026-06-22): the 'c=8 ambiguity' flagged open in v277/v344 (E8 vs SO(16) both have c=8) is RESOLVED, non-circularly, by the seam's … · `v351_continuum_realisation_sharpened.py`
- **SEAM.EQUIV.CONTINUUM.03** [exakt, Teil bedingt oder offen] DIRECTION A (2026-06-23): the continuum chain of SEAM.EQUIV.01 assembled theorem-by-theorem, with the residual reduced to ONE lattice-realization input tied to the seam's one-sidedness (the … · `v356_continuum_mmst_applicability.py`
- **SEAM.MMST.INCLASS.01** [exakt, Teil bedingt oder offen] Direction 7 (2026-06-23): the seam collar is IN MMST's free-lattice-fermion scaling-limit class, so the chiral scaling limit is a chiral CFT of central charge 8 pinned to (E8)_1 -- reducing … · `v366_mmst_seam_collar.py`
- **SEAM.S3.LATTICE.01** [exakt, Teil bedingt oder offen] Track 1 (2026-06-23): the EXPLICIT lattice realisation of the S3 input -- a concrete gapped chiral-Majorana (p+ip) lattice model whose edge is the 16-Majorana SO(16)_1 carrier, … · `v367_seam_s3_lattice.py`
- **SEAM.S3.INFLOW.01** [exakt, Teil bedingt oder offen] Track 1 (2026-06-23): the bulk-edge / anomaly-inflow leg (S4 of the v356 chain) made CONCRETE on the v367 lattice model -- a finite STRIP of the gapped p+ip model has a REALISED, … · `v368_seam_s3_inflow.py`
- **SEAM.S3.CENTRALCHARGE.01** [exakt] S3 closure stack (2026-06-23): the central charge c=8 extracted NUMERICALLY from the lattice via the Calabrese-Cardy finite-size entanglement-entropy scaling of the critical free-fermion … · `v376_seam_s3_centralcharge.py`
- **SEAM.S3.E8CHARACTER.01** [exakt] S3 closure stack (2026-06-23): the chiral c=8 content of the seam is the SPECIFIC net (E8)_1, not just 'some c=8' -- shown by the exact affine character q-series and the mu4/GSO promotion … · `v377_seam_s3_e8character.py`
- **SEAM.S3.MODULAR.01** [exakt] S3 closure stack (2026-06-23): the GENUS-1 discriminator -- the torus ground-state degeneracy of the seam net is 1 (holomorphic = (E8)_1), via the KLM index condensation and the modular … · `v378_seam_s3_modular.py`
- **SEAM.S3.RP.01** [exakt] S3 closure stack / C7 (2026-06-23): reflection positivity (the Osterwalder-Schrader axiom) of the EXPLICIT gapped collar measure, verified numerically -- the remaining constructive-QFT … · `v379_seam_s3_rp.py`
- **SEAM.EQUIV.BOOTSTRAP.01** [Regel oder Synthese] THE HONEST CORRECTION OF v349's FRAMING (prompted by the right objection 'the inputs are not axioms -- they are back-determined by the theory', 2026-06-22). v349 concluded 'the raw seam … · `v350_bootstrap_inputs_correction.py`
- **SEAM.EQUIV.GOLDEN.01** [exakt, Teil bedingt oder offen] THE HONEST TEST OF 'DOES THE RAW SEAM CARRY THE GOLDEN RATIO?' (the v348 question, 2026-06-22) -- answer, computed and honest: NO. The bare carrier (D5, A3) and the dynamical seam data are … · `v349_raw_seam_golden_test.py`
- **SEAM.EQUIV.RIGID.01** [exakt, Teil bedingt oder offen] ROUTE B (THE RIGIDITY ROUTE) FOR THE ONE OPEN ARROW L2, CARRIED OUT (2026-06-22) -- it lands on a strikingly SIMPLE statement: the whole keystone reduces to ONE question, does the raw seam … · `v348_seam_rigidity_route.py`
- **SEAM.EQUIV.CLOSURE.01** [Regel oder Synthese] THE CLOSURE-MODE CLASSIFICATION OF THE ONE OPEN ARROW L2 (the seam-orbifold realisation, v346) -- a direct honest answer to 'is there ANY other way to FULLY solve it?' (2026-06-22). [E] THE … · `v347_seam_closure_modes.py`
- **SEAM.EQUIV.GEOM.01** [Regel oder Synthese] THE GEOMETRIC BRIDGE AS AN EXPLICIT SIX-LINK CHAIN, WITH THE ONE OPEN ARROW PINNED (2026-06-22) -- the honest capstone of the SEAM.EQUIV.01 investigation (v335/v336/v344/v345). The chain … · `v346_seam_geometric_bridge.py`
- **SEAM.DETK.02** [exakt, Teil bedingt oder offen] EXECUTING THE COMBINATORIAL CORE OF ROUTE R3 (the hypergraph-homotopy route to det K=1, v344) -- the most promising fresh angle, now carried out as far as it computably goes (2026-06-22). … · `v345_hypergraph_homotopy.py`
- **SEAM.DETK.01** [Regel oder Synthese] THE BIRD'S-EYE SYNTHESIS OF THE ONE OPEN KEYSTONE BIT det K=1 (2026-06-22): a full-spectrum map of SEAM.EQUIV.01's only open residual -- it is ONE statement with six equivalent faces; it is … · `v344_detk_synthesis.py`
- **SEAM.AREACOEFF.02** [exakt] Central-theorem RESIDUAL RESOLVED (honestly) as an anchor, not a gap [I-symbolic + L]. The Seeley-DeWitt coefficient for the carrier Dirac operator is a2 = (1/(4pi)^2)(1/6)(R d + 6E d) with … · `v68_seeley_dewitt_residual.py`
- **SEAM.AREACOEFF.01** [exakt] Central theorem (1/(8pi) area-law coefficient) - STRUCTURAL CLOSURE + reduction to ONE sharp coefficient (honest partial solution of Alessandro's central target). CLOSED [I]/[L]: by the … · `v67_area_law_coefficient.py`
- **SEAM.AREALAW.01** [exakt] Area-law evidence for the open Seam-Horizon step (necessary condition met, NOT a closure): (1) NUMERICAL [Python-only] - the simplest reflection-positive Gaussian system (harmonic chain) … · `v59_area_law_evidence.py`
- **SEAM.THEOREM.01** [exakt, Teil bedingt oder offen] Seam-Horizon chain (precise framing) + the OPEN Seam-Horizon Theorem. PRECISE: the seam is NOT identical to an event horizon - it is the abstract boundary normaliser whose LOCAL … · `v58_seam_horizon_chain.py`
- **SEAM.CROSSLINK.01** [exakt, Teil bedingt oder offen] Horizon cross-links - published black-hole results meeting compiler atoms (exact arithmetic [I]; physical seam=horizon identification is [P]): (1) JACOBSON/EINSTEIN: c3=1/(8pi) is exactly … · `v57_horizon_crosslinks.py`
- **SEAM.ATTRACTOR.01** [exakt] Unique attractor: a gapped boundary transport selects ONE fixed point [I]/[L]. The boundary transport spectrum {1,(2/3)^6,(1/3)^6} has spectral gap Delta=-log(2/3)^6=6 log(3/2)=2.4328>0. By … · `v56_unique_attractor.py`
- **SEAM.CYCLE.01** [exakt] A cyclic element INSIDE the compiler hull (computed, not asserted): built from the E8 Cartan matrix, the Coxeter element c=s1...s8 has order EXACTLY h(E8)=30=|Z2|*N_fam*g_car=2*3*5 … · `v55_coxeter_cycle.py`
- **SEAM.HORIZON.01** [exakt] Seam=horizon [I] keystones (pure identities, NO cosmological narrative): (1) the integer 8 in c3=1/(8pi) is OVERDETERMINED and gravitationally aligned: geometry 8=2|mu4| (8pi=|mu4|*2pi seam … · `v54_seam_horizon_keystones.py`
- **SEAM.KHALF.01** [exakt] k=c3/2 FORCED - the central seam-area coefficient (Alessandro central target met at the dimensionless level): k=c3/2=(1/2)x(1/(|Z2|*2pi*chi(S2))) splits into TWO cutoff-independent pieces, … · `v73_k_c3_half.py`
- **SEAM.AREACOEFF.03** [exakt] Seam-Horizon programme: the Fursaev-Solodukhin factor 4pi is now DERIVED, not imported [I] (removes the literature import of v67): (1) smoothed-cone Gauss-Bonnet - regularising the 2D cone … · `v90_conical_defect_chain.py`
- **SEAM.EHMODEL.01** [exakt, Teil bedingt oder offen] THE EH FORM FROM THE REPLICA VARIATION, GAPPED-MODEL LEVEL (the first movement on R3/SEAM.THEOREM.01; the gate narrows, does NOT close): (1) CONICAL DEFICIT EXACT [I] - on C/Z_N the … · `v150_replica_eh_model.py`
- **SEAM.EHMODEL.02** [exakt, Teil bedingt oder offen] THE BFK SPLIT (the Calderon-transfer step of SEAM.EHMODEL.01/v150 ANSWERED; SEAM.THEOREM.01 narrows again, stays open): (1) DOUBLING => C_N = C_D [I] - reflecting a wedge across its edges … · `v151_bfk_split.py`
- **SEAM.EHMODEL.03** [exakt] THE R3 NORMALISATION IS THE DIMENSIONFUL ANCHOR (not a separate open item; SEAM.THEOREM.01's residual list shrinks to the kernel premise alone): (1) THE COEFFICIENT IS A LOG-RATIO [I] - the … · `v152_norm_is_anchor.py`
- **SEAM.EQUIV.01** [exakt, Teil bedingt oder offen] THE SEAM EQUIVALENCE CONTRACT (2026-06-18): concentrates the whole open structure on ONE named arrow (does NOT prove it). THE NAMED THEOREM: the raw reflection-positive seam state (gap, … · `v286_seam_equivalence_contract.py`
- **SEAM.EQUIV.A01** [exakt, Teil bedingt oder offen] ROUTE A (AQFT) DEPENDENCY CHECKER (2026-06-18): for the open Gral arrow 'free Gaussian RP bulk => holomorphic single-sector boundary net', types the 5-lemma AQFT chain and isolates the ONE … · `v287_free_rp_bulk_to_holomorphic_boundary.py`
- **SEAM.EQUIV.B01** [exakt, Teil bedingt oder offen] ROUTE B (DtN) FULL-L^2 SUB-PRINCIPAL Z4 BLOCK-DIAGONALITY (2026-06-18): proves (numerically, on the full Toeplitz operator) the advisor's 'probably provable' step, lifting v201/v284 from … · `v288_full_l2_subprincipal_z4.py`
- **SEAM.ADVERSARY.01** [exakt, Teil bedingt oder offen] ROUTE-B RED-TEAM: Z4 BLOCK-DIAGONALITY != MARK-LOCALITY (2026-06-18). Constructs the strongest adversary against Route B -- a SMOOTH Z4-symmetric off-mark curvature f(theta)=eps*cos(4 … · `v290_z4_smooth_curvature_adversary.py`
- **SEAM.EQUIV.A.INV.01** [Regel oder Synthese] DISCHARGE ROUTE A'S 'INVERTIBLE BULK' VIA THE FREE-FERMION CLASSIFICATION (2026-06-19): the single open hypothesis of Route A (v297 LIT-A, 'the raw quasi-free seam bulk is invertible/SRE'); … · `v301_route_a_invertible.py`
- **SEAM.EQUIV.GAP.01** [exakt, Teil bedingt oder offen] THE LAST LEVER -- THE SEAM GAP IS THE DERIVED RECOVERY GAP (2026-06-19): after v300 (Route B = Route A's rationality) and v301 (Route A's invertibility discharged), the WHOLE open content … · `v302_seam_gap.py`
- **SEAM.EQUIV.A.LIT.01** [exakt, Teil bedingt oder offen] ROUTE A LITERATURE THEOREM STACK (2026-06-18): writes Route A's single import (v287 L3) as a precise, citable composition of three established theorems, isolating the one open hypothesis. … · `v297_route_a_literature_stack.py`
- **SEAM.EQUIV.CHAIN.01** [exakt, Teil bedingt oder offen] SEAM.EQUIV.01 COMPOSITION CERTIFICATE (2026-06-21, TOE attack point 1): assembles the v297 literature stack into ONE well-typed logical chain gap>0 -> invertible/SRE (Kitaev) -> holomorphic … · `v308_seam_equiv_chain.py`
- **SEAM.EQUIV.UNIFY.01** [Regel oder Synthese] THE KEYSTONE UNIFICATION (2026-06-22): there is ONE open theorem (SEAM.EQUIV.01), QGEO.SYM.01 is its COROLLARY (a conformal-net axiom), and QG.AMB.01 is a decoupled general QG problem, NOT … · `v335_seam_equiv_unify.py`
- **SEAM.EQUIV.CONTINUUM.01** [offen] THE ONE OPEN LEMMA OF SEAM.EQUIV.01, stated precisely and mapped to the recent rigorous literature (2026-06-22): the continuum-limit and OS-reconstruction legs are now CITABLE, the … · `v336_continuum_limit.py`
- **SEAM.S3.SCALINGLIMIT.01** [exakt] An HONEST attack on the one open lemma of SEAM.EQUIV.01 -- the continuum scaling-limit existence (2026-06-24). It does NOT close v336 (the abstract existence is the cited MMST theorem); it … · `v392_seam_s3_scalinglimit.py`
- **SEAM.RIGIDITY.01** [exakt, Teil bedingt oder offen] THE SEAM STATE RIGIDITY THEOREM AS ONE NAMED TARGET (2026-06-24, Paper A / the G_net closure contract): consolidates the scattered state-invariance reductions … · `v398_seam_state_rigidity.py`
- **SEAM.EQUIV.02** [exakt, Teil bedingt oder offen] The dual keystone at tau=omega (2026-06-24): the family/flavor sector as the order-3 Eisenstein/A2 face of E8, dual to the seam=(E8)_1 at tau=i (SEAM.EQUIV.01/v404). A NAMED CONTRACT (genre … · `v405_seam_equiv_omega.py`
- **SEAM.GALOIS.01** [exakt, Teil bedingt oder offen] The seam mu4 IS the carrier's Galois group (2026-06-25): a positive resolution of the v409/RES.COXETER.SYMMETRY.01 cyclic/Galois asymmetry. In the cyclotomic picture the flavor clocks … · `v419_seam_galois_carrier.py`
- **SEAM.GALOIS.NET.01** [exakt, Teil bedingt oder offen] The Galois<->Net bridge (2026-06-25): the seam mu4 = Gal(Q(zeta5)) is the SAME cyclic Z/4 as the (E8)_1 simple-current glue -- NOT a mere order-4 coincidence -- connecting the Galois … · `v422_galois_net_bridge.py`
- **SEAM.EQUIV.LTORP.01** [exakt, Teil bedingt oder offen] The (LTO-RP) reduction of the keystone (2026-06-25): an HONEST MMST-style reduction of SEAM.EQUIV.01's one open lemma -- the raw collar's intrinsic Bisognano-Wichmann property (the residual … · `v424_seam_lto_rp_reduction.py`
- **SEAM.EQUIV.LTORP.KMS.01** [exakt, Teil bedingt oder offen] The invertible beta=1-KMS variant of (LTO-RP) (2026-06-25): the one corner Naaijkens-Penneys-Wallick [NPW26, arXiv:2605.10693] leave open -- they prove (LTO-RP) only for TOPOLOGICALLY … · `v426_seam_lto_rp_kms.py`
- **SEAM.EQUIV.BW.HSMI.01** [exakt, Teil bedingt oder offen] The half-sided modular inclusion (Borchers/Wiesbrock) as an INTRINSIC Bisognano-Wichmann route (2026-06-30) for face (ii) of SEAM.EQUIV.01 -- geometric modular flow WITHOUT presupposing … · `v438_seam_hsmi_borchers.py`
- **SEAM.EQUIV.CONTINUUM.LR.01** [exakt, Teil bedingt oder offen] The Lieb-Robinson/Nachtergaele-Sims continuum-existence MOTOR for the gapped seam collar (2026-06-30) -- grounds the BULK half of the cited MMST existence (SEAM.MMST.INCLASS.01/v336) in the … · `v439_seam_lieb_robinson.py`
- **SEAM.EQUIV.LTORP.BETA.01** [exakt, Teil bedingt oder offen] The beta-interpolation connecting the NPW26-PROVEN trace case (beta=0) to the seam's beta=1 KMS case along an obstruction-free homotopy (2026-06-30) -- upgrading the single-point v426 … · `v440_seam_lto_rp_beta.py`
- **SEAM.EQUIV.APPLICABILITY.01** [exakt, Teil bedingt oder offen] The applicability ledger (2026-06-30) -- a hypothesis-by-hypothesis audit of the THREE external theorems the seam keystone cites (MMST scaling limit, NPW26 (LTO-RP), … · `v441_seam_applicability_ledger.py`
- **SEAM.RIGIDITY.UNIFORM.01** [exakt, Teil bedingt oder offen] The uniform-in-N hardening of the Seam State Rigidity Theorem (2026-06-30; hardens SEAM.RIGIDITY.01/v398) -- the band-limited mu4-character rigidity is NOT a small-N artifact. v398 verified … · `v442_seam_rigidity_uniform.py`
- **SEAM.E8.DISCRIMINATOR.01** [exakt] The empirical/structural DISCRIMINATOR bundle (2026-06-30; mode D of the v347 closure-mode classification) -- concrete FALSIFIABLE kill tests separating the claimed seam=(E8)_1 from the … · `v443_seam_e8_discriminator.py`
- **SEAM.EQUIV.CONTINUUM.EDGE.01** [exakt, Teil bedingt oder offen] The CHIRAL EDGE scaling-limit MOTOR (2026-06-30) -- the edge half of the one SEAM.EQUIV.01 continuum residual, attacking exactly the step v439 left open. v439 grounded the BULK … · `v444_seam_edge_scaling.py`
- **SEAM.RIGIDITY.FORCING.01** [exakt, Teil bedingt oder offen] The rigidity FORCING theorem (2026-06-30; upgrades v442/SEAM.RIGIDITY.UNIFORM.01 from 'permit' to 'force') -- v442 (and v398) established that a mu4-block-diagonal K commutes with the … · `v445_seam_rigidity_forcing.py`
- **SEAM.RIGIDITY.CLOCKINV.01** [exakt, Teil bedingt oder offen] DERIVING clock-invariance ('transfer in the commutant') from RP-symmetry (2026-06-30; closes the single residual of v445/SEAM.RIGIDITY.FORCING.01 one notch further) -- v445 reduced rigidity … · `v446_seam_clock_invariance.py`
- **SEAM.EQUIV.EDGE.CHERN.01** [exakt, Teil bedingt oder offen] An INDEPENDENT (topological) reading of the edge central charge c_-=8 (2026-06-30) -- corroborating v444/SEAM.EQUIV.CONTINUUM.EDGE.01 from a different observable. v444 read c_- from the … · `v447_seam_edge_chern.py`
- **SEAM.EQUIV.EDGE.FOURPOINT.01** [exakt, Teil bedingt oder offen] Pinning the ONE external MMST fact empirically (2026-06-30) -- the edge FOUR-point function is quasi-free and conformally convergent, extending v444's two-point scaling limit to the next … · `v448_seam_edge_fourpoint.py`
- **SEAM.EQUIV.MMST.UNIFORM.01** [exakt, Teil bedingt oder offen] UNIFORM-IN-N control of the edge scaling limit (2026-06-30) -- the empirical content of the cited MMST uniform-convergence theorem (the single residual of SEAM.EQUIV.01 after … · `v449_seam_mmst_uniform.py`
- **SEAM.EQUIV.EDGE.ENTANGLEMENT.01** [exakt, Teil bedingt oder offen] A THIRD independent reading of the edge central charge c_-=8 (2026-06-30) -- the Calabrese-Cardy entanglement scaling, after v444 (correlator) and v447 (topology). For a critical 1D … · `v450_seam_edge_entanglement.py`
- **SEAM.EQUIV.EDGE.CARDY.01** [exakt, Teil bedingt oder offen] The Cardy finite-size conformal TOWER (2026-06-30) -- the edge spectrum reproduces the EXACT Ising primary dimensions {1/16,1/2}, naming the edge CFT. v450 fixed the central charge c=1/2; … · `v451_seam_edge_cardy_tower.py`
- **SEAM.EQUIV.E8.MODULAR.01** [exakt, Teil bedingt oder offen] The (E8)_1 torus MODULAR DATA (2026-06-30) -- one primary, S-invariant, T-phase e^{-2pi i/3}, the holomorphic c=8 signature, extending v156's character q-expansion to the modular … · `v452_seam_e8_modular.py`
- **SEAM.RIGIDITY.MU4FROMMARKS.01** [exakt, Teil bedingt oder offen] DERIVING the mu4-symmetry of the seam RP data (QGEO.SYM.01) from the four marks (2026-06-30) -- the LAST structural premise of the rigidity chain folded into the existing realisation … · `v453_seam_mu4_from_marks.py`
- **SEAM.EQUIV.EDGE.VIRASORO.01** [exakt, Teil bedingt oder offen] The lattice level-1 current algebra with Sugawara c=8 (2026-06-30, G2 of the post-F next steps) -- testing the cited MMST 'Koo-Saleur lattice generators -> continuum Virasoro with the … · `v454_seam_edge_virasoro.py`
- **SEAM.EQUIV.BW.UNIFORM.01** [exakt, Teil bedingt oder offen] A UNIFORM-IN-N Tomita-Takesaki tower for the intrinsic Bisognano-Wichmann condition u_Theta=J (2026-06-30, G1 of the post-F next steps) -- advancing v424 sub-step (i) (the continuum lift), … · `v455_seam_bw_uniform.py`
- **SEAM.S3.FROM-P1.01** [exakt, Teil bedingt oder offen] The edge chirality (c_- != 0, the S3 topological phase) is FORCED by the one-sidedness that defines c3=1/(8pi) (P1) (2026-06-30, G6 of the post-F next steps) -- so S3 is a CONSEQUENCE of … · `v456_seam_chirality_from_c3.py`
- **SEAM.EQUIV.KILLTEST.01** [exakt] An explicit (E8)_1-vs-SO(16)_1 character/sector kill-test for the seam edge (2026-06-30, G7 of the post-F next steps) -- a FALSIFIABLE consistency test of the (E8)_1 identification of … · `v457_seam_character_killtest.py`
- **SEAM.EQUIV.MMST.AUDIT.01** [offen] The EXACT hypothesis-by-hypothesis citation audit of the MMST scaling-limit theorem against the collar (2026-06-30, G3 of the post-F next steps, the 'exakter Zitat-Abgleich') -- 'does the … · `v458_seam_mmst_citation_audit.py`
- **SEAM.EQUIV.LATTICEVOA.01** [exakt, Teil bedingt oder offen] The SECOND, independent route to (E8)_1 -- the lattice-VOA / conformal-net A_Q(E8) construction supplies EXACTLY the 128-spinor extension MMST (v458) leaves open (2026-06-30, G5 of the … · `v459_seam_lattice_voa_route.py`
- **SEAM.S3.LTO.01** [exakt, Teil bedingt oder offen] The EXPLICIT flat-band commuting-projector (local-topological-order) realisation of the S3 input (2026-06-30) -- advancing v424 sub-step (i) ('realise the raw seam as a Z/2-reflection … · `v460_seam_s3_lto.py`
- **SEAM.S3.LOCALITY.01** [exakt, Teil bedingt oder offen] The Kapustin-Fidkowski/Thouless strict-locality OBSTRUCTION made EXPLICIT (2026-06-30) -- sharpening v460's verdict from a cited statement to an exhibited topological obstruction on the … · `v461_seam_strict_locality.py`
- **SEAM.EQUIV.SPINOR.01** [exakt, Teil bedingt oder offen] The 128-spinor extension SO(16)_1 -> (E8)_1 exhibited at CHARACTER level and as a finite-L -> continuum convergence (2026-06-30) -- the constructive companion of the v458 (MMST audit) / … · `v462_seam_spinor_continuum.py`
- **SEAM.EQUIV.UNIQ.01** [exakt, Teil bedingt oder offen] The c=8 HOLOMORPHIC-UNIQUENESS pin (2026-06-30, B of the A+B 'actually close it' round) -- the positive lever making the IDENTIFICATION half of SEAM.EQUIV.01 classification-forced rather … · `v463_seam_c8_holomorphic_uniqueness.py`
- **SEAM.EQUIV.ONEPARTICLE.01** [exakt, Teil bedingt oder offen] The one-particle RIGIDITY + scaling limit attacking R1 (2026-06-30, A of the A+B round) -- the REALISATION input of SEAM.EQUIV.01. The seam is quasi-free (v113/v160/v161), so by Araki … · `v464_seam_oneparticle_rigidity.py`
- **SEAM.EQUIV.CROSSEDPRODUCT.01** [exakt, Teil bedingt oder offen] The NET-LEVEL CROSSED-PRODUCT certification of the 128-spinor extension + the R1 -> R1' invariant-level reduction (2026-07-02, closure-route round) -- re-founds BOTH halves of the … · `v469_seam_crossedproduct_route.py`
- **SEAM.EHMODEL.04** [exakt, Teil bedingt oder offen] SEAM-HORIZON REPLICA CHAIN AT COLLAR LEVEL (the kernel-identification premise of SEAM.THEOREM.01, exercised numerically; SEAM.EHMODEL.03's residual list sharpens again, the gate stays … · `v471_seam_horizon_replica.py`
- **SEAM.KRONHEIMER.01** [exakt, Teil bedingt oder offen] THE R3 COMBINATORIAL-TO-GEOMETRIC BRIDGE IS KRONHEIMER'S THEOREM (2026-07-14, bird's-eye round): SEAM.DETK.01/v344 named R3's residual 'the rewrite attractor (a graph) IS the du Val … · `v479_kronheimer_quiver_bridge.py`
- **SEAM.CONTACT.UNIT.01** [exakt, Teil bedingt oder offen] THE SEAM CONTACT UNIT -- THE MERGED RESIDUAL OF THE TWO ANALYTIC TARGETS (2026-07-14, deep-dive round II): the shared 'c3 per boundary insertion' rule of ALPHA.QUILLEN.EXACT.01 (c3-ladder … · `v484_seam_contact_unit.py`
- **SEAM.CONTACT.UNIT.02** [exakt, Teil bedingt oder offen] THE DIAGONAL CHANNEL SETTLED -- THE MERGED TARGET COLLAPSES ONTO THE KEYSTONE (2026-07-14, deep-dive round III; completes the computable side of SEAM.CONTACT.UNIT.01/v484): [E] DIAGONAL … · `v485_contact_diagonal_closed.py`
- **SEAM.CCC.01** [exakt, Teil bedingt oder offen] THE CHIRAL CENTRAL CHARGE c_- = 8 MEASURED FROM THE BULK GROUND STATE (2026-07-20; ground-state witness round): the Kim-Shehab-Kim modular commutator J(A,B,C) = i Tr(rho_ABC [K_AB, K_BC]) = … · `v489_seam_modular_commutator.py`
- **SEAM.MU.01** [exakt, Teil bedingt oder offen] THE SPIN-STRUCTURE PARITY CENSUS: mu_pre = 4 WITNESSED ON THE LATTICE (2026-07-20; ground-state witness round, companion of SEAM.CCC.01/v489): the T^2 ground-state fermion parity of the … · `v490_seam_parity_census.py`
- **SEAM.CLOCK.RIGIDITY.01** [exakt] SEAM CLOCK RIGIDITY -- ROUTE B OF THE v503 CLASSIFICATION EXECUTED (2026-07-21; two exact theorems in one module, Part A Moebius + Part B Fock, 34 checks; NO marker moves -- QGEO.REALIZE.01 … · `v506_seam_clock_rigidity.py`
- **SEAM.BIT.ORIGIN.01** [exakt] THE TAUTOLOGY ATTACK ON THE v506 ALIGNMENT BIT, EXECUTED AND REFUTED (2026-07-21; one module, Part A torus origin + Part B fermion dichotomy, 26 checks; NO marker moves -- the bit stays … · `v507_seam_bit_origin.py`
- **SEAM.BIT.FREEDOM.01** [exakt] THE FREEDOM ATTACK ON THE v506/v507 ALIGNMENT BIT -- THE POSITION HALF OF THE BIT IS TOPOLOGY (2026-07-21; one module, Part A dihedral/fermionic + Part B Moebius/RP-geometry, 30 checks; NO … · `v510_seam_bit_freedom.py`
- **SEAM.TAU.FLAG.01** [exakt] THE FLAG-TRANSITIVITY EQUIVALENCE WEB -- THE tau = i ATTACK EXECUTED ON THE FREE RP SEAM CIRCLE (2026-07-22; one module, 29 checks; NO marker moves -- the bit stays [C], QGEO.REALIZE.01 and … · `v512_seam_tau_flag.py`
- **SEAM.EQUIV.MMST.01** [bedingt] THE MMST ROUTE OF THE SEAM EQUIVALENCE TARGET (2026-07-22 route split, carved out of SEAM.EQUIV.01 -- an ID restructuring, NO status change in substance; this row carries exactly the … · `v366_mmst_seam_collar.py`
- **SEAM.EQUIV.TWISTOR.01** [offen] THE TWISTOR ROUTE OF THE SEAM EQUIVALENCE TARGET (2026-07-22 route split, the second continuum route named in its own right -- an ID restructuring, NO status change in substance): the open … · `v492_celestial_z4_orbifold.py`
- **SEAM.BIT.RPBLIND.01** [exakt] THE RP/THETA ATTACK ON THE ALIGNMENT BIT DECIDED (2026-07-22; one module, 25 checks, ~25 s; an honest, decided NEGATIVE result in the v508 style): does free reflection positivity (v519 RP) … · `v521_seam_bit_rp_blind.py`
- **SEAM.BIT.TWISTBLIND.01** [exakt] THE TWIST-STATE ATTACK ON THE ALIGNMENT BIT DECIDED (2026-07-23; one module, 20 checks, ~80 s; an honest, decided NEGATIVE result in the v508 style -- the NINTH side-blind test on the v512 … · `v525_seam_bit_twist_blind.py`
- **SEAM.THERMAL.KMS.01** [exakt, Teil bedingt oder offen] THE THERMAL SEAM: KMS -> NARIAI/BH BRIDGE (2026-07-23; one module, 22 checks, ~8.4 s; verdict SUCCESS): the reconstructed free OS quotient of the 16-Majorana NS seam circle is KMS, and its … · `v526_seam_thermal_kms_nariai_bridge.py`
- **SEAM.CLOCK.SILVER.01** [exakt] THE SILVER CLOCK EIGENVALUE EXPLAINED AND DEMYSTIFIED (2026-07-23; one module, 22 checks, ~10 s; double verdict, anti-numerology style): STRUCTURE = SUCCESS (closed-form origin of {1, … · `v527_seam_clock_silver_spectrum.py`
- **SEAM.BIT.TWISTCLASS.01** [exakt, Teil bedingt oder offen] THE ALIGNMENT BIT AS AN EXACT MEASURABLE TWIST-CLASS CHOICE (2026-07-23; one module, 18 checks, ~46 s; verdict SUCCESS -- the physical definition): the alignment bit is the choice of which … · `v528_seam_bit_twist_class_definition.py`
- **SEAM.INT.FKTOY.01** [exakt, Teil bedingt oder offen] THE MINIMAL INTERACTING FIDKOWSKI-KITAEV SEAM TOY (2026-07-23; one module, 24 checks, ~97 s; triple-verdicted -- the FIRST firing Kill-Test-2 shadow). F1 = KILL, F2 = KILL, F3 = ERFOLG. (1) … · `v529_seam_interacting_toy_fk.py`
- **SEAM.STRADDLE.CONE.01** [exakt, Teil bedingt oder offen] THE RP CONE ON THE INTERACTING FK QUARTETS -- THE FIRST DYNAMICAL SELECTION OF THE ALIGNMENT BIT (2026-07-24; one module, 15 checks, ~7 s; double-verdicted; the named v534 candidate of the … · `v534_seam_straddle_cone.py`
- **SEAM.CLIFFORD.MODULAR_S.01** [exakt] the missing Clifford index-2 bit IS the metaplectic lift of the modular S transformation at its fixed point tau = i -- coset decision positive, deck placed, anchor named; the frozen … · `v798_seam_clifford_modular_s.py`
- **SEAM.CODE.TYPEII.01** [exakt] Type II (doubly-even self-dual) is FORCED from three physical seam axioms -- the direct attack on residual R1 of ARF.BOUNDARY.CODE.01, decided exhaustively (19/19 checks, ~5 s; discovery … · `v799_seam_code_typeii.py`
- **SEAM.CFIN.WICKFUNCTOR.01** [offen] THE WICK FUNCTOR ARC [O] (constructed/obstructed/resolved 2026-08-09, rounds CVIII-CX; promoted round CXI in v896): (1) THE SIX ROLES EXIST in the deployed seam system (D5 slot pairs + A3 … · `v896_wick_block_functor.py`
- **SEAM.CFIN.KMSMIX.01** [exakt, Teil bedingt oder offen] THE CHANNEL-MIXING CANDIDATE STATE (constructed 2026-08-09, round CXIII; promoted round CXV in v898): the v896 demand -- a channel-mixing state realizing the block Wick functor's value … · `v898_kms_schur_mixing.py`
- **SEAM.STATE.DERIVATION.01** [exakt] THE RP/MODULAR CARVING OF THE SEAM STATE SPACE (measured 2026-08-10, round CXXIV; promoted in v903): on the v898 family the reflections and the modular route CARVE for the first time. (1) … · `v903_seam_rp_exclusion.py`
- **SEAM.CFIN.MIXING.NORMALIZATION.01** [exakt] DEAD READING (negative result, first-class; measured 2026-08-10, round CXVI; promoted in v903): 'THE MIXING STRENGTH t = 1/8 = 2pi c3/|Z2| IS A COMPILER VALUE' IS A DECORATION, NOT A … · `v903_seam_rp_exclusion.py`
- **SEAM.CFIN.MINIMAL.MEDIATOR.01** [exakt] DEAD READING (negative result, first-class; measured 2026-08-10, round CXVII; promoted in v903): 'N_fam = 3 IS THE MINIMAL BOUNDARY-MEDIATION RANK' IS REFUTED at theorem grade for this … · `v903_seam_rp_exclusion.py`
- **SEAM.CFIN.TWISTED.RP.01** [exakt] THE COMPLETE TWISTED-OS CENSUS IS EXCLUSIONARY (measured 2026-08-10, round CXXXI; promoted in v903): the COMPLETE census of C6-twisted reflections Theta_g = U_g o theta (2 bases x 6 C6 … · `v903_seam_rp_exclusion.py`
- **SEAM.CFIN.RP.DILATION.01** [exakt] THE DILATION ROUTE SPLITS (measured 2026-08-10, round CXXXII; promoted in v903): on the explicit C6-covariant parent family over the v898 split 16 = 10 (+) 6 the dilation mechanism is … · `v903_seam_rp_exclusion.py`
- **SEAM.CFIN.GAP.PENCIL.01** [exakt] t_gap IS THEOREM-SHAPED (measured 2026-08-10, round CXXXIII; promoted in v903): t_gap = 0.2309488708... is the UNIQUE positive root of the irreducible cubic q(t) = 9t^3 + 21t^2 - t - 1 … · `v903_seam_rp_exclusion.py`
- **SEAM.CFIN.MARGINAL.STABILITY.01** [exakt] STRICT 2-CYCLE RP IS IMPOSSIBLE ON THE WHOLE COVARIANT CLASS (measured 2026-08-10, round CXXXVIII(60); promoted in v908): the closure-vs-isolated dilemma of the v903 marginal witness is … · `v908_seam_equilibrium_wiring.py`
- **SEAM.CFIN.NESS.PARENT.01** [exakt] NO NESS NEEDED -- THE PRICE OF THE MIXING IS EXACTLY ZERO ENTROPY PRODUCTION (demand closed as NOT-NEEDED; measured 2026-08-10, round CXXXIX(60); promoted in v908): strict theta_S (collar) … · `v908_seam_equilibrium_wiring.py`
- **SEAM.STATE.WIRING.SELECTOR.01** [exakt] THE WIRING FREEDOM THEOREM -- IS PURE-I COMPILER-FORCED? ANSWERED: NO -- DEPLOYMENT-CHOICE, NOT COMPILER THEOREM (the v908-registered contract, closed 2026-08-12, rounds 69-71 notes … · `v911_wiring_freedom.py`
- **SEAM.BULK4D.RECON.01** [offen] THE SEAM-TO-BULK 4D RECONSTRUCTION CONTRACT [O] (registered 2026-08-27): the theorem previously only implicit between the SEAM.EQUIV.01 conclusion (a chiral c=8 boundary object, 1+1D) and … · `-`
- **SEAM.INT.ODDSPLIT.01** [exakt] THE ODD-SPLIT LAW (measured 2026-08-27; promoted in v972): the v529 straddle law refined on the 16-Majorana FK seam toy via the Jaffe-Pedrocchi Majorana-RP class made DECIDABLE -- the … · `v972_seam_interaction_front.py`
- **SEAM.INT.OSAVOID.01** [exakt] OS VIABILITY ON THE QUARTET-AVOIDING FAMILY (measured 2026-08-27; promoted in v972): the v529 kill-test-2 counterweight -- (1) the straddle census 24/24 reproduced; (2) the clock-invariant … · `v972_seam_interaction_front.py`
- **SEAM.STATE.RPMIXING.01** [exakt, Teil bedingt oder offen] THE RP MIXING CENSUS (measured 2026-08-27; promoted in v972): the v903 exclusivity is RAY- and REFLECTION-relative, NOT slice-wide. The full 33-dim C6-covariant mixing slice censused under … · `v972_seam_interaction_front.py`
- **SEAM.EQUIV.SKELETON.01** [exakt] THE CROSS-ROUTE OS SKELETON (2026-08-27, v973): the two continuum routes of the seam equivalence target -- the MMST scaling-limit route (the p+ip collar edge) and the twistor route's Z4 … · `v973_seam_route_narrowing.py`
- **SEAM.SIMPLECURRENT.GENERATOR.01** [offen] THE SIMPLE-CURRENT GENERATOR CONTRACT [O] (registered 2026-08-27 from the external bird's-eye review; the seam-gap compression conjecture). CONTEXT: v973 carved the cited-theorem coverage … · `v983_simple_current_generator.py`
- **SEAM.DETLINE.UNIFICATION.01** [offen] THE DETERMINANT-LINE UNIFICATION HYPOTHESIS [O] (registered 2026-08-27 from the external bird's-eye review; a precise research hypothesis, NOT a result). OBSERVATION: three separately open … · `v991_detline_bulk_edge.py`
- **SEAM.MMST.TYPEIII.CHARGED.01** [offen] THE MASTER SCALING-LIMIT THEOREM FOR THE CHARGED/TYPE-III IDENTIFICATION [O] (registered 2026-08-29, review wave 5). DEMAND: ONE scaling-limit theorem Psi_lambda,N -> Psi_lambda strongly on … · `v1006_mmst_lemma_battery.py`
- **SEAM.MILNOR.LOCALRING.01** [offen] THE MILNOR LOCAL-RING IDENTIFICATION CONTRACT [O] (registered 2026-08-30, review wave 6, v1004). REQUIRED: a CANONICAL identification between the raw seam collar / E8 transition bus and the … · `v1004_milnor_frobenius_gray.py`
- **SEAM.MMST.TELB.C1C2A.01** [exakt] TEL-B C1/C2a SUBCLAIM CONTRACT (registered 2026-09-05, v1025; no parent marker move). C1: a native 100-bit Arb/Acb certificate integrates 3888 analytic rectangles, adds a closed axis-strip … · `v1025_telb_c1_c2a_certificates.py`
- **SEAM.MMST.TELB.ROUND3.01** [exakt] TEL-B ROUND-3 MIXED CERTIFICATE (registered 2026-09-05, v1022; no parent marker move). Native 256-bit Arb/Acb over 5812 rational cells gives TV_p(F'') <= 15.817695962988924... < 15.928674 … · `v1022_telb_round3_certificates.py`
- **SEAM.MMST.TELB.NORM.01** [exakt] For the fixed M=1 Ny=8 collar, CF<0.065697, DG<0.042959, C2b<0.466902 and one-sided CROSS<0.282637/N imply Bound-B squared<0.714386. Validated 512-node Acb Fourier core plus the complete … · `v1026_telb_round4_closure.py`
- **SEAM.MMST.ALG2.FRAME.01** [exakt, Teil bedingt oder offen] Given full neutral block algebras, nonzero charged corners generate the full finite matrix algebra iff their sector graph is connected. A noncollapsing actual-word frame with uniform word … · `v1030_alg2_sector_frame.py`
- **SEAM.MMST.DISORDER.REGISTER.01** [exakt, Teil bedingt oder offen] Adding an auxiliary local Z4 flux register constructs order-four charged corners in the common QWZ matter Fock space. Scalar matter gauge transformations alone cannot change seam holonomy. … · `v1033_charged_disorder.py`

### FTRANSFER (8)

- **FTRANSFER.SUITE.01** [exakt, Teil bedingt oder offen] THE F_transfer FUNCTOR CERTIFICATE (2026-06-24, Paper E / the physical-interface contract): the four frontier sectors are ONE functor F_transfer = F_observable o F_threshold o F_RG, … · `v402_ftransfer_suite.py`
- **FTRANSFER.STT.01** [exakt] PHYS moonshot slice (T3 / F_transfer) -- can thermal time replace the external clock contract? The conjugacy-class KILL test (one module, 11 checks, ~1 s; discovery probe … · `v723_phys_modular_clock.py`
- **FTRANSFER.COCYCLE.01** [exakt] T3-B -- conditioned modular flows from ONE internal mother state, the second thermal-time KILL (one module, 12 mechanical checks, ~1 s; discovery probe phys_t3b_conditioned_modular_flows.py … · `v724_phys_t3b_modular_flows.py`
- **FTRANSFER.CLOCKS.01** [exakt] Preregistered theory-contract search target (frozen 2026-08-04, experiments/ftransfer-clocks/hypotheses/ftransfer_clocks_v1.yaml; EXECUTED 2026-08-05, v777): the external-clock contract for … · `v777_ftransfer_clock_jets.py`
- **FTRANSFER.KMSCOSET.01** [exakt] the coset assignment of the fibered clock functor is a THEOREM from KMS/modular typing (18/18 checks incl. the three must-fire controls, ~3 s; discovery probe ftransfer_kms_coset_probe.py, … · `v792_ftransfer_kms_coset.py`
- **FTRANSFER.GENERATING.01** [offen] THE GENERATING FUNCTIONAL CONTRACT [O] (registered 2026-08-27 from the external bird's-eye review; the constructive successor named by the FR.TRANSFER.01 guard). The four … · `v990_wj_transduction_mechanism.py`
- **FTRANSFER.SK.RHO0.01** [offen] THE SCHWINGER-KELDYSH INITIAL-STATE CONTRACT [O] (registered 2026-08-28, review wave 3; mechanism executed 2026-08-29). DEMAND: cosmological transfers (F_relic, F_Boltzmann, and any … · `v1000_spin2_sk_witnesses.py`
- **FTRANSFER.PARENT.Z5CAP.01** [exakt, Teil bedingt oder offen] FINITE PARENT-INTERNAL Z5 CAP AND ENLARGEMENT BOUNDARIES (registered 2026-09-05, v1024; no parent marker move). A 12-dimensional internal cap has gap 0.197423715893 and preserves the finite … · `v1024_parent_internal_boundaries.py`

### DYN (6)

- **DYN.SEMIGROUP.01** [exakt, Teil bedingt oder offen] STATIC->DYNAMIC FLIP (2026-06-17): the seam's modular flow (reversible) and its recovery channel read as a Markov/Lindblad generator (irreversible), with the dissipative gap EQUAL to the OS … · `v238_modular_lindblad_dynamics.py`
- **DYN.KMS.01** [exakt, Teil bedingt oder offen] KMS / THERMAL TIME (2026-06-17): the modular flow of v238 makes the seam state a KMS (thermal) state at the canonical inverse temperature beta = 1, and the geometric … · `v239_kms_thermal_time.py`
- **DYN.GNS.01** [exakt, Teil bedingt oder offen] GNS / OS RECONSTRUCTION (2026-06-17, finite toy): from the discrete seam data (a *-algebra M + the state omega) one reconstructs a Hilbert space, a cyclic vector, the modular flow … · `v240_gns_os_reconstruction.py`
- **DYN.TRANSFER.UNIVERSAL.01** [exakt, Teil bedingt oder offen] The single-flow reduction of F_transfer (2026-06-25): the four frontier transfers (F_pole, F_Boltzmann, F_relic, F_QCD) do NOT each need a separate external solver -- their DYNAMICS is the … · `v425_dyn_transfer_universal.py`
- **DYN.MARKOV.EMBED.01** [exakt, Teil bedingt oder offen] THE MARKOV EMBEDDING CONTRACT (registered 2026-08-27; 'three meanings of dynamics' split, step (a)->(b)). FINITE KERNEL HALF [E] (v971, sympy-exact): for the v221 seam transfer kernel T … · `v971_markov_embedding_generator.py`
- **DYN.UNITARY.DILATION.01** [offen] THE UNITARY DILATION CONTRACT [O] (registered 2026-08-27; 'three meanings of dynamics' split, step (b)->(c)). The continuous semigroup e^{tQ} (DYN.MARKOV.EMBED.01) is irreversible … · `v984_markov_qca_dilation.py`

### FTR (6)

- **FTR.PATH.01** [exakt, Teil bedingt oder offen] DIAMOND F_TRANSFER PATH (2026-06-16): F_transfer is ONE path on the Sheet Diamond, not four external interfaces. Family M(s,t)=R+Q diag(s,t,t): the transfer corners K=M(1,-1), C=M(1,0), … · `v224_diamond_ftransfer_path.py`
- **FTR.DISC7.NORM.01** [exakt] THE SHEET TRANSFER IS A NORM LINE OF DISCRIMINANT -7 (2026-07-23; 12 checks, <1 s, sympy exact): the branch-preserving F_transfer path J -> K -> C -> F = M(1, -2..1) (FTR.PATH.01/v224, … · `v533_ftransfer_disc7_norm.py`
- **FTR.CRDISC.01** [exakt, Teil bedingt oder offen] The preregistered cross-ratio check (self-code review Priority 4 / v566 S9), executed on its DISCRETE half (one module, 9 checks, <1 s; discovery probe cr_discrete_corners_probe.py 9/9, … · `v571_cr_discrete_corners.py`
- **FTR.CRCONT.01** [exakt, Teil bedingt oder offen] The CONTINUOUS half of the preregistered cross-ratio check, executed against the in-repo solver dynamics (one module, 7 checks, <1 s; discovery probe cr_continuous_uniflow_probe.py 7/7, … · `v575_cr_continuous_uniflow.py`
- **FTR.JETS.01** [exakt, Teil bedingt oder offen] The jet/Schwarzian contract, executed natively (one module, 6 checks, <1 s; discovery probe native_jets_probe.py 6/6, verdict MOEBIUS-IN-OWN-CLOCK; the second review's Priority 6): the … · `v578_native_jets.py`
- **FTR.PGL2.01** [exakt, Teil bedingt oder offen] The F_transfer PGL2 round (one module, 25 checks, ~1 s; discovery probe ftransfer_pgl2_probe.py 25/25, verdict COMMON-ACTION-DISCRETE / UNIFLOW-CONTINUOUS -- NOT KILLED). The v533 … · `v632_ftransfer_pgl2.py`

### GNET (5)

- **GNET.RAMIFIED.01** [exakt] PHYS moonshot slice (T2 / G_net) -- is the Ramond projection (1+i)-adic? (one module, 9 checks, ~1 s; discovery probe phys_ramified_ns_r_probe.py 9/9, verdict RAMIFIED-SEES-NSR). THE … · `v722_phys_ramified_ns_r.py`
- **GNET.PPINDEX.01** [exakt] PHYS T2 main slice -- Pimsner-Popa / Watatani index 4 on the CAR ladder, with the mu4 action DERIVED from the clock and the Ramond zero-mode test (one module, 57 checks, ~7 s; discovery … · `v726_phys_car_pp_index.py`
- **GNET.LOCALNET.01** [exakt, Teil bedingt oder offen] The G_net local net functor -- net axioms exact, Watatani index locally exactly 4, Ramond as half-line/bond defect (one module, 16 checks, ~2 s; discovery probe … · `v746_phys_gnet_local_functor.py`
- **GNET.GNSLIMIT.01** [exakt] The G_net continuum strand -- finite-level precursors of the GNS limit net COHERENT, the exact C[Z4] obstruction dim law, and the honest Arf negative; one module from two probes (27 + 27 … · `v779_gnet_gns_arf.py`
- **GNET.MARTINGALE.LIMIT.01** [Regel oder Synthese] the noncommutative-martingale route to the G_net GNS limit -- its finite hypotheses MEASURED, all met (15/15 checks, controls 4/4, ~7 s; discovery probe gnet_martingale_probe.py, … · `v802_gnet_martingale.py`

### QFT (4)

- **QFT.OSMOMENT.01** [exakt, Teil bedingt oder offen] ATOMIC OS MOMENT LEMMA + SUGAWARA GAP SAFETY (the admissible transfer sector's OS positivity as a positive moment problem, and its gap as an E8 Sugawara inequality): (1) … · `v171_os_moment_cluster.py`
- **QFT.FOCK.01** [exakt] SEAM FOCK-SPACE READINGS (three exact audit bridges that sharpen the QFT/topological reading of already-used objects): (1) CAR CONE COMPRESSION [I] -- the even CAR algebra of the 16 seam … · `v174_seam_fock_readings.py`
- **QFT.NETRP.01** [exakt, Teil bedingt oder offen] NET EXISTENCE + FULL-CONE RP DISCHARGED TO [F] (the heavy artillery; the 'three open QFT theorems' are two [F] + one irreducible [O], nothing fabricated): (1) FULL-CONE REFLECTION … · `v175_net_existence_full_cone.py`
- **QFT.MSC.01** [exakt, Teil bedingt oder offen] MODULAR SPECTRAL CLOSURE (2026-06-17): the assembly/reduction certificate turning the whole QFT layer (v238-v260) into ONE relative object … · `v261_modular_spectral_closure.py`

### DHR (3)

- **DHR.SECTORS.01** [exakt, Teil bedingt oder offen] PARTICLES = DHR SUPERSELECTION SECTORS (2026-06-17): the discriminant/anyon content of the carrier and the condensation tower (v235/v237) made explicit as an abelian modular tensor category … · `v241_dhr_sectors.py`
- **DHR.MODULAR.01** [exakt, Teil bedingt oder offen] SPIN-STATISTICS + MODULAR (S,T) DATA (2026-06-17): the exchange statistics, modular S/T matrices, Verlinde fusion and Gauss-Milgram central charge of the carrier MTC -- the … · `v242_spin_statistics_modular.py`
- **DHR.SMATRIX.01** [exakt, Teil bedingt oder offen] SCATTERING: HAAG-RUELLE SKELETON + BRAIDING S-MATRIX (2026-06-17): the S-matrix layer of the emergent QFT. The seam net has a mass gap (v64/v238), so Haag-Ruelle scattering theory … · `v243_haag_ruelle_braiding.py`

### TRANSFER (2)

- **TRANSFER.HIDDEN.CIRCULATION.01** [exakt, Teil bedingt oder offen] THE HIDDEN CIRCULATION CONTRACT, executed half (registered 2026-08-27; external-review round, companion of the DYN.MARKOV.EMBED.01/DYN.UNITARY.DILATION.01 split). [E] (v976, … · `v976_seam_lift_birkhoff.py`
- **TRANSFER.COHERENT.WILSON.01** [exakt, Teil bedingt oder offen] THE COHERENT WILSON CONTRACT, executed half (registered 2026-08-27; external-review round, companion of TRANSFER.HIDDEN.CIRCULATION.01). [E] (v977, sympy-exact): (1) UNISTOCHASTICITY: the … · `v977_transfer_unistochastic_wilson.py`

### BOUNDARY (1)

- **BOUNDARY.01** [Regel oder Synthese] the two gates are the discrete<->continuous HANDOFF, not missing formulas: quark cross-ratio cleanness tracks generation (t/b=3/26 clean err 3e-5; c/s fits 21/29 AND 34/47; u/d fits 8/17 … · `v35_quark_qcd_boundary.py`

### CMNORM (1)

- **CMNORM.DUAL.01** [exakt] CM-NORM DUALITY (2026-06-16): the two exceptional CM moduli give two central TFPT integers as ideal NORMS. SQUARE seam (Gaussian Z[i], j=1728): N(g_car+i|mu4|)=5^2+4^2=41=10 b1 -- the EM … · `v222_cm_norm_duality.py`

### DECOUPLING (1)

- **DECOUPLING.THEOREM.01** [exakt] THE DECOUPLING THEOREM (consolidated, citable, 2026-06-22): every TFPT dimensionless readout factors through the gapped admissible sector and is INDEPENDENT of the ambient QG measure -- the … · `v337_decoupling_theorem.py`

### DYNAMICS (1)

- **DYNAMICS.UNIVERSAL.01** [exakt, Teil bedingt oder offen] The UNIVERSAL spectral-gap / Perron-Frobenius principle (2026-06-23, bird's-eye synthesis): every TFPT sector is the SAME structural object -- a gapped operator with a UNIQUE leading … · `v383_dynamics_universal.py`

### MARKLOCAL (1)

- **MARKLOCAL.RAW.01** [exakt, Teil bedingt oder offen] RAW MARK-LOCALITY DECOMPOSITION (2026-06-18): decompose the ONE remaining Route-B residual after v288 ('why is the raw seam sub-principal term mark-local?') into an explicit 5-lemma chain, … · `v289_marklocal_raw.py`

### RATE (1)

- **RATE.TRANSLATE.01** [exakt, Teil bedingt oder offen] DO THE DISCRETE-KERNEL AND DYNAMIC RATES TRANSLATE? AN EXACT NUMBER-FIELD SPLIT (2026-06-21, v312/v313 follow-up). [E] FIELD SPLIT: the affine-E8 adjacency spectrum is the rational part … · `v314_rate_translation.py`

### RES (1)

- **RES.COXETER.SYMMETRY.01** [exakt, Teil bedingt oder offen] CLOSED AS A STRUCTURAL LEMMA (2026-06-24, v409): prime-2 (the seam sheet involution) is NECESSARY but NOT autonomous, so no prime-2-only attractor exists -- resolving the v394 open … · `v394_coxeter_atoms.py`

### RESIDUAL (1)

- **RESIDUAL.CERTIFICATION.01** [exakt, Teil bedingt oder offen] Historical compiler-residual audit (2026-06-23; scope corrected 2026-09-05): v384 classifies the then-listed compiler residual into three external kinds, (A) external math proof, (B) … · `v384_residual_certification.py`

### SHEET (1)

- **SHEET.GEN.BINARY.01** [exakt, Teil bedingt oder offen] The sheet generator V is a binary internal compiler (2026-06-24): the rank-2 sheet axis V=Q diag(0,1,1) of the centered cross (v95/v218) has powers that print the binary carrier spine. [E] … · `v410_sheet_generator_binary.py`

## E8, Gitter, Uhren und Codes

### E8 (42)

- **E8.GLU.01** [exakt] E8 = (D5 (+) A3) + mu4 glue: explicit lattice certificate (240 roots norm^2=2; even; unimodular Gram det=1; glue index 4=|mu4|); disc=Z4; q(D5)+q(A3)=5/4+3/4=2 · `v1_e8_glue.py`
- **E8.NUM.01** [exakt] 240 = 16*5*3 and 248 = 240+8 as carrier traces · `v1_e8_glue.py`
- **E8.UNMAPPED.BANDWIDTH.01** [exakt, Teil bedingt oder offen] A DISCIPLINED BANDWIDTH SEARCH OF THE UNMAPPED E8 REGION (2026-06-22, the next step after v354). The natural follow-up -- 'search that region for hits and structures we overlooked' -- is … · `v355_e8_unmapped_bandwidth.py`
- **E8.REVERSE.AUDIT.01** [exakt, Teil bedingt oder offen] THE REVERSE NUMEROLOGY AUDIT + WHAT THE GOLDEN RATIO MEANS (2026-06-22). TFPT's forward discipline (v305) is one-directional: every physical READOUT must map onto an E8 structure … · `v354_e8_reverse_audit.py`
- **E8.OTHERSIDE.AUDIT.01** [exakt, Teil bedingt oder offen] THE 'OTHER SIDE OF THE SEAM' AUDITED AGAINST THE UNMAPPED E8 STRUCTURE (2026-06-27, v430): the sheet/deck complement of the reverse audit v354/v355. The reverse audit found FIVE Casimir … · `v430_other_side_reverse_audit.py`
- **E8.DEGREE.LADDER.01** [exakt, Teil bedingt oder offen] THE REVERSE-AUDIT 'UNMAPPED' E8 CASIMIR DEGREES ARE THE FORCED TWO-FAMILY DECOMPOSITION (2026-06-27, v431): the structural complement of v354/v355/v430. The reverse audit found FIVE … · `v431_e8_degree_ladder.py`
- **E8.DEGREE.JOINT.01** [exakt] THE TWO STRUCTURAL INTEGERS ARE JOINTLY FORCED BY E8 DEGREE-INVARIANTS, AND E8 IS FIXED BY ITS DEGREES (2026-06-29, v437): a CONSOLIDATION of v6/v66/v355/v431 deriving more from the E8 … · `v437_e8_degree_joint.py`
- **E8.GRADE.01** [exakt, Teil bedingt oder offen] THE GRADED HULL + THE ZERO-MODE MULTIPLICITY (both named computations gain their structural anchors): R2 HALF - (1) EXPLICIT COSET CONSTRUCTION [I]: the 240 E8 roots, built as explicit … · `v128_graded_hull.py`
- **E8.SLICE.01** [exakt] E8 SLICE COMPRESSION LEMMA (the seven E8 audit slices are ONE projection of two small alphabets; anti-numerology compression, NOT new physics): (1) the anchor power sums p_n(a)=2+2^n give … · `v170_diamond_slice_compression.py`
- **E8.CHAN.01** [exakt] DEGREE/EXPONENT CHANNEL SPLIT (2026-06-16): 248=120+128 as a magnitude/phase channel typing (replaces the over-strong 'S^- is dark matter'). E8 exponents {1,7,11,13,17,19,23,29} … · `v227_degree_exponent_channel_split.py`
- **E8.NETWORK.01** [bedingt] E8 AS A NETWORK ATTRACTOR / WOLFRAM-PROGRAM READING (2026-06-18): an HONEST experimental audit of whether TFPT's structure is naturally a hypergraph/rewrite system, built on the McKay … · `v298_e8_network_attractor.py`
- **E8.CMDUAL.01** [exakt, Teil bedingt oder offen] The (E8)_1 vacuum character at the two CM points (2026-06-24): chi_E8 = E4/eta^8 = j^{1/3} equals 12 at the seam CM point tau=i and 0 at the family CM point tau=omega, extending v282 … · `v404_e8_char_cm_duality.py`
- **E8.DEGREE.MODULAR.01** [exakt] THE E8 DUAL-DEGREE MODULAR CHECKSUM 744 + THE PRIME-TOTATIVE FINGERPRINT (2026-07-23; AUDIT-TYPED under the no-free-pattern rule -- feeds no derivation, kills no alternative alone, … · `v532_e8_degree_modular_checksum.py`
- **E8.CODE.01** [exakt, Teil bedingt oder offen] The E8 code round (one module, 7 checks, ~2 s; discovery probe e8_code_probe.py 7/7, verdict E8-IS-A-CODE). Construction A on the self-dual extended Hamming code [8,4,4] yields the even … · `v626_e8_code.py`
- **E8.INCIDENCE.01** [exakt, Teil bedingt oder offen] The root incidence census (one module, 6 checks, ~1 s; discovery probe root_incidence_probe.py 6/6, verdict INCIDENCE-CENSUS-DECIDED). The naive 48 x 5 incidence to the 240 E8 roots is … · `v629_root_incidence.py`
- **E8.ZIQUOTIENT.01** [exakt] The mu4-quotient round (one module, 17 checks, ~3 s; discovery probe orbit60_quotient_probe.py 17/17, verdicts ZI-E8-QUOTIENT-CONFIRMED + CASCADE-KILLED). The quotient of the 240 E8 roots … · `v633_orbit60_quotient.py`
- **E8.ST31.01** [exakt] The ST31 structure round (one module, 56 checks, ~10 s; discovery probe st31_structure_probe.py 56/56). G31 = <the 60 unitary reflections> has order 46080 (BFS + independent Schreier-Sims) … · `v634_st31_structure.py`
- **E8.CODESEM.01** [exakt] The code-semantics rounds, kill + stage 2 (one module, 27 checks, ~3 s; discovery probes code_semantics_kill_probe.py 13/13 CODE-SEMANTICS-RELABELLED and code_semantics_stage2_probe.py … · `v638_code_semantics.py`
- **E8.RM13REV.01** [exakt] The RM(1,3) reverse reading -- are compiler numbers CODE operations? (one module, 20 checks, ~1 s; discovery probe rm13_reverse_probe.py 20/20; exact integer arithmetic). An honestly … · `v646_rm13_reverse.py`
- **E8.ST31DEG.01** [exakt, Teil bedingt oder offen] The ST31 degree-24/-20 round (one module, 25 checks, ~4 s; discovery probe st31_degree24_probe.py 25/25; exact Fraction/Z[i] arithmetic). THE REGULAR 24-CLOCK [E, structural]: an explicit w … · `v647_st31_degree24.py`
- **E8.ST31D8.01** [exakt] The ST31 degree-8 round + the unified d/4-theorem (one module, 30 checks, ~10 s; discovery probe st31_degree8_probe.py 30/30; exact Fraction/Z[i] arithmetic). ORDER-8 LANDSCAPE [E]: 10560 … · `v654_st31_degree8.py`
- **E8.GAUSSCODE.01** [exakt] The Gaussian code bridge (review contract E8.GAUSSIAN.CODE.01) -- E8(Z[i])/(1+i) IS the information space of RM(1,3): the mu4 complex structure J makes E8 a rank-4 Z[i]-lattice L; reduction … · `v689_gaussian_code_bridge.py`
- **E8.QUARTICHALF.01** [exakt] The quartic invariant restriction (review contract E8.QUARTIC.HALF.01) -- is G31 the holomorphic mu4 shadow of the E8 Coxeter compiler? ALIVE, with the honest split (one module, 22 checks, … · `v690_quartic_half.py`
- **E8.GAUSSIAN.CODE.01** [Regel oder Synthese] Review contract (external review, theorem candidate, registered 2026-08-03): E8(Z[i])/(1+i) is the information space of RM(1,3) -- the mu4 complex structure J makes E8 a rank-4 … · `-`
- **E8.QUARTIC.HALF.01** [Regel oder Synthese] Review contract (external review, theorem candidate, registered 2026-08-03): restricting the W(E8) invariant root power sums to the holomorphic mu4 eigenspace V^{1,0} kills exactly the half … · `-`
- **E8.ORBIT60.01** [exakt] The ORBIT60 census -- the 60 mu4-lines of E8 as a combinatorial object: quotient graph, sigma/clock/Coxeter dynamics on the ORBIT SPACE, and a STRICT park test for the cascade 60 -> 6 (one … · `v700_orbit60_census.py`
- **E8.ORBITPACKET.01** [exakt] The Gaussian ORBIT PACKET lemma (one module, 26 checks, ~1 s; discovery probe orbit_packet_probe.py 26/26, verdict ORBIT-PACKET-EXACT). BURNSIDE EXACT: the sigma orbit census on the 15 … · `v736_orbit_packet.py`
- **E8.SEAM48.01** [exakt] The seam-48 intertwiner test -- the sharpest structural kill of the review morning (one module, 24 machinery checks, ~4 s; discovery probe seam48_intertwiner_probe.py 24/24 machinery, … · `v737_seam48_intertwiner.py`
- **E8.PROJHAMMING.01** [exakt] The projective Hamming incidence lemma on V = L/(1+i)L -- the symplectic F2^4 of the unimodular hermitian E8 lattice, and the v737 resolution (one module, 36 checks, ~5 s; discovery probe … · `v752_projective_hamming_incidence.py`
- **E8.ONEOBJECT.01** [exakt] The d = 4 forcing census + the 60-line clock synchronization kill, one module from two probes (21/21 + 23/23 checks, ~5 s; discovery probes one_object_reduction_probe.py ONE-OBJECT-PARTIAL … · `v781_one_object_clock_census.py`
- **E8.TRANSITIONBUS.01** [exakt] The transition-bus structure theorem (35/35 checks, ~2 s; discovery probe e8_transition_bus_probe.py TRANSITION-BUS-TORSOR, 2026-08-05). The 240 roots of the unimodular hermitian Gaussian … · `v782_e8_transition_bus.py`
- **E8.CLIFFORD2Q.01** [exakt] The two-qubit stabilizer identification of the compiler carrier (32/32 checks, ~4 s; discovery probe two_qubit_clifford_probe.py CLIFFORD-PARTIAL(H4), 2026-08-05). H1 EXACT: under the … · `v783_two_qubit_clifford.py`
- **E8.TORSOR.FOURIER.01** [exakt] the projective Fourier modes of the 15 Gaussian root torsors -- carrier-block purity at CHARACTER level with the mutually-unbiased honest cap (30/30 checks, ~60 s; discovery probe … · `v800_e8_torsor_fourier.py`
- **E8.RAMIFIED.JETCODE.01** [exakt] E8.RAMIFIED.JETCODE.01 (+ E8.AFFINE.NSR.01): the jet-code normal form of the Gaussian E8 compiler + the NS/R bit as an affine 8-point geometry, ONE module from two probes (18/18 + 15/15 … · `v803_e8_jetcode_affine_nsr.py`
- **E8.SYNDROME.ALGEBRA.01** [exakt] E8.SYNDROME.ALGEBRA.01: the honest kill -- E8 is NOT an error-correction hull in the sharpest sense, one probe (13/13 checks incl. the 3 must-fire controls, ~6 s; discovery probe … · `v805_e8_syndrome_algebra.py`
- **E8.SIMPLEX.FOURIER.01** [exakt] v857 -- E8.SIMPLEX.FOURIER.01 + FLAVOR.WINDING.QUADRATIC.01: the E8 audits, part I -- the Gaussian census IS a Fourier law and v817's per-prime measurement is the prime restriction of a … · `v857_simplex_fourier_winding.py`
- **E8.G31.CLOCK_ALPHABET.01** [exakt] v858 -- E8.G31.CLOCK_ALPHABET.01: the E8 audits, part II -- the compiler clock alphabet is a NORMAL FORM of G31: Deg(G31)/4 = {2, 3, 5, 6} = {|Z2|, N_fam, g_car, |R+(A3)|} with the audit … · `v858_g31_clock_alphabet.py`
- **E8.WINDING.DECODER.01** [exakt] THE WINDING DECODER OF THE CLOCK ALPHABET (2026-08-08 round 36, v874; discovery probe g31_winding_decoder_probe.py 11/11, verdict WINDING-DECODER-EXACT; exact arithmetic). THEOREM: the v857 … · `v874_feedback_decoder_lean.py`
- **E8.SHELL.GLOBAL.01** [Regel oder Synthese] THE GLOBAL GAUSSIAN SHELL THEOREM (2026-08-08 round 36, v875; discovery probe gaussian_shell_global_probe.py 15/15, verdict SHELL-GLOBAL-THEOREM) -- the census->theorem upgrade: the v857 … · `v875_shell_global_lean.py`
- **E8.DIRECTED.READOUT.01** [exakt, Teil bedingt oder offen] THE SEVEN EXACT E8 READOUT CELLS [E] Identity (round r609, promoted 2026-09-02; module v1018, re-derived, no probe imports). (1) Seifert/Euler: S+S^T=A_E8, charpoly(-S^{-1}S^T)=Phi_30, … · `v1018_e8_directed_readout.py`
- **E8.COXETER.EULER.COMPLETION.01** [exakt] THE COXETER-EULER COMPLETION OF E8 [E] Identity + Numerical (round r617, promoted 2026-09-02; module v1019, re-derived, no probe imports). Exact: C=-S^{-1}S^T, det(I-xC)=Phi_30, Tr C=-1; … · `v1019_coxeter_euler_completion.py`
- **E8.STANDARDFORM.CHART.01** [offen] THE STANDARD-FORM CHART CONTRACT [O] (registered 2026-11-20; exploration candidate, source experiments/exploration/2026-11-20-kleinstes-verbundsystem, scripts b1-b4; no registered module -- … · `-`

### FORM (23)

- **FORM.P2.02** [exakt] P2 carrier algebra (hypercharge/family) machine-verified · `lean4-carrier-rigidity`
- **FORM.SEAM.MMST.01** [exakt] LEAN 4 FORMALISATION of the S3 / SEAM.EQUIV.01 continuum leg (the MMST scaling-limit theorem applied to the seam collar), 2026-06-23 -- new module TfptCarrier/SeamScalingLimit.lean, lake … · `lean4-carrier-rigidity/TfptCarrier/SeamScalingLimit.lean`
- **FORM.CAR.01** [exakt] Theorem A (carrier rank uniqueness g_car=5) machine-formalised in Lean 4: the Pascal carrier condition 2^(g-1)=C(g,0)+C(g,1)+C(g,2), division-free 2^g=g^2+g+2, has the UNIQUE solution g=5 … · `lean4-carrier-rigidity`
- **FORM.GLUE.01** [exakt] Lean 4 formalisation of the v89/v92 arithmetic cores - AUDIT: PASS (machine-formalised [F]): new module TfptCarrier/GlueUniqueness.lean proves by kernel decide (NO native_decide, no new … · `lean4-carrier-rigidity/TfptCarrier/GlueUniqueness.lean`
- **FORM.LADDER.01** [exakt] LEAN FORMALISATION of the anchor ladder and the Pascal ladder (hardening round, 2026-06-11): new modules TfptCarrier/AnchorLadder.lean and TfptCarrier/PascalLadder.lean, built and audited … · `lean4-carrier-rigidity (AnchorLadder.lean, PascalLadder.lean)`
- **FORM.QGEO.01** [exakt, Teil bedingt oder offen] LEAN 4 FORMALISATION of the QGEO.SYM.01 CONDITIONAL THEOREM (2026-06-15; the algebraic core of v201/v210, machine-proved) -- AUDIT: PASS, no sorry/admit, no domain axioms (#print axioms = … · `lean4-carrier-rigidity`
- **FORM.QGEO.02** [exakt, Teil bedingt oder offen] LEAN 4 FORMALISATION of QGEO.UNIFORM.01 (the geometric normal form of v177; the UNIFORM node of the QGEO proof tree) -- AUDIT: PASS, no sorry/admit, no domain axioms (only propext, … · `lean4-carrier-rigidity`
- **FORM.QGEO.03** [exakt, Teil bedingt oder offen] LEAN 4 FORMALISATION of QGEO.COHOM.01 (the character node of v177) + the MODULE reflection parity -- AUDIT: PASS, no sorry/admit, no domain axioms (only propext, Classical.choice, … · `lean4-carrier-rigidity`
- **FORM.QGEO.04** [exakt] FLAT-METRIC ALL-ORDERS CLOSURE FORMALISED IN LEAN (2026-06-18): TfptCarrier/SeamDeckClosure.lean (theorems diagonal_commutesClock, flat_all_orders_clock) machine-proves the v276 all-orders … · `lean4-carrier-rigidity`
- **FORM.FLATAWAY.01** [exakt] LEAN: FLAT-AWAY HEAT POSITIVE-DEFINITENESS (2026-06-18). Lean 4 formalisation (TfptCarrier/SeamDeckClosure.lean Part 5) of the algebraic core of the v292/v295 heat route: the heat-trace … · `v295_flataway_a2_exact.py`
- **FORM.SEAMEQUIV.01** [exakt] LEAN: SEAM.EQUIV.01 COMPOSITION CHAIN (2026-06-22). Lean 4 formalisation (TfptCarrier/SeamEquivChain.lean) of the v308 closing chain: the four nodes GapPositive -> InvertibleSRE -> … · `v308_seam_equiv_chain.py`
- **FORM.QGEO.BW.01** [exakt] LEAN: QGEO.SYM.01 FROM SEAM.EQUIV.01 VIA BISOGNANO-WICHMANN (2026-06-22). Lean 4 formalisation (TfptCarrier/BWKeystone.lean) of the v323 linkage: the nodes RotationInvariantVacuum -> … · `v323_bw_geometric_modular.py`
- **FORM.CARTAN.DET.01** [exakt] LEAN: THE CARRIER CARTAN DETERMINANTS AS REAL Matrix.det COMPUTATIONS (2026-06-22). Lean 4 (TfptCarrier/CartanDeterminants.lean) grounds the carrier piece of the holomorphy/anyon … · `v325_pillowcase_keystone.py`
- **FORM.MU4COMM.01** [exakt] LEAN: THE mu4 MARK-LOCALITY COMMUTATION AS REAL LINEAR ALGEBRA (2026-06-22). Lean 4 (TfptCarrier/Mu4Commutation.lean) turns the v201/v309/v323 finite-model fact (the mu4 clock rho=diag(i^n) … · `v331_necessity_of_H.py`
- **FORM.SPECTRALGAP.01** [exakt] LEAN 4 FORMALISATION of the SCALAR CORE of the universal spectral-gap principle (door #4, 2026-06-24) -- new module TfptCarrier/SpectralGapAttractor.lean, lake build OK (target … · `lean4-carrier-rigidity (SpectralGapAttractor.lean)`
- **FORM.COXETER.PRIME2.01** [exakt] LEAN 4 FORMALISATION of the RES.COXETER.SYMMETRY.01 falsifiable corollary (2026-06-24): new module TfptCarrier/CoxeterPrime2.lean machine-proves 'no prime-2-only attractor exists', … · `v409_coxeter_prime2_lemma.py`
- **FORM.SEAM.BW.HSMI.01** [exakt] LEAN 4 FORMALISATION of the Borchers/Wiesbrock standard-pair algebra (2026-06-30; the kernel-hardening of v438/SEAM.EQUIV.BW.HSMI.01) -- new module TfptCarrier/SeamStandardPair.lean, lake … · `lean4-carrier-rigidity/TfptCarrier/SeamStandardPair.lean`
- **FORM.SEAM.APPLICABILITY.01** [exakt] LEAN 4 FORMALISATION of the applicability ledger (2026-06-30; the kernel-hardening of v441/SEAM.EQUIV.APPLICABILITY.01) -- new module TfptCarrier/SeamApplicabilityLedger.lean, lake build … · `lean4-carrier-rigidity/TfptCarrier/SeamApplicabilityLedger.lean`
- **FORM.SEAM.RIGIDITY.FORCING.01** [exakt] LEAN 4 FORMALISATION of the rigidity FORCING converse (2026-06-30; the kernel-hardening of v445/SEAM.RIGIDITY.FORCING.01, also hardening v442/SEAM.RIGIDITY.UNIFORM.01) -- new module … · `lean4-carrier-rigidity/TfptCarrier/SeamRigidityForcing.lean`
- **FORM.SEAM.EDGE.CHERN.01** [exakt] LEAN 4 FORMALISATION of the edge central-charge arithmetic (2026-06-30; the kernel-hardening of the edge readings v447/v450/v451/v452) -- new module TfptCarrier/SeamEdgeChern.lean, lake … · `lean4-carrier-rigidity/TfptCarrier/SeamEdgeChern.lean`
- **FORM.SEAM.RESIDUAL.01** [exakt] LEAN 4 FORMALISATION reducing the SEAM.EQUIV.01 residual to a SINGLE named realization axiom (2026-06-30, G4 of the post-F next steps) -- new module TfptCarrier/SeamResidualAxiom.lean, lake … · `lean4-carrier-rigidity/TfptCarrier/SeamResidualAxiom.lean`
- **FORM.FLAVOR.FEEDBACK.01** [exakt] LEAN 4 FORMALISATION of the feedback normal form (2026-08-08 round 36; TfptCarrier/FlavorFeedback.lean, 16 declarations, 11 theorems) -- AUDIT: PASS, no sorry, no native_decide, no domain … · `lean4-carrier-rigidity/TfptCarrier/FlavorFeedback.lean`
- **FORM.E8.SHELL.GLOBAL.01** [exakt] LEAN 4 FORMALISATION of the global shell theorem's consequence chain (2026-08-08 round 36; TfptCarrier/GaussianShells.lean) -- AUDIT: PASS, no sorry, no native_decide, no domain axioms … · `lean4-carrier-rigidity/TfptCarrier/GaussianShells.lean`

### ARF (5)

- **ARF.SPINORCOMPILER.01** [exakt] The Arf/spinor compiler on V = L/(1+i)L = F2^4 -- Priority 1 of the Arf compiler programme, the exact finite-algebra layer (one module, 46 checks, ~2 s; discovery probe … · `v774_arf_spinor_compiler.py`
- **ARF.ROOTCLASS.01** [exakt] THE FIRED PHYSICAL KILL of the Arf compiler programme -- gaussian_class_d5_weight_purity, Priority 2, preregistered (one module, 18 checks, ~2 s; discovery probe … · `v775_gaussian_class_d5_purity.py`
- **ARF.BOUNDARY.CODE.01** [exakt, Teil bedingt oder offen] Boundary code uniqueness -- the non-circular selection chain boundary -> code -> E8 -> F2^4 -> 5 slots EXISTS with g_car = 5 as OUTPUT, and the residual premise is named (one module, 31 … · `v776_boundary_hamming_uniqueness.py`
- **ARF.VACUUM.SYNDROME.01** [exakt] v852 -- ARF.VACUUM.SYNDROME.01 + OVOID.DECODER.01: the Arf-vacuum syndrome theorem and the integral five-mode decoder -- the corpus's flat counts, orbit splits, Hecke matches and code … · `v852_arf_vacuum_ovoid.py`
- **ARF.BENT.CSS.01** [exakt] ARF.BENT.CSS.01: exact bent phase and CSS code identities; Mangoldt incidence experiment moved to Prime Front. Original module ID retained; original complete mixed record is preserved in … · `v853_bent_css_mangoldt.py`

### ARITH (4)

- **ARITH.CAPSTONE.01** [Regel oder Synthese] THE CAPSTONE: THE SM STRUCTURAL SECTOR IS Q(zeta30) + GALOIS mu4 x Z2, COMPLETE INPUT {a, pi, v_geo} (2026-06-21, synthesis of the v313-v317 arc). [E] ARITHMETIC SKELETON: the SM STRUCTURAL … · `v318_arithmetic_capstone.py`
- **ARITH.HULL.01** [exakt, Teil bedingt oder offen] The three number-field facets COMPOSE to a degree-8 = rank(E8) abelian field (2026-06-24): K=Q(i,sqrt-3,sqrt5), the arithmetic hull dual to the E8 lattice hull, composing the v390/v394 … · `v403_facet_compositum.py`
- **ARITH.TRICHOTOMY.01** [exakt, Teil bedingt oder offen] The atom trichotomy (2026-06-25): {2,3,5} are ramified/inert/split in the two exceptional CM rings, and each atom is the unique ramified prime of one quadratic facet of K=Q(i,sqrt-3,sqrt5) … · `v416_atom_trichotomy.py`
- **ARITH.NORMTRIPLE.01** [exakt, Teil bedingt oder offen] The cyclotomic norm triple (2026-06-25): the carrier split (3,2) read in the three atom-rings gives (7,13,55)=(scalaron, Delta_Q, quark numerator), and the missing CARRIER-5 clock of v417 … · `v418_cyclotomic_norm_triple.py`

### GALOIS (4)

- **GALOIS.READOUT.01** [exakt, Teil bedingt oder offen] DOES THE GALOIS mu4 x Z2 (v315) ORGANIZE THE PHYSICAL READOUTS? YES -- THE CP PHASES (2026-06-21, v315 follow-up to phenomenology). [E] FAMILY FACTOR: the hexagonal CP unit … · `v316_galois_readout.py`
- **GALOIS.FAMILY.01** [exakt, Teil bedingt oder offen] ARE THE 3 GENERATIONS A mu3/GALOIS ORBIT? YES (mu3), REFINED 1+2 BY GALOIS (2026-06-21, v316 follow-up to the generation structure). [E] GENERATIONS = CUBE ROOTS: the 3 cusp weights … · `v317_galois_family.py`
- **GALOIS.CP.PREDICT.01** [exakt, Teil bedingt oder offen] A NEW FALSIFIABLE PREDICTION FROM THE GALOIS STRUCTURE: THE TWO CP PHASES ARE LOCKED, delta_PMNS = delta_CKM,lead + pi (2026-06-22, the v316 follow-up turned testable). [E] ONE UNIT TWO … · `v320_galois_cp_relation.py`
- **GALOIS.CP.SHARPEN.01** [exakt, Teil bedingt oder offen] THE v320 CP-LOCK PREDICTION SHARPENED INTO A QUANTITATIVE PRE-REGISTERED KILL TEST (2026-06-22). [E] DISCRETE LADDER: the leading CP phase sits on the six hexagonal nodes arg(rho^k)=k*60 … · `v322_cp_lock_sharpened.py`

### CFIN (3)

- **CFIN.UNIQUE.01** [exakt] CFIN.UNIQUE.01: finite compiler uniqueness, one orbit and C6 automorphisms; cofinal Weil mirror moved to Prime Front. Original module ID retained; original complete mixed record is … · `v849_cfin_unique_cofinal_lean.py`
- **CFIN.QSTAR.NORMALFORM.01** [exakt] v880 (part 1) -- the closed census-free normal form of the distinguished quadratic refinement: q*(x) = sum_{i<j} x_i x_j + sum_i x_i = C(|x|+1, 2) mod 2, WEIGHT-ONLY with value table … · `v880_finite_anchor_theorems.py`
- **CFIN.AUT.FLAVOR.C6.01** [exakt] THE FLAVOR HEXAGON IS THE COMPILER AUTOMORPHISM ORBIT (closed 2026-08-09, round CV; promoted round CVIII in v888; probe cfin_aut_flavor_probe.py 19/19 C6FLAVOR-MEASURED + INTERTWINED): … · `v888_finite_pfaffian_signature.py`

### COX (3)

- **COX.CLOCK.01** [exakt] COXETER TOTATIVE CLOCK (2026-06-16): the mu4 clock is the order-4 character of the E8 Coxeter cycle. (Z/30)^x={1,7,11,13,17,19,23,29}=the E8 exponents=live Coxeter phases; phi(30)=8=rank … · `v223_coxeter_totative_clock.py`
- **COX.COUPLE.01** [exakt, Teil bedingt oder offen] COUPLE OR FACTORIZE? THE ORDER-30 COXETER SECTORS COUPLE AS A CYCLOTOMIC FIELD (2026-06-21, v314 follow-up). The v314 question -- does 30=g_car*(2 N_fam)=5*6 couple the carrier (5-fold, … · `v315_coxeter_coupling.py`
- **COX.STAGEA.TENSOR.01** [exakt, Teil bedingt oder offen] COXETER STAGE A CLOSED ON THE PRIMITIVE-CHARACTER SPACE (2026-07-23; 15 checks, <1 s, sympy exact): the tfpt_1 two-stage Coxeter audit asked for a natural order-30 operator with … · `v531_coxeter_tensor_stage_a.py`

### CURVE (3)

- **CURVE.CODE.2TORSION.01** [exakt] The curve/code 2-torsion decider (33/33 checks, ~9 s; discovery probe curve_code_2torsion_probe.py CURVE-CODE-PARTIAL, 2026-08-05). The frozen CM-model route through the seam curve C: y^3 = … · `v784_curve_code_2torsion.py`
- **CURVE.CODE.OUTERSPIN.01** [exakt] the canonical OUTER bridge between the code's duad side and the curve's syntheme side of the S6 duality (24/24 checks, ~11 s; discovery probe duad_syntheme_bridge_probe.py, 2026-08-05, … · `v796_duad_syntheme_bridge.py`
- **CURVE.CODE.DOILY.01** [exakt] CURVE.CODE.DOILY.01: the certified outer bridge sharpened into an incidence-operator statement with a spectral payoff, one probe (21/21 checks incl. the 3 must-fire controls, ~1 s; … · `v809_curve_code_doily.py`

### COXETER (2)

- **COXETER.PRIME2.01** [exakt, Teil bedingt oder offen] The prime-2 facet of the order-30 Coxeter clock, completing v383 (2026-06-24). The clock 30=2*3*5=h(E8) had two named number-field facets (the golden phi in Q(sqrt5), prime-5/carrier; … · `v390_prime2_facet.py`
- **COXETER.ATOMS.01** [exakt, Teil bedingt oder offen] The clockwork coherence (2026-06-24, bird's-eye): the three Coxeter primes ARE the three anchor atoms ARE the three number-field facets, and the prime-2 seam is the common dynamical factor. … · `v394_coxeter_atoms.py`

### CAS (1)

- **CAS.01** [exakt] cascade D_n=60-2n: endpoints; exponent rungs -> 240; IR tail 46080; variance 6552 · `v5_e8_cascade.py`

### DUAL (1)

- **DUAL.FRAME.01** [exakt, Teil bedingt oder offen] DUAL NORMAL FRAME (2026-06-16): the pair (d,n) is one canonical boundary compass for flavor and horizon, with CP as orientation. d=a^T R^{-1}=a^T L^{-1}=(-1/2,-1/2,1)=-1/2 (1,1,-2) … · `v225_dual_normal_frame.py`

### GAUSSIAN (1)

- **GAUSSIAN.RAMLADDER.01** [exakt] GAUSSIAN.RAMLADDER.01: the (1+i) ramification ladder -- the four documented roles of the ramified Gaussian prime pi2 = 1+i are ONE machine-checked object with four rungs, each re-certified … · `v833_gaussian_ramification_ladder.py`

### GOLD (1)

- **GOLD.ATOMS.01** [exakt, Teil bedingt oder offen] THE GOLDEN RATIO IS THE g_car=5 SIGNATURE -- THE (2,3,5) ATOMS ARE THE NETWORK SPECTRAL ANGLES (2026-06-21, v312 follow-up made exact). [E] EXACT: the affine-E8 adjacency charpoly factors … · `v313_golden_atoms_spectral.py`

### H2 (1)

- **H2.BRIDGE.01** [exakt] H2 bridge attempt: explicit per-puncture monodromies M_k computed (cusp class, prod M_k=I to 1e-5). HONEST NEGATIVE: |diag M_k|=(0,1/2,1/2), natural resolvent-dressed extraction … · `v34_h2_bridge_attempt.py`

### MCKAY (1)

- **MCKAY.E8.01** [exakt] McKAY BEDROCK (2026-06-16): WHY the three atoms are {2,3,5}. E8 is the exceptional TOP of the McKay tower of finite SU(2) subgroups (2T->affine E6 h=12, 2O->affine E7 h=18, 2I binary … · `v219_icosahedral_mckay.py`

### NORMALFORM (1)

- **NORMALFORM.CFIN.01** [exakt] v845 -- NORMALFORM.CFIN.01: the finite compiler normal form C_fin = (V, hbar, q*, sigma, iota) assembled as ONE machine-checked object from the Gaussian quotient, with the code-to-matter … · `v845_cfin_normal_form.py`

### QEC (1)

- **QEC.SEAM.01** [exakt, Teil bedingt oder offen] SEAM FINITE RECOVERABILITY CODE (2026-06-16): the seam as a finite quantum-error-correction-style code with recovery rate (2/3)^6=64/729. Code dimension = dim S^+ = 16 = 2^(g_car-1); the … · `v221_seam_qecc.py`

## Standardmodell, Eichstruktur und 4D

### GATE (46)

- **GATE.NSTAR.01** [bedingt] RESIDUAL GATE (5.0 review, B5): N_star is a reheating INPUT (observed n_s=0.9649 -> N*=57, continuous), not a compiler consequence; kept [P] · `v22_open_gates_audit.py`
- **GATE.METRIC.01** [offen] RESIDUAL GATE (5.0 review, B6/C7) -- RESCOPED 2026-06-24: the LOCAL covariant metric-sector field equation is now CLOSED parameter-free (G_ab+Lambda g_ab=(1/c3)T_ab, both coefficients … · `v22_open_gates_audit.py`
- **GATE.WOLFRAM.01** [exakt] second [C] path VERIFIED RUNNING: Wolfram Engine 14.3 reproduces alpha^-1, E8 glue, M=41, K/Q ladder, lepton c, solar angle, anchor generator, quark ratios, G2 heat-kernel … · `wolfram/tfpt_readouts.wl`
- **GATE.UWALL.02** [exakt] (U_wall) U2 kill-switch HARDENED: kill-switch fires iff D4-cusp rep reducible; across 400 random cusp-class generators M_k=U^{k-1}MU^{1-k} commutant is scalar at every point (irreducible … · `v38_uwall_killswitch.py`
- **GATE.UWALL.03** [exakt] (U_wall) SELECTOR side closed on the explicit point (algebraic, no PDE): exponents at inf eig(sum A_k)={2,1,1}=anchor a => splitting O(-2)+O(-1)^2; pardeg=0 (deg E=-4=-|mu4| + 4 weight-sums … · `v39_uwall_selectors.py`
- **GATE.UWALL.04** [exakt] (U_wall) harmonic metric = FINITE linear algebra, NOT a transcendental Hitchin PDE: polystable point is unitary (Mehta-Seshadri) => Phi=0, harmonic metric = constant invariant form H. … · `v40_harmonic_metric.py`
- **GATE.UWALL.05** [exakt] (U_wall) FINAL leg-assignment test (honest result): (i) delta=1/2 resolvent reproduces lepton c=(16/7,4/3,7/6) exactly => Lambda=(0.475,1.107,0.917); (ii) Lambda_mu=1.107>1 RULES OUT … · `v41_leg_assignment.py`
- **GATE.UWALL.15** [exakt, Teil bedingt oder offen] EXTERIOR LEG LEMMA (re-types the v41 negative as a category error): u,d share hexagon-radial (r,w)=(1,1) and the C6 resolvent reads only (r,w) => ANY scalar leg gives c_u/c_d=1 identically … · `v42_exterior_leg.py`
- **GATE.UWALL.18** [exakt] the '11' is PASCAL (third, simplest origin; no special math): Pl(K) is the exterior algebra of the family fundamental 4=mu4 (A3=SU(4) fundamental). Lambda^k(4)=C(4,k)=Pascal row … · `v45_family_exterior.py`
- **GATE.UWALL.17** [exakt, Teil bedingt oder offen] LIE-LEVEL GROUNDING of the exterior leg (no special math): the carrier half-spinor is the even exterior algebra 16=Lambda^even(5)=C(5,0)+C(5,2)+C(5,4)=1+10+5 (already in CAR.SM.01). Under … · `v44_carrier_exterior.py`
- **GATE.UWALL.16** [exakt] Lambda^2 F bridge test (honest): (1) POSITIVE - for SU(3) the exterior square is the conjugate rep, Lambda^2(M~)=conj(M~) exactly => the exterior leg IS the 3-bar (confirms the Exterior Leg … · `v43_exterior_bridge.py`
- **GATE.QGEOM.02** [exakt] Integer-lift of Q - residual sharpened to ONE named lattice invariant det Q=N_fam (honest reduction, not closure) [I]/[L]. Continuing v69: (1) DERIVED - in the mu4-puncture homology basis … · `v70_q_integer_lift.py`
- **GATE.QGEOM.01** [exakt] D4-equivariant derivation of the Q-geometry (Alessandro residual gate #2) - advances Theorem Q from [P] (algebraically forced) to [L] (geometrically derived). The 4 punctures … · `v69_d4_q_geometry.py`
- **GATE.GMETRIC.XI.01** [exakt] Seam-Engineering Index - an exact closed form for the IR-stability margin [I], plus the four-class possibility taxonomy [I]/[A]. Validates+sharpens the boundary-engineering synthesis: the … · `v63_seam_engineering_index.py`
- **GATE.METRIC.02** [exakt, Teil bedingt oder offen] Gate 2 attack - two clean tiers (was: diffuse 'full QG missing'). TIER A (physics, CLOSED under RP+gap): Decoupling Theorem - Delta_eff=Delta-2||V||=6log(3/2)-31/(4pi^2)=1.648>0 with … · `v76_gmetric_reduction.py`
- **GATE.METRIC.03** [exakt, Teil bedingt oder offen] G6 route - the seam BOUNDARY measure IS a rigorously constructed conformal net (imports G6 into existing rigor): level-1 central charges c=k dim/(k+h^v) give c(E8)_1=248/31=8=rank E8 … · `v77_e8_conformal_net.py`
- **GATE.METRIC.04** [exakt, Teil bedingt oder offen] Red-team Target A residual 1 CLOSED [L] (holomorphy is necessary AND sufficient to pin (E8)_1): at level 1 the number of primaries = det(Cartan) = |fundamental group|, giving A3=4, D5=4, … · `v83_e8net_holomorphic_uniqueness.py`
- **GATE.METRIC.05** [exakt, Teil bedingt oder offen] Red-team Target A residual (ii) CONDITIONALLY CLOSED [L] -> Target A = ONE residual: bulk-reconstruction uniqueness is NOT independent of residual (i). By the AQFT classification of full 2D … · `v87_bulk_uniqueness_reduction.py`
- **GATE.METRIC.06** [exakt, Teil bedingt oder offen] CARRIER INDEX LEMMA [L] - the Gate-A theorem gets an equivalent index-theoretic form (sharper as a target): by Kawahigashi-Longo-Mueger (mu_A=[B:A]^2 mu_B; lattice nets … · `v89_carrier_index_lemma.py`
- **GATE.METRIC.07** [exakt, Teil bedingt oder offen] GLUE UNIQUENESS [I]/[L] - the carrier extension tower completely classified (Gate A loses its last structural freedom): extensions of the carrier net/lattice correspond to isotropic … · `v92_glue_uniqueness.py`
- **GATE.WOLFRAM.02** [exakt] Second [C] path EXTENSION VERIFIED RUNNING: verification/wolfram/tfpt_readouts_extension.wl mirrors the v84-v110 round (frozen registry, master cover, reheating arithmetic, SO(16)_1 modular … · `wolfram/tfpt_readouts_extension.wl`
- **GATE.UWALL.06** [exakt, Teil bedingt oder offen] TORSION NORMAL FORM + THE DELTA = 1/2 THEOREM (executes the GATE.UWALL.05 'genericity future work'; the distinguished lepton transport value delta = 1/2 stops being a point observation): … · `v114_torsion_delta.py`
- **GATE.UWALL.07** [exakt, Teil bedingt oder offen] THE ANCHOR PINS THE RESIDUE (exact normal form for the U_wall flat-bundle residue; answers the GATE.UWALL.06 branch question): (1) MU4-AVERAGE LEMMA [I, symbolic generic X] - sum_k U^k X … · `v115_anchor_residue.py`
- **GATE.UWALL.08** [exakt, Teil bedingt oder offen] RESONANCE THEOREM (the GATE.UWALL.07 [N] findings promoted to [I] by a LINEAR computation - the anticipated Groebner step turned out unnecessary): (1) TWISTED-AVERAGE LEMMA [I, symbolic … · `v116_resonance_uniqueness.py`
- **GATE.UWALL.09** [exakt, Teil bedingt oder offen] THE MONODROMY IS THE WEYL GROUP OF A3 (the U_wall holonomy is an EXACT 24-element representation; the harmonic diagonal becomes a theorem): (1) THE EXACT MONODROMY [I]: M0 = [[0, -(1+i)/2, … · `v117_monodromy_weyl_a3.py`
- **GATE.METRIC.08** [exakt, Teil bedingt oder offen] THE GLUE Q-SYSTEM (the index-4 mu4 extension exists as an explicit Longo Q-system; the R2 gate sharpens from 'find a Q-system' to 'match the constructed one'): (1) Q-SYSTEM AXIOMS [I] on A … · `v125_glue_qsystem.py`
- **GATE.UWALL.10** [exakt, Teil bedingt oder offen] THE DUAL-NORMAL SELECTOR ((d,n) pins R columnwise; the R4' residue is bounded; external review 2026-06-12, validated and sharpened): (1) COLUMN NORMAL FORM [I]: d.c_j = a_j and n.c_j = 8 … · `v136_dual_normal_selector.py`
- **GATE.UWALL.11** [exakt, Teil bedingt oder offen] THE SELECTOR TRIANGLE (anchor + A0* zero mode + three atoms => R; the R4' residue shrinks to three pairings): (1) d IS ANCHOR DATA [I]: d = (3/2)a - 2*1 exactly - the dual anchor … · `v139_selector_triangle.py`
- **GATE.UWALL.12** [exakt, Teil bedingt oder offen] FRAME INTEGRALITY (R4' narrowed: three pairings -> two): the pairing triples (m.1, m.a, m.sigma) of INTEGER covectors against the frame (1,a,sigma) form an index-11 sublattice of Z^3, cut … · `v142_frame_integrality.py`
- **GATE.METRIC.09** [exakt, Teil bedingt oder offen] GRADED FROBENIUS (R2 identification at the finite Lie level): the Z4-graded hull (E8.GRADE.01/v128) carries (1) EXACT COSET DUALITY coset(-r) = -coset(r) mod 4 for all 240 roots (-C1=C3 … · `v143_graded_frobenius.py`
- **GATE.UWALL.13** [exakt, Teil bedingt oder offen] PAIRING VALUES ARE ATOM IDENTITIES (R4' collapses from two values to ONE map): given the lift n = w0(sigma) + |mu4| e3 (the v139 audit identity read structurally: w0 = the A3 EXPONENT … · `v145_pairing_atoms.py`
- **GATE.METRIC.10** [exakt, Teil bedingt oder offen] NS/R SECTOR CENSUS (R2 scoped honestly; corrected same day - the first version graded the Ramond labels by the sum q_D+q_A, which is not the glue label): (1) NS CENSUS [I] - every … · `v148_fock_census.py`
- **GATE.UWALL.14** [exakt, Teil bedingt oder offen] BOTH DUAL NORMALS ARE ANCHOR DATA (R4' at one discrete assignment): the pairings of the torsion normal n = (5,-9,6) with the cusp eigenbasis (v98) are EXACTLY (n.e1, n.e2, n.e3) = (6,3,5) = … · `v149_cusp_normal.py`
- **GATE.METRIC.11** [exakt, Teil bedingt oder offen] SIMPLE-CURRENT EXTENSION THEOREM (the central G_net closing theorem, assembled and stated as Definition+Theorem; external review 2026-06-13): the carrier net A = (D5)_1 (x) (A3)_1 extended … · `v154_simple_current_theorem.py`
- **GATE.METRIC.12** [exakt, Teil bedingt oder offen] QUASI-FREE BOUNDARY (Python-only; the seam-side identification of GATE.METRIC.11/v154 reduced to ONE shared physical premise): (1) RESTRICTION OF A QUASI-FREE STATE IS QUASI-FREE [I] - a … · `v155_quasifree_boundary.py`
- **GATE.METRIC.13** [exakt, Teil bedingt oder offen] THE SEAM-NET CONSTRUCTION (Route 3 carried to its exact verifiable endpoint; B=(E8)1 turned from assertion into checked computation): (1) DtN = FREE CHIRAL DISPERSION [I] - the … · `v156_seam_net_construction.py`
- **GATE.METRIC.14** [exakt, Teil bedingt oder offen] FREENESS AS THE RIGID BOUNDARY FIXED POINT (simplicity-first reframing of premise A; the free-bulk postulate becomes a fixed-point property): (1) UNIVERSAL DtN SYMBOL [I]/[L] - the … · `v157_rigid_fixed_point.py`
- **GATE.METRIC.15** [exakt, Teil bedingt oder offen] THE FREE CHIRAL c=8 FIXED POINT IS STABLE (exact operator-dimension count; the finite core of premise (A)'s reaches-the-fixed-point residual): (1) DIMENSIONS [I] - in the 16-Majorana chiral … · `v158_fixed_point_stable.py`
- **GATE.METRIC.16** [exakt, Teil bedingt oder offen] SEAM GAUSSIANITY AS A CONDITIONAL FIXED-POINT THEOREM (the proposed closure of premise (A), formalised and honestly typed; the load relocates to ONE internal statement): (1) QUASI-FREE => … · `v160_seam_gaussianity_from_pf.py`
- **GATE.METRIC.17** [exakt, Teil bedingt oder offen] THE CONE REDUCES TO ONE PARTICLE (discharges the v160/GATE.METRIC.16 scope premise via Bogoliubov second quantisation; the infinite-dimensional Schwinger cone is eliminated, a finite … · `v161_one_particle_reduction.py`
- **GATE.METRIC.18** [exakt, Teil bedingt oder offen] THE FINAL IDENTIFICATION, FORCED AS FAR AS IT GOES (the v161/GATE.METRIC.17 one-particle statement carries NO new open content; premise (A) factors into the two already-open items A2 + R1): … · `v162_seam_transport_identification.py`
- **GATE.METRIC.19** [exakt, Teil bedingt oder offen] PREMISE (A) MAXIMAL HONEST CLOSURE (the two residual gates A2 + R1 pinned to their floor; (A) carries ZERO independent open content): (1) A2 (NET EXISTENCE) CLOSES BY ASSEMBLY [C] - every … · `v165_premise_a_closure.py`
- **GATE.HOLO.01** [exakt, Teil bedingt oder offen] SEAM-HOLOMORPHY SELECTION CERTIFICATE (2026-06-17): the entire structural residual of TFPT is ONE condition with THREE provably-equivalent faces, all forcing E8 -- so P2/G_net/Target A/QGEO … · `v234_seam_holomorphy_selection.py`
- **GATE.HOLO.02** [exakt, Teil bedingt oder offen] SEAM CHERN-SIMONS REALISATION (2026-06-17): the closing step of GATE.HOLO.01 realised in standard abelian Chern-Simons / quantum-Hall language, where 'holomorphic <=> det K = 1' is a … · `v235_seam_chern_simons.py`
- **GATE.HOLO.03** [exakt, Teil bedingt oder offen] SEAM SRE CLOSURE (2026-06-17): the closing step of GATE.HOLO.01/02 recast as a PHYSICAL, falsifiable condition -- no topological ground-state degeneracy <=> det K=1 <=> the seam bulk is … · `v237_seam_sre_closure.py`
- **GATE.LEDGER.01** [exakt] THE STATUS-LEDGER INTEGRITY GUARD (2026-08-08 round 36, v878_ledger_integrity.py 10/10; a CI sentinel, no new physics). THE TRIGGER: the round-36 ledger audit found the first PRIMARY-KEY … · `v878_ledger_integrity.py`

### PS (17)

- **PS.E8BRANCH.01** [exakt] E8 HULL FORBIDS THE 126 (2026-06-17): the E8 audit hull supplies exactly {1,10,16,45} for SO(10) -- turning the scalaron-match '16-Higgs over 126-Higgs' selection into an algebraic … · `v247_e8_branching_no126.py`
- **PS.ALGEBRA.01** [exakt, Teil bedingt oder offen] PATI-SALAM NCG ALGEBRA + EXACT SM MATCHING (2026-06-17): the finite algebra A_F = H_L (+) H_R (+) M_4(C) and its consequences (Reviewer-A Theorem 2). (1) [E] unimodular gauge group SU(2)_L … · `v248_ps_algebra.py`
- **PS.RGTEST.01** [bedingt] CARRIER-NATIVE PATI-SALAM/SO(10) UNIFICATION KILL-TEST + NEGATIVE CONTROLS (2026-06-17): the suite claim-check of experiments/gauge-unification. Two-step M_Z->(SM)->M_PS->(PS)->M_GUT … · `v249_ps_unification.py`
- **PS.DIRAC.01** [exakt, Teil bedingt oder offen] FIRST-ORDER CONDITION VIOLATED BY THE MAJORANA/sigma (2026-06-17): the Chamseddine-Connes-van Suijlekom 'beyond first order' mechanism exhibited on a faithful minimal model -- advances … · `v251_first_order_sigma.py`
- **PS.DIRAC.02** [exakt, Teil bedingt oder offen] FULL 96-DIM FINITE SPECTRAL TRIPLE BUILT + NCG AXIOMS VERIFIED (2026-06-17): advances CONTRACT.QFT4D.DIRAC.01 (v250) from a contract listing 7 open obligations to a concrete object on which … · `v252_full_finite_triple.py`
- **PS.KAPPA.01** [exakt, Teil bedingt oder offen] kappa IN M_PS = kappa M_s PINNED TO O(1) ~ 1 (2026-06-17): discharges the open residual 6(b) of v249 ('exact kappa, now ~1.0-1.2, not pinned') as far as is honest. The scalaron scale M_s = … · `v253_kappa_scalaron_ps.py`
- **PS.NCG.FLUCT.01** [exakt, Teil bedingt oder offen] INNER FLUCTUATIONS Omega^1_D -> SM HIGGS DOUBLET (no 126, no junk) + sigma (2026-06-17): closes the open obligation v250 #5 (Higgs content, NOT 126), the NO-JUNK obligation, and fixes the … · `v254_inner_fluctuations.py`
- **PS.SPECACT.01** [exakt, Teil bedingt oder offen] SPECTRAL-ACTION HEAT-KERNEL EXPANSION (a0,a2,a4) + EXPLICIT kappa(f2/f0) (2026-06-17): closes the open obligation v252 #11 (spectral-action expansion) and makes the kappa-dependence on the … · `v255_spectral_action_expansion.py`
- **PS.NCG.ORIENT.01** [exakt, Teil bedingt oder offen] ORIENTABILITY + POINCARE DUALITY -- HONEST STRUCTURAL TREATMENT (2026-06-17): the two NCG obligations of v250 that do NOT close to a clean [E] on the simplified v252 representation. NO … · `v256_orientability_poincare.py`
- **PS.NCG.ORIENT.02** [exakt, Teil bedingt oder offen] ORIENTABILITY + POINCARE DUALITY -- HONEST LITERATURE-CONFIRMED RESOLUTION (2026-06-17): corrects the slightly-too-generous typing of v256. The published result (Cacic-Stephan 'Moduli … · `v257_quasiorient_modpd.py`
- **PS.DIRAC.03** [exakt, Teil bedingt oder offen] DIRAC AS MODULAR/COVARIANCE INDUCTION (2026-06-17): the finite Dirac operator D_F (v252) is the modular/covariance induction of the boundary seam quasi-free (KMS, beta=1) state restricted … · `v258_dirac_covariance_induction.py`
- **PS.SPECACT.02** [exakt, Teil bedingt oder offen] MODULAR/KMS SPECTRAL-ACTION CUTOFF (2026-06-17): the spectral-action cutoff f is the seam's own modular/KMS weight (beta=1 by Tomita-Takesaki/BW, v239; 2pi=1/(4c3), v58), so the heat-kernel … · `v259_modular_cutoff_kappa.py`
- **PS.UVBRANCH.01** [Regel oder Synthese] PATI-SALAM AS A CONDITIONAL UV BRANCH (2026-06-18): carrier-native Pati-Salam is permitted as a falsifiable UV completion, NOT the default -- E8-allowed content {1,10,16,45} (no 126), sin^2 … · `v265_qft4d_fork_freeze.py`
- **PS.THRESHOLD.01** [Regel oder Synthese] PATI-SALAM THRESHOLD MODEL EXPLICIT (2026-06-18): the PS UV branch is admissible only with an explicit two-step threshold model meeting the scalaron scale (M_PS=kappa M_s, kappa in … · `v265_qft4d_fork_freeze.py`
- **PS.PROTON.KILL.01** [Kill-Test] PROTON-DECAY SAFETY RULE (2026-06-18): no discovery tau_p is predicted; PS is killed if its thresholds need forbidden reps, leave the kappa band, or predict tau_p below the current bound; … · `v265_qft4d_fork_freeze.py`
- **PS.THRESHOLD.02** [exakt, Teil bedingt oder offen] EXPLICIT PATI-SALAM THRESHOLD MODEL (2026-06-18; recount 2026-11-24, v1036): the two-step carrier-native PS chain made explicit for the E8-allowed contents, at one and two loops. Minimal … · `v266_ps_threshold_proton.py`
- **PS.PROTON.02** [Kill-Test] PROTON DECAY EXCLUDES THE v266 CONTENTS OF BRANCH B (2026-06-18; corrected 2026-11-24, v1036): d=6 tau_p(p->e+pi0) = 1e36 (M_GUT/1e16)^4 (alpha_U^-1/40)^2 yr with the x3 hadronic band, … · `v266_ps_threshold_proton.py`

### QFT4D (13)

- **QFT4D.MATTER.01** [exakt, Teil bedingt oder offen] SPECTRAL-ACTION MATTER CONTENT (2026-06-17): the computable core of CONTRACT.QFT4D.01 (v244) -- the NCG matter representation, discharging the representation-theory + tree-relation core to … · `v245_spectral_matter_content.py`
- **QFT4D.RGTEST.01** [bedingt] UNIFICATION DATA CROSS-CHECK (2026-06-17): the honest data-confrontation kill-test for the spectral-action prediction (CONTRACT.QFT4D.01, v244/v245). Does the boundary condition g1=g2=g3, … · `v246_unification_data_crosscheck.py`
- **QFT4D.FORK.01** [Regel oder Synthese] 4D-QFT FORK FROZEN AS A BRANCH POLICY (2026-06-18): boundary-only is the canonical 4D reading; carrier-native Pati-Salam is a falsifiable UV branch; SM-only 4D-GUT unification is killed. A … · `v265_qft4d_fork_freeze.py`
- **QFT4D.BOUNDARY.DEFAULT.01** [Regel oder Synthese] CANONICAL 4D READING = BOUNDARY ONLY (2026-06-18): TFPT delivers the boundary QFT, the relative spectral action and the SM readout; a 4D-GUT is not claimed by default (E8 is the audit hull, … · `v265_qft4d_fork_freeze.py`
- **QFT4D.SMONLYGUT.01** [Kill-Test] SM-ONLY 4D-GUT UNIFICATION KILLED (2026-06-18): 1-loop SM running with b=(41/10,-19/6,-7) does not unify -- at the alpha1=alpha2 scale alpha3^-1 is off by ~5.6, spread ~8.8 at 2e16 GeV. The … · `v265_qft4d_fork_freeze.py`
- **QFT4D.NO_OVERCLAIM.01** [exakt] 4D OVER-CLAIM WORDING GUARD (2026-06-18): bare 4D over-claims ('TFPT solves 4D QFT', 'SM couplings unify', '4D GUT prediction', 'unbroken E8 gauge') are absent from the active docs, and the … · `v265_qft4d_fork_freeze.py`
- **QFT4D.SPERT.01** [bedingt] PERTURBATIVE 4D S-MATRIX (S_pert) AS AN EPSTEIN-GLASER / pAQFT SKELETON (2026-06-18): the honest middle layer of the three-layer S-matrix -- S_top (2d DHR braiding monodromy |M|=1, … · `v269_spert_paqft_skeleton.py`
- **QFT4D.SPERT.02** [exakt, Teil bedingt oder offen] CONCRETE EPSTEIN-GLASER ONE-LOOP QUARTIC RENORMALIZATION (2026-06-18): takes the v269 S_pert pAQFT skeleton from 'exists' to an actual EG number (no path integral). [E] T_1=:L_int: needs no … · `v271_eg_oneloop_quartic.py`
- **QFT4D.SPERT.03** [exakt, Teil bedingt oder offen] EG ONE-LOOP GAUGE SELF-ENERGY -> (b1,b2,b3) = (41/10,-19/6,-7) (2026-06-18): takes v271's quartic mechanism to the gauge sector. [E] the one-loop gauge 2-point function is the order-2 EG … · `v273_eg_gauge_running.py`
- **QFT4D.SPERT.04** [exakt, Teil bedingt oder offen] S_pert -> S_phys LSZ BRIDGE + ONE-LOOP UNITARITY (2026-06-18): connects the EG perturbative S-matrix (v269/v271/v273) to the physical asymptotic S-matrix (S_phys = LSZ on the … · `v278_lsz_bridge_unitarity.py`
- **QFT4D.EG.ALLORDER.01** [bedingt] The all-order Epstein-Glaser/BRST contract for the perturbative 4D S-matrix S_pert (2026-06-23, review Point 3): the 'all-order closing statement' for the 4D perturbative QFT leg, typing … · `v381_qft4d_eg_allorder.py`
- **QFT4D.OS.RECON.01** [offen] THE 4D OSTERWALDER-SCHRADER PACKAGE CONTRACT [O] (registered 2026-08-27; step (3) of the DIMENSION.UPLIFT.FIREWALL.01 route ordering). INPUT: TFPT fixes the discrete data (representation … · `-`
- **QFT4D.LATTICE.FUNDAMENTAL.01** [Regel oder Synthese] THE LATTICE-FUNDAMENTAL READING (registered 2026-08-28 as a typed DECISION, not an [E] claim and not a marker move). DECISION: TFPT's PHYSICAL completeness criterion is the finite, local, … · `v989_lattice_gate_battery.py`

### CHIRAL4D (5)

- **CHIRAL4D.NOMIRROR.01** [offen] THE CHIRAL 4D LATTICE MATTER CONTRACT [O] (registered 2026-08-27; step (4) of the firewall route ordering). The seam's c_- = 8 and Chern/inflow structure are a boundary-ANOMALY signal, not … · `v989_lattice_gate_battery.py`
- **CHIRAL4D.MIRROR.DET16.01** [bedingt] THE DET16 RANK-ONE MIRROR PROJECTOR CANDIDATE (registered 2026-08-29, review wave 5, v1002). Candidate [C]: exact finite-cluster gapping mechanism on the full 2^16 chiral Spin(10) Fock … · `v1002_spin10_mirror_projector.py`
- **CHIRAL4D.MIRROR.DFILTER.BLOCK.01** [exakt, Teil bedingt oder offen] DERIVATIVE-FILTERED WALL BLOCK LEMMA (registered 2026-09-05, v1023; no parent marker move). For H_df=Delta b^dagger b+lambda(b+g D_GW p)^dagger(b+g D_GW p), every singular-value block has … · `v1023_derivative_filtered_wall.py`
- **CHIRAL4D.MIRROR.SIGNED.CAR.01** [exakt, Teil bedingt oder offen] Even onsite unitary conjugation constructs canonical DET excitations under the determinant-singlet premise. The additive Delta*N_a replaces the flat Delta*Q. Signed wall det … · `v1027_signed_det_car_wall.py`
- **CHIRAL4D.MIRROR.FACTORIZED.BOUND.01** [exakt, Teil bedingt oder offen] The actual v1008 factorized onsite model at n=4 J=2 and |d|<=3/5 has uniform gap >=1+sqrt(3)/2. An ordinary hopping bond has zero vacuum expectation but nonzero action and therefore cannot … · `v1034_parent_mirror.py`

### DIMENSION (2)

- **DIMENSION.UPLIFT.FIREWALL.01** [Regel oder Synthese] THE DIMENSION-UPLIFT FIREWALL (registered 2026-08-27 from the external deep review: 'chiral 1+1D seam net does NOT imply a local, interacting 3+1D quantum world'). STANDING RULE, in force … · `-`
- **DIMENSION.SELECTOR.4D.01** [bedingt] THE CONDITIONAL 4D DIMENSION SELECTOR [C] (registered 2026-08-27 from the external bird's-eye review, the master-route wave; machine-checked in v975, 16/16 sympy-exact). STATEMENT: under … · `v975_dimension_selector_4d.py`

### SM (2)

- **SM.RG.01** [bedingt] dimensionful EW/QCD masses (m_W; m_Z; m_H; sin^2 thetaW; alpha_s) · `-`
- **SM.QFT.01** [exakt, Teil bedingt oder offen] admissible-sector QFT measure (gap 6 log 3/2; RP from cusp; OS) · `-`

### GAUGE (1)

- **GAUGE.DETLINE.FIXPOINT.01** [offen] THE GAUGE-COUPLING FIXPOINT CONTRACT [O] (registered 2026-08-28, review wave 3). DEMAND: all three gauge couplings (g1, g2, g3 / equivalently alpha, sin^2 theta_W, alpha_s) arise as … · `v996_gauge_grammar_kill.py`

### HIGGS (1)

- **HIGGS.FREESEAM.01** [exakt, Teil bedingt oder offen] HIGGS QUARTIC FROM THE FREE SEAM (a derived boundary condition tying m_H to premise (A); validated at two loops): the seam UV is the free chiral c=8 fixed point (v157/v158, Gaussian, no … · `v166_higgs_free_seam.py`

### PARENT (1)

- **PARENT.GAUSS.D4.ROUND4.01** [exakt, Teil bedingt oder offen] Full 192-mode site exterior character and exact mod5 Frobenius reduction give Gauss dimension3 mod5 for one site and0 mod5 for the open 2x2x2 cube with the specified nine-irrep link … · `v1028_gauss_d4_cap.py`

### QCD (1)

- **QCD.LAMBDA.01** [exakt, Teil bedingt oder offen] QCD CONFINEMENT SCALE PARAMETER-FREE FROM THE CARRIER b3 = -7 (the m_p/m_e numerator is independently sourced; the ratio itself stays an honest non-claim): (1) running alpha_s(M_Z) = 0.1179 … · `v164_qcd_scale.py`

### TFPT4D (1)

- **TFPT4D.LATTICE.ACTION.01** [offen] THE EXPLICIT 4D LATTICE ACTION CONTRACT [O] (registered 2026-08-27 from the external bird's-eye review; route step (2) of the DIMENSION.UPLIFT.FIREWALL.01 ordering made concrete -- THE … · `v989_lattice_gate_battery.py`

## Feinstruktur und Elektromagnetismus

### ALPHA (9)

- **ALPHA.QUILLEN.01** [exakt, Teil bedingt oder offen] THE ALPHA FIXED POINT AS THE STATIONARITY OF A U(1) DETERMINANT LINE (2026-06-22): reformulates the EM closure from 'a very cleverly built fixed point' into 'the stationarity of a named … · `v341_alpha_quillen.py`
- **ALPHA.QUILLEN.EXACT.01** [exakt, Teil bedingt oder offen] Names the Quillen determinant-line VARIATION as a tracked [O] target (2026-06-23, review Point 10): delta_tau(log det_zeta Delta_U(1) + 8 b1 c3^6 log phi_seam)=0 on the determinant line … · `v382_alpha_quillen_exact.py`
- **ALPHA.QUILLEN.PROGRESS.01** [exakt, Teil bedingt oder offen] An HONEST attempt at the external target ALPHA.QUILLEN.EXACT.01 (door #5, 2026-06-24): it REDUCES the open step to a sharper named sub-target via a solvable model, but does NOT close it -- … · `v391_alpha_quillen_progress.py`
- **ALPHA.QUILLEN.PROGRESS.02** [exakt, Teil bedingt oder offen] A SECOND HONEST step on the external target ALPHA.QUILLEN.EXACT.01 (2026-06-27, v433): it GROUNDS and CONNECTS, it does NOT close it -- ALPHA.QUILLEN.EXACT.01 stays [O] (external math, type … · `v433_alpha_quillen_heatkernel.py`
- **ALPHA.QUILLEN.PROGRESS.03** [exakt, Teil bedingt oder offen] A THIRD HONEST step on the external target ALPHA.QUILLEN.EXACT.01 (2026-06-28, v434): it settles the status of the three residuals named after v433 and shows they are NOT three independent … · `v434_alpha_quillen_betafunction.py`
- **ALPHA.QUILLEN.PROGRESS.04** [exakt, Teil bedingt oder offen] A FOURTH HONEST step on the external target ALPHA.QUILLEN.EXACT.01 (2026-06-29, v435): it attacks the SINGLE remaining [O] left after v434 -- residual (2), the cubic alpha^3 Chern/Maxwell … · `v435_alpha_quillen_chernlevel.py`
- **ALPHA.QUILLEN.INFLOW.01** [exakt, Teil bedingt oder offen] THE ALPHA^3 LEVEL = THE COMPUTED BULK CHERN INVARIANT + THE SEAM F-NORMALISATION = THE AFFINE EMBEDDING INDEX (2026-07-02, closure-route round) -- upgrades BOTH residuals left after … · `v470_alpha_inflow_level.py`
- **ALPHA.QUILLEN.DETLINE.01** [exakt] THE OCCUPIED BERRY LINE OVER THE U(1)-TWIST MODULI HAS CHERN NUMBER = THE BLOCH RESPONSE LEVEL (2026-07-03) -- the finite Berry-Chern response, COMPUTED at the collar model; the … · `v472_quillen_detline_moduli.py`
- **ALPHA.QUILLEN.FACES.01** [exakt, Teil bedingt oder offen] THE TWO NAMED FACES OF ALPHA.QUILLEN.EXACT.01 CARRY COMPUTED WITNESSES (2026-08-27, seam/alpha promotion wave; probes seam_bfk_conical_det_probe 21/21 SPEC d06b3bbec70e7efc / … · `v974_alpha_faces_computed.py`

### EM (5)

- **EM.FP.01** [exakt] alpha^-1 = 137.0359992168 unique root of F_U1=0 (existence + uniqueness + interval-arithmetic + ablation) · `v3_em_alpha.py`
- **EM.BUDGET.01** [exakt] EM budget M=41 is FORCED = 10 b1 = sum L + N_Phi (U(1) index, not a free fit); (4/5)M c3^6 = 8 b1 c3^6 · `v13_open_gates.py`
- **EM.WARD.01** [exakt, Teil bedingt oder offen] EM Boundary Ward decomposition (Theorem C): F_U(1)=0 read as reduced boundary U(1) Ward identity in three terms - Maxwell alpha^3 + Calderon -2c3^3 alpha^2 (2=|Z2|) + transport -8b1 c3^6 … · `v48_em_ward.py`
- **EM.WARD.02** [exakt, Teil bedingt oder offen] HEAT-KERNEL SHARPENING OF EM.WARD.01 (2026-06-22): derives the ORIGIN and STRUCTURE of the F_U(1) determinant-line terms from textbook Seeley-DeWitt/Gilkey heat-kernel coefficients, … · `v342_em_ward_heatkernel.py`
- **EM.BUDGET.02** [exakt, Teil bedingt oder offen] PYRATE GAUGE CROSS-CHECK (independent third-party confirmation of EM.BUDGET.01/CAR.SM.01 and an F_transfer gauge-input check): the carrier/SM field content gives, by the textbook one-loop … · `v159_pyrate_gauge_crosscheck.py`

## Flavor, Massen und CP

### FLAV (63)

- **FLAV.R.01** [exakt] residue matrix R: det=8; minors {2;3;5}; char poly; SNF diag(1;1;8) (positive invariant factors); ||R||^2=78; sum L=40 · `v4_flavor_matrix.py`
- **FLAV.KQ.01** [exakt] Q-Sigma algebra: K=R+Q*Sigma; L=R+Q*(I+Sigma); chi_Q=(t-1)(t^2-5t+3); chi_K=(t-1)(t^2-8t+4); det ladder (3;4;8;20) product 1920=|W(D5)|; a^T K a=41 · `v10_projection_involution.py`
- **FLAV.QGEO.01** [bedingt] geometric origin of Q reduced: Q_+ = A3 exponent grading (spec 1;2;3); Q_-^2 = N_fam on support; open = realise Q_+/Q_- on parabolic P1\mu4 UPDATE (2026-07-31, v567, QGEO.QRECON.01): the … · `v13_open_gates.py`
- **FLAV.H2.01** [exakt] H2 splitting type 4=1+1+2; coker R = Z8 (branch fixed) · `v4_flavor_matrix.py`
- **FLAV.TH12.01** [exakt] solar angle sin^2 th12 = 1/3 - phi0/2 = 0.30675 via explicit mu-tau texture · `v9_neutrino_texture.py`
- **FLAV.SOLAR.01** [exakt] dual anchor a^T R^-1 = a^T L^-1 = (-1/2,-1/2,1) exact; eps=(N_fam/|mu4|)phi0=3/4 phi0 = c3 + 36 c3^4 (leading term c3, fourth-order puncture correction exact); full PMNS readout · `v16_solar_dual_anchor.py`
- **FLAV.RESOLV.01** [exakt] finite hexagonal resolvent D_y^-1 = (sum y^{5-m} delta^m U6^m)/(y^6-delta^6) proven; quark c-digits need U_f* holonomy (open finite step) · `v17_hexagonal_resolvent.py`
- **FLAV.QUARK.01** [exakt] charged-fermion Yukawas y=lambda_Y^L*Lambda reproduce quark source ratios exactly (mu/md=0.470085, mc/ms=13.61, mt/mb=40.80), lepton c=(16/7,4/3,7/6) to Lambda rounding, full hierarchy; … · `v18_quark_yukawa.py`
- **FLAV.UFSTAR.01** [exakt] exact pole delta_ph=z_*^(1/6), z_*=(794-7sqrt9961)/2187 (2187=3^7); bare resolvent fails lepton ratios 3-11x (holonomy essential); required Lambda=(0.475,1.107,0.917), Lambda_mu>1 => … · `v19_monodromy_moduli.py`
- **FLAV.LEPTONC.01** [exakt] lepton c's DERIVED: non-anchor c=|mu4|^w/(5/4-cos(r pi/3)) at delta=1/2,y=1 gives c_e=16/7,c_mu=4/3 exactly; anchor from product rule prod c=2^gcar/Nfam^2=32/9 gives c_tau=7/6; reproduces … · `v20_lepton_c_derivation.py`
- **FLAV.SOLAR.02** [exakt] solar coefficient is the A3 glue-norm: eps=q(A3)*phi0=(3/4)phi0, q(A3)=3/4 the same norm in q(D5)+q(A3)=2; leading term = c3 exact; sin^2 th12=1/3-phi0/2=0.30675. Product rule … · `v21_solar_product_quark.py`
- **FLAV.QRATIO.01** [exakt] quark RATIO closure: c_u/c_d=g_car*11/(Nfam^2 DeltaQ)=55/117, c_c/c_s=p5(a)/(Oadm-1)=34/47, c_t/c_b=Nfam/(2 DeltaQ)=3/26; mass ratios 0.470085/13.605/40.81 match targets <0.03%; rational c … · `v24_quark_ratio_closure.py`
- **FLAV.CKM.01** [exakt] CKM core s12=lambda_C; s23; s13; delta=pi/3+3lC^2 · `-`
- **FLAV.PLUCKER.01** [exakt] anchor-plane Plucker apparatus (anti-numerology): left/right 2x2 minors of B_M=(1^T M; a^T M) and C_M=(M1 Ma) for R,Q,K,L are compiler atoms; ||Pl(K)||_1=11 (localizes the '11' in K), … · `v37_plucker_anchor.py`
- **FLAV.PENCIL.01** [exakt] K+xQ pencil: L=K+Q => mass and transport are one pencil P(x)=K+xQ with det(K+xQ)=3x^3+7x^2+6x+4=(N_fam, scalaron 7, |R^+(A3)|, |mu4|); chi_K=t^3-9t^2+12t-4; [K,L]=[L,Q]=[K,Q], … · `v37_plucker_anchor.py`
- **FLAV.RIGID.01** [exakt] Readout Rigidity (Theorem U2): c_u/c_d is CONSTANT on the discrete selector stratum S_{8,Q}={det R=8, SNF(R)=(1,1,8), Spec(Q+)={1,2,3}}, independent of the continuous D4 position => … · `v49_readout_rigidity.py`
- **FLAV.QGEO.02** [exakt, Teil bedingt oder offen] Q geometry (Theorem Q): Sigma-split Q_+=(Q+SigQSig)/2=[[3,0,0],[0,2,0],[0,2,1]] (chi=(t-1)(t-2)(t-3), Spec={1,2,3}=A3 exponents), Q_-=[[0,1,0],[3,0,0],[3,0,0]] (chi=t(t^2-3), … · `v50_q_geometry.py`
- **FLAV.MASSVOL.01** [exakt] Grand Mass Volume = the absolute mass-scaling normalisation [I] (NOT a [P] product rule): the individual c_q are a gauge, but the sector mass-matrix determinants are clean phi0-powers = … · `v46_grand_mass_volume.py`
- **FLAV.DELTA.01** [exakt] delta=1/2 glue-norm origin + (g_car,N_fam)/|mu4| unification: the two E8-glue discriminant norms ARE the carrier rank and family count over the glue index, q(D5)=g_car/|mu4|=5/4 (rank … · `v51_boundary_half_step.py`
- **FLAV.KQ.02** [exakt] K+xQ pencil endpoints (audit ext. of v37): P(x)=det(K+xQ)=3x^3+7x^2+6x+4 anchored at four integer nodes by compiler atoms: P(-1)=2=|Z2| (K-Q sheet endpoint), P(0)=4=det K=|mu4|, P(1)=20=det … · `v52_pencil_endpoints.py`
- **FLAV.RIGID.02** [exakt] SIMPLE R-bridge closure (the much-simpler solution, no transcendental machinery): the selector stratum S_{8,Q}={det R=8, SNF(R)=(1,1,8), Spec(Q+)={1,2,3}} is now FULLY DERIVED - det … · `v71_simple_r_bridge.py`
- **FLAV.QGEO.03** [exakt] det Q = N_fam DERIVED from the cusp class (Alessandro Q-residual closed to the same geometry): the integer-lift datum det Q=|coker Q|=3 is the ORDER of the cusp monodromy class (eigenvalues … · `v72_q_det_from_cusp.py`
- **FLAV.UPOINT.01** [exakt] Gate 1 COMPLETE - U_point reduces to v_geo (the absolute flavor amplitude is not an independent gate): quark ratios closed (Plucker, v49/v71) + lepton c exact (16/7,4/3,7/6, v20) + each … · `v75_upoint_to_vgeo.py`
- **FLAV.PENCIL.02** [exakt, Teil bedingt oder offen] Operator-pencil geometry of (K,Q,R,L) - the anchor plane as a small area in operator space (validated review block; all in admissible classes). STRONG: [2] ANCHOR SINGULARITY THEOREM det … · `v80_operator_pencil_geometry.py`
- **FLAV.PENCIL.03** [exakt, Teil bedingt oder offen] The mass-transport pencil P(x)=K+xQ as a projective line with two anchor-block singularities, and the DOUBLE COVER branched over them (validated review block; all in admissible classes). … · `v81_singular_pencil_matrices.py`
- **FLAV.PENCIL.04** [exakt, Teil bedingt oder offen] Anchor-placement splitting trichotomy - the clean rational double cover is NON-GENERIC (anchor-first hardening) [I] + [P] (sector id). The splitting-type exponents (2,1,1) (= the anchor … · `v82_koide_attractor_splitting.py`
- **FLAV.PENCIL.05** [exakt, Teil bedingt oder offen] MASTER-COVER THEOREM (answers next.txt P3; resolves P6 negatively; types P4) [I]/[L]: the anchor-block determinant is GL(2)-COVARIANT on the pencil plane span{K,Q} - det B((aK+bQ)+x(cK+dQ)) … · `v85_master_cover.py`
- **FLAV.CP.01** [exakt, Teil bedingt oder offen] Red-team Target D follow-up: CP-phase residual QUANTIFIED [N] (not closed; freeze respected): frozen leading reading delta_CKM=pi/3+3 lam^2=1.198232 rad=68.654 deg (registry DELTA_CKM_RAD … · `v88_cp_phase_audit.py`
- **FLAV.DIAMOND.01** [exakt] SHEET DIAMOND & WINDING LINE (second review round 2026-06-11; arithmetic fully validated, interpretations corrected where needed) [I]: (1) all four flavor operators are distinguished points … · `v94_sheet_diamond.py`
- **FLAV.PLUCKER.02** [exakt] RIGHT-PLUECKER CANONICALISATION of the quark-ratio integers (load-bearing improvement; supersedes the loose-26 readings as the canonical source): with Pi_R(M) = 2x2 minors of (M1 | Ma): … · `v94_sheet_diamond.py`
- **FLAV.DIAMOND.02** [exakt] CENTERED FLAVOR DIAMOND (third review round 2026-06-11; all arithmetic validated incl. SNF) [I]: the v94 sheet diamond is a CENTERED CROSS - with the center … · `v95_centered_diamond.py`
- **FLAV.SHEET.01** [exakt, Teil bedingt oder offen] BRANCH-KERNEL SELECTION (next.txt P1 unblocked): at each branch point of the anchor-block cover the 2x2 anchor block is RANK 1 with canonical INTEGER kernels (Koide -2/3: right … · `v96_branch_kernel_selection.py`
- **FLAV.SHEET.02** [exakt, Teil bedingt oder offen] SHEET-CONJUGATION BRIDGE (P1 -> ONE [P] residual): the family-side conjugation (cusp negation 0 fixed, 1/3<->2/3 on the Q_+ eigenspaces e1=(0,0,1), e2=(0,1,2), e3=(1,0,0)) has a UNIQUE … · `v97_sheet_conjugation_bridge.py`
- **FLAV.SHEET.03** [exakt] DISCRIMINANT DICTIONARY DERIVED (P1 closed modulo GATE.QGEO): the generation space carries an INTEGER mu4 - G = T_A*Sigma acts on the cusp basis as G e1=-e1, G e3=e2, G e2=-e3 (the v69 B1+E … · `v98_discriminant_dictionary.py`
- **FLAV.H2.02** [exakt, Teil bedingt oder offen] THE HEXAGON IS THE SIGN-TWISTED FAMILY SPECTRUM (first exact piece of the H2 dictionary; the v20 lepton coefficients are resolvent determinants of the exact W(A3) monodromy GATE.UWALL.09): … · `v118_hexagon_family_dictionary.py`
- **FLAV.H2.03** [exakt, Teil bedingt oder offen] THE ADDRESS TABLE (the open FLAV.H2.02 item structured; exact arithmetic on FROZEN v18/v20 integers): (1) THE LEPTON WORDS ARE THE COMPILER ATOMS [I]: (L_e, L_mu, L_tau) = (8, 5, 3) = (rank … · `v120_address_table.py`
- **FLAV.H2.04** [exakt, Teil bedingt oder offen] THE ADDRESS TABLE IS PINNED (the FLAV.H2.03 per-fermion item closed at the information level; exact census theorem): (1) THE WORD TABLE IS THE RESIDUE OPERATOR [I]: the established identity … · `v121_address_pinning.py`
- **FLAV.H2.05** [exakt, Teil bedingt oder offen] THE MARGINS ARE THEOREMS (the FLAV.H2.04 residue closed: the established frozen selectors pin R uniquely - the atom margins lose their input status): THREE ESTABLISHED SELECTORS, all frozen … · `v122_margin_theorem.py`
- **FLAV.DUAL.01** [exakt, Teil bedingt oder offen] THE DUAL ANCHOR LEMMA (the inverse flavor response of the anchor IS the Nariai root; the third purely algebraic leg of the flavor-horizon bridge; external review 2026-06-12, validated): (1) … · `v134_dual_anchor.py`
- **FLAV.SURF.01** [exakt, Teil bedingt oder offen] THE DETERMINANT SURFACE (the sheet diamond is a 2d surface with iso-volume walls at the discrete atoms; external review 2026-06-12, validated): (1) SURFACE + WALLS [I]: on the established … · `v135_det_surface.py`
- **FLAV.QGEO.04** [exakt, Teil bedingt oder offen] THE Q+ GRADING IS COHOMOLOGY (the R5/GATE.QGEO dictionary at the spectrum level; external review 2026-06-12, validated and sharpened): (1) BASIS + CHARACTERS [I]: b_1(P^1 minus mu_4) = 3 = … · `v137_qplus_cohomology.py`
- **FLAV.QGEO.05** [exakt, Teil bedingt oder offen] THE CANONICAL MAP EXISTS AND IS RIGID (the GATE.QGEO naturality residue collapses to ONE Z_3 choice): (1) THREE mu_4 ACTIONS ONE MATCH [I]: the deck U = diag(1,i,-i) has characters (0,1,3) … · `v140_canonical_map.py`
- **FLAV.QGEO.06** [exakt, Teil bedingt oder offen] DECK SELECTION THEOREM (R5's Z_3 choice DERIVED, executing the v140 residue): the established integer deck G=T_A*Sigma (v97 anchor-forced conjugation; v98 integer mu4 on generation space) … · `v141_deck_selection.py`
- **FLAV.QGEO.07** [exakt, Teil bedingt oder offen] MOEBIUS D4 REALISATION (R5 at its floor): the FULL dihedral symmetry group of the seam curve <delta: z->iz, iota: z->1/z> = D4 (order 8) acts on H^1(P^1 minus mu4) EXACTLY as the integer … · `v146_moebius_d4.py`
- **FLAV.RGSTAB.01** [exakt, Teil bedingt oder offen] F_TRANSFER SCALE-TAGGING (the lepton-mixing predictions are RG-STABLE outputs; the quark mass-ratio readouts carry a scale tag): (1) the PMNS angles run only through the charged-lepton Y_e … · `v163_rg_stability_flavor.py`
- **FLAV.STRONGCP.02** [exakt] STRONG CP AS PFAFFIAN REALITY (theta_eff=0 as a model-level CAR/Pfaffian theorem; formalises the tfpt_2 structural chain): in the self-dual CAR (16-Majorana) seam model, (1) the sheet-odd … · `v173_pfaffian_cp_car.py`
- **FLAV.PMNS.02** [exakt, Teil bedingt oder offen] M_nu AS THE TYPE-I SEESAW OF D_F (2026-06-18): the light neutrino matrix is the seesaw of the finite Dirac operator's own blocks, M_nu=-m_D M_R^{-1} m_D^T, advancing the v9 solar-angle … · `v263_mnu_from_df.py`
- **FLAV.TH13.01** [exakt, Teil bedingt oder offen] REACTOR ANGLE EXPONENT = CARRIER HYPERCHARGE TRACE (2026-06-18): the theta_13 value sin^2 theta_13 = phi0 e^{-5/6} = 0.0231 (v16) has its EXPONENT origin closed exactly -- gamma = tr_E Y^2 … · `v268_theta13_carrier_trace.py`
- **FLAV.PMNS.03** [offen] COMPLETE COMPLEX PMNS MATRIX + LEPTONIC JARLSKOG J_PMNS (CP STRENGTH) (2026-06-18): the four TFPT channels (sin^2 th12 = 1/3 - phi0/2 v9/v263, sin^2 th23 = 1/2 v9/v263, sin^2 th13 = phi0 … · `v270_pmns_jarlskog_assembly.py`
- **FLAV.NUSCALE.01** [exakt, Teil bedingt oder offen] ABSOLUTE NEUTRINO-MASS SCALE TYPED HONESTLY (2026-06-18): does NOT fabricate a Sigma m_nu number (would contradict the v184 firewall); closes the PMNS bookkeeping instead. [E] the absolute … · `v272_nu_mass_scale.py`
- **FLAV.TH13.PRESSURE.01** [exakt] THE theta13 PRESSURE POINT -- THE HONEST WEAK SPOT, QUANTIFIED AND PRE-REGISTERED (2026-06-22). [N] TENSION: sin^2 theta13 = phi0 e^{-5/6} = 0.023108 vs NuFIT 6.0 NO (0.02195 +- 0.00058) is … · `v328_theta13_pressure.py`
- **FLAV.PMNS.CLOSE.01** [exakt, Teil bedingt oder offen] The PMNS dynamics closed modulo the tau=omega keystone (2026-06-24): the CP-phase ORIGIN is the family CM point (v405) and the absolute scale is one v_geo-class unit, so PMNS has zero free … · `v406_pmns_phase_origin.py`
- **FLAV.SELECTOR.CLOSE.01** [exakt, Teil bedingt oder offen] The R4' selector residue folds into the tau=omega keystone (2026-06-24): the three (d,n) selector pairings ARE the tau=omega family-slice atoms, so deriving the selector pair reduces to the … · `v407_dn_pairings_omega.py`
- **FLAV.UD.VPOWER.01** [exakt, Teil bedingt oder offen] The quark u/d ratio is a pure V-power readout (2026-06-24): the established c_u/c_d=g_car*||Pl(K)||_1/(N_fam^2 Delta_Q)=55/117 (v94) re-expresses EXACTLY via the sheet generator V (v410). … · `v411_ud_ratio_vpower.py`
- **FLAV.THIRDGEN.PATTERN.01** [exakt, Teil bedingt oder offen] THIRD-GENERATION LADDER-BASE PATTERN (2026-07-02, v467): the suite's THREE ~2 sigma mixing tensions (sin^2 theta13 +2.0 sigma, sin^2 theta23 +1.8 sigma, V_cb +1.8 sigma; v307) are ONE … · `v467_thirdgen_mixing_pattern.py`
- **FLAV.DM2RATIO.01** [exakt, Teil bedingt oder offen] SPLITTING RATIO = LEPTONIC JARLSKOG (2026-07-02, v468): the neutrino mass-squared splitting ratio equals the Jarlskog invariant of the FROZEN PMNS matrix -- dm21^2/dm31^2 = |J_PMNS| = … · `v468_dm2_ratio_jarlskog.py`
- **FLAV.NUSCALE.02** [bedingt] ABSOLUTE NEUTRINO SCALE UNDER THE CARRIER YUKAWA NORMALISATION (2026-07-14, bird's-eye round; CANDIDATE class, genre v467/v468): FLAV.NUSCALE.01 typed the absolute scale as one open UV … · `v481_seesaw_carrier_ladder.py`
- **FLAV.NUSCALE.03** [Regel oder Synthese] THE v481 DECISION COMPUTATION, EXECUTED (2026-07-14, deep-dive round; bracketed / assumption-robust form): the integer rung M_R = c3^3 Mbar is EXCLUDED as an UNSTRUCTURED exact pin, … · `v482_seesaw_rung_decision.py`
- **FLAV.NUSCALE.04** [exakt, Teil bedingt oder offen] THE 3/5 CLEBSCH DOOR ENUMERATED AND CLOSED (2026-07-15; executes the mechanism search behind FLAV.NUSCALE.03's named escape): v482 recorded that only a third-generation Majorana structure … · `v488_majorana_clebsch_door.py`
- **FLAV.NUSCALE.05** [bedingt] THE SCALARON-TRACE RUNG CANDIDATE (registered 2026-08-28, external master-object review wave 2): M_R = N_fam * M_scal = 3 c3^{7/2} Mbar = 9.180e13 GeV -- a STRUCTURED point inside the v481 … · `v986_nu_scalaron_rung.py`
- **FLAV.NUSCALE.06** [bedingt] THE PENTAGON-CLASS MISALIGNMENT CANDIDATE (registered 2026-08-29, review wave 4, v1001). Candidate [C]/[N]; mechanism = the Q_+-to-flavor operator [O]; NO closure; FLAV.NUSCALE.05 unmoved. … · `v1001_nu_pentagon_chain.py`
- **FLAV.NU.TEXTURE.MECHANISM.01** [offen] THE NEUTRINO-TEXTURE MECHANISM CONTRACT [O] (registered 2026-08-29, review wave 5, v1003). Both structural nulls are BINDING CONSTRAINTS. [E] charged-sector data h=(1,1,1), d=(64,60,64) … · `v1003_nu_structural_nulls.py`
- **FLAV.NU.PARENT.INTERTWINER.01** [exakt, Teil bedingt oder offen] FINITE PARENT-INTERNAL A3-D4 INTERTWINER (registered 2026-09-05, v1024; no parent marker move). A positive, fermion-even, gauge-compatible, non-additive seam/neutrino-pair term is … · `v1024_parent_internal_boundaries.py`

### DIAMOND (13)

- **DIAMOND.AXIS.01** [exakt] DIAMOND AXIS CURVATURE LEMMA (2026-06-15): along the two canonical axes of the centered cross (C center, U winding, V sheet, v95) the Sheet-Diamond determinant is LINEAR on the winding axis … · `v218_diamond_axis_geometry.py`
- **DIAMOND.PLUCKER.01** [exakt] PLUCKER TRANSFER LADDER (2026-06-15): the anchor-Plucker coordinates lift K->C->F in two exact steps: Pl(K)=(-1,6,4), Pl(C)=(0,14,14), Pl(F)=(1,22,30); Pl(C)-Pl(K)=(1,8,10)=(N_Phi, rank … · `v218_diamond_axis_geometry.py`
- **DIAMOND.SPECTRAL.01** [exakt] SPECTRAL RAMIFICATION LADDER (2026-06-15): for Q,K,C,F the characteristic polynomial splits as (t-r)(t^2-bt+c) with r rational, and the full cubic discriminant factors as Disc(chi)=q(r)^2 … · `v218_diamond_axis_geometry.py`
- **DIAMOND.SHEET.SOURCE.J.01** [exakt, Teil bedingt oder offen] The sheet source corner J (2026-06-24): the determinant surface has two s-independent iso-volume walls (v135) t=-1 (det=4=|mu4|, carries K) and t=-2 (det=2=|Z2|); the operator on the Z2 … · `v412_sheet_source_corner_J.py`
- **DIAMOND.SHEET.CALCULUS.01** [exakt, Teil bedingt oder offen] The sheet axis is a discrete characteristic-difference calculator (2026-06-24): on M_t=C+tV (v95/v218) the elementary-symmetric coefficients encode atoms as DIFFERENCE ORDERS. [E] … · `v413_sheet_characteristic_calculus.py`
- **DIAMOND.CENTER.RESOLVENT.01** [exakt, Teil bedingt oder offen] The center C is a resolvent portal G2 -> F4 -> E8 (2026-06-24): the center C=M(1,0) (v95) is the t=0 evaluation of V; its shifted determinants climb the exceptional ladder. [E] CENTER … · `v414_center_resolvent_portal.py`
- **DIAMOND.GAUSS.01** [exakt, Teil bedingt oder offen] The square CM-norm dictionary is operator-realised and the mu4 deck is intrinsic (2026-06-25): the Gaussian arithmetic of v222/v230 lifted from numbers to OPERATORS via J=[U,V]/3, the … · `v415_gaussian_operator.py`
- **DIAMOND.EISEN.01** [exakt, Teil bedingt oder offen] The Eisenstein/CP operator -- the hexagonal dual of the Gaussian seam operator v415 (2026-06-25): where v415 put the SQUARE seam (tau=i, mu4) on operators via J=[U,V]/3, this puts the … · `v417_eisenstein_cp_operator.py`
- **DIAMOND.CENTER.QUOTIENT.01** [exakt] THE CENTER QUOTIENT COMPILER (2026-07-23; external review round, all identities machine-confirmed before promotion; 14 checks, <1 s, sympy exact): the Sheet-Diamond center C = M(1,0) (v95) … · `v530_center_quotient_compiler.py`
- **DIAMOND.PARA.SELFCODE.01** [exakt] THE PARABOLIC ANCHOR SELF-CODE (2026-07-31; external structural review audited by discovery probe parabolic_selfcode_audit_probe.py 23/23, verdict SELFCODE-VERIFIED, and promoted; one … · `v566_parabolic_anchor_selfcode.py`
- **DIAMOND.BIT.SELECT.01** [exakt] THE RESIDUAL RECONSTRUCTION BIT: LOCATED, SEPARATED, AND GIVEN A CANONICAL EXACT WITNESS (2026-07-31; Priorities 2+3 of the self-code review fused with the v567 residual bit; one module, 11 … · `v568_bit_selector_ladder.py`
- **DIAMOND.BIT.RPFORM.01** [exakt, Teil bedingt oder offen] The RP form on the 81-twins -- spectrally BLIND, integrally ANCHOR-VISIBLE (one module, 11 checks, ~1 s; discovery probe rp_bit_form_probe.py 11/11, verdict ANCHOR-VISIBLE): with the corpus … · `v572_rp_bit_form.py`
- **DIAMOND.BIT.EXTSOURCE.01** [exakt, Teil bedingt oder offen] The interacting follow-up of v572, executed -- RANK-ONE EXT-SOURCE COUPLING ON THE GIBBS PENCIL SEPARATES THE 81-TWINS DIRECTIONALLY INSIDE THE UNIT COUPLING BOX (one module, 8 checks, ~1 … · `v584_ext_source_response.py`

### CP (5)

- **CP.MOD.01** [exakt, Teil bedingt oder offen] CP HEXAGONAL MODULUS (2026-06-16): CP lives in the HEXAGONAL phase fiber, not the square seam deck (a geometric sharpening of red-team Target D, the magnitude bijection misses the CP … · `v220_cp_hexagonal_modulus.py`
- **CP.MU6.01** [exakt, Teil bedingt oder offen] CP PHASES AS MU6 POWERS (2026-06-17): both CP phases are mu6 powers of ONE hexagonal CM unit, split by the sheet -- a structural reduction of red-team Target D, sharper than v220/v225 … · `v231_cp_mu6_phases.py`
- **CP.TRIALITY.01** [exakt, Teil bedingt oder offen] CP = THE UNIVERSAL TRIALITY PHASE, SHEET-SPLIT (2026-06-17): one step past CP.MU6.01/v231 in the canonical mu6 = mu3 (family/triality) x mu2 (sheet) factorisation. (1) mu6=mu3 x mu2 … · `v233_cp_triality_phase.py`
- **CP.CHANNEL.IDENT.01** [exakt] THE CP CHANNEL / INVARIANT IDENTITIES (2026-07-29; one module, 25 checks, ~1 s; promotion of the exact cores of discovery T177 (contract CP.INVARIANT, probe cp_invariant_probe.py, 49/49, … · `v561_cp_channel_identities.py`
- **CP.FIBER.PIN.01** [exakt] THE FRAME DOOR OF THE CP MIRROR, DECIDED (2026-07-30; T178, contract FRAME.CARRIER, probe cp_frame_carrier_probe.py 31/31, verdict FIBER-PINNED; one module, 22 checks, ~2 s, sympy-exact … · `v564_cp_fiber_pin.py`

### SEED (4)

- **SEED.ANOMALY.01** [exakt] TRACE-ANOMALY ORIGIN OF THE SEED PREFACTOR 4/3 + the shared integer 7 (b3 = scalaron exponent): (1) the retained seed phi0=(4/3)c3+48 c3^4 has base (4/3)c3=1/(6pi); the prefactor 4/3 is the … · `v172_trace_anomaly_seed.py`
- **SEED.CROSSVAL.01** [exakt] LEAVE-ONE-OUT SEED CROSS-VALIDATION (2026-06-21): turns the 'many observables, one seed phi0' liability into a falsifiable out-of-sample test [N] (Python-only, external data). [N] … · `v306_seed_crossval.py`
- **SEED.CROSSSECTOR.01** [exakt] CROSS-SECTOR ONE-PARAMETER SEED TEST (2026-07-01, v465): the theta13-independent, externally-runnable SUBSET of SEED.CROSSVAL.01/v306 -- the sharpest shared-origin signature is a … · `v465_seed_crosssector_joint.py`
- **SEED.LEPTONMASS.01** [exakt, Teil bedingt oder offen] SIXTH SEED CHANNEL, NEW MEASUREMENT SECTOR (2026-07-01, v466): the charged-lepton mass ratio m_e/m_mu=(12/7)phi0^2 is a sixth independent readout of the ONE axiom seed phi0=1/(6pi)+48 c3^4, … · `v466_seed_leptonmass_channel.py`

### FLAVOR (2)

- **FLAVOR.GRAPH.FILTRATION.01** [exakt] FLAVOR.GRAPH.FILTRATION.01: the flavor-matrix filtration Q -> K -> L as an incidence-compiler fingerprint -- budgets exact, entries NOT reconstructible, one probe (16/16 checks incl. the 4 … · `v812_flavor_graph_filtration.py`
- **FLAVOR.FEEDBACK.NORMALFORM.01** [exakt] THE FEEDBACK NORMAL FORM OF THE WINDING LINE (2026-08-08 round 36, v874; discovery probe flavor_feedback_probe.py 14/14, verdict FEEDBACK-NORMALFORM-EXACT; exact sympy/integer arithmetic in … · `v874_feedback_decoder_lean.py`

### LEP (1)

- **LEP.FROB.01** [exakt, Teil bedingt oder offen] LEPTON ETALE FROBENIUS ALGEBRA (2026-06-16): the three charged-lepton coefficients form a genuine commutative Frobenius algebra over the mu6 resolvent. Exact coefficients … · `v229_lepton_frobenius_algebra.py`

### THETA13 (1)

- **THETA13.BUDGET.01** [exakt] THE THETA13 SUB-LEADING BUDGET vs the current global fit (2026-06-22): the +2 sigma pressure point (v328) quantified against BOTH NuFIT 6.0 normal-ordering variants and JUNO; neither … · `v338_theta13_budget.py`

## Twistor und Himmelskugel

### CELEST (23)

- **CELEST.WP1.01** [exakt, Teil bedingt oder offen] WP1 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-20; CELEST.SEAM.01 work package 1, verdict B): the E8 mu4-glue grading (v128) IS the flat Z4 monodromy of the glue-equivariant … · `v492_celestial_z4_orbifold.py`
- **CELEST.SEAM.01** [offen] THE CELESTIAL-HOLOGRAPHIC RESEARCH CONTRACT (2026-07-20; the fourth research contract alongside CONTRACT.U.01, CONTRACT.G.01 and CONTRACT.F.01): HYPOTHESIS (as corrected by WP1) -- the seam … · `v492_celestial_z4_orbifold.py`
- **CELEST.WP2.01** [exakt, Teil bedingt oder offen] WP2 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-20; CELEST.SEAM.01 work package 2, verdict B; companion of CELEST.WP1.01/v492): the clock-invariant deformation XY = Z^4 + a0 of … · `v493_celestial_wp2_clock_deformation.py`
- **CELEST.WP3.01** [exakt, Teil bedingt oder offen] WP3 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 3, verdict B; the pre-registered 'alignment test only'): the Green-Schwarz axion coefficient … · `v495_celestial_wp3_gs_alignment.py`
- **CELEST.WP4.01** [exakt, Teil bedingt oder offen] WP4 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 4, verdict B(ii); the 'type mismatch, head-on'): the (E8)_1 character E4/eta^8 against the … · `v496_celestial_wp4_character_shadow.py`
- **CELEST.WP5A.01** [exakt, Teil bedingt oder offen] WP5a OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5a, the first constructive milestone of WP5 -- 'the null ideal from the limit'; the v496 … · `v497_celestial_wp5a_null_ideal.py`
- **CELEST.WP5B.01** [exakt] WP5b OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5b, the second constructive milestone of WP5 -- 'the singular vector as an operator'; WP5a/v497 … · `v498_celestial_wp5b_singular_vector.py`
- **CELEST.WP5C.01** [exakt] WP5c OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5c, the third constructive milestone of WP5 -- 'the GNS limit state'; WP5a/v497 made the … · `v500_celestial_wp5c_gns_limit_state.py`
- **CELEST.WP5D.01** [exakt] WP5d-ALPHA OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5d, alpha stage -- 'the two-interval index from the lattice'; the local-net question): … · `v501_celestial_wp5d_two_interval_index.py`
- **CELEST.WP5E.ALPHA.01** [exakt] WP5e-ALPHA OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5e, alpha stage -- 'prefactor and level bookkeeping', the CFT-side anchor of the … · `v502_celestial_wp5e_prefactor_level.py`
- **CELEST.WP5DB.01** [exakt] WP5d-BETA OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5d, beta stage -- 'split + strong additivity from the lattice'; the two remaining … · `v504_celestial_wp5d_klm_completeness.py`
- **CELEST.WP5E.BETA.01** [exakt] WP5e-BETA OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5e, beta stage -- 'the equivariant anomaly ledger on twistor space', the twistor-side … · `v505_celestial_wp5e_beta_equivariant_ledger.py`
- **CELEST.WP5E.GAMMA.01** [exakt] WP5e-GAMMA OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5e, gamma stage -- 'the sphere-axion pairing check', an honest rigid NEGATIVE result in … · `v508_celestial_wp5e_gamma_pairing.py`
- **CELEST.WP5E.EPS2.01** [exakt] WP5e-EPSILON-2 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5e, epsilon-2 stage -- 'the CPS level-from-flux dial', verdict B; exact sympy/Fraction … · `v509_celestial_wp5e_eps2_level_from_flux.py`
- **CELEST.WP5E.DELTA2.01** [exakt] WP5e-DELTA-2 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-21; CELEST.SEAM.01 work package 5e, delta-2 stage -- 'the full-tensor ledger', the named escape route of the v508 … · `v511_celestial_wp5e_delta2_tensor_ledger.py`
- **CELEST.DTERM.NONDERIV.01** [exakt, Teil bedingt oder offen] THE c_d NEGATIVE CERTIFICATE -- CAN GEOMETRY DERIVE c_d = 1920? (2026-07-22; anti-numerology certificate in the v508 style, honest rigid NEGATIVE result; 24 checks, ~2 s; NO marker moves -- … · `v513_celestial_dterm_nonderivation.py`
- **CELEST.WP5E.EPS1.01** [exakt] WP5e-EPSILON-1 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-22; CELEST.SEAM.01 work package 5e, epsilon-1 stage -- 'the O(-2) bulk-axion slot', verdict B; exact sympy/Fraction … · `v514_celestial_wp5e_eps1_axion_slot.py`
- **CELEST.WP5E.M1.01** [exakt] MILESTONE M1 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-22; CELEST.SEAM.01 WP5e back-reaction programme, milestone M1 -- 'the A3 Omega_N', SUCCESS on the preregistered v514 … · `v515_celestial_wp5e_m1_omega_n.py`
- **CELEST.WP5E.M2.01** [exakt] MILESTONE M2 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-22; CELEST.SEAM.01 WP5e back-reaction programme, milestone M2 -- 'the twisted Kodaira-Spencer measure', SUCCESS on the … · `v516_celestial_wp5e_m2_twisted_ks_measure.py`
- **CELEST.WP5E.M3.01** [exakt] MILESTONE M3 OF THE CELESTIAL-HOLOGRAPHIC CONTRACT EXECUTED (2026-07-22; CELEST.SEAM.01 WP5e back-reaction programme, milestone M3 -- 'the a0 uplift', SUCCESS on the preregistered v514 S8.3 … · `v517_celestial_wp5e_m3_a0_uplift.py`
- **CELEST.WP5E.DELTA1.01** [exakt] THE DELTA-1 CHAIN OF THE CELESTIAL-HOLOGRAPHIC CONTRACT DECIDED (2026-07-22; CELEST.SEAM.01 WP5e, delta-1b + delta-1c + delta-1d consolidated -- 'the derived chiral measure kills the … · `v518_celestial_wp5e_delta1_derived_kill.py`
- **CELEST.WP5E.MEASURE.01** [exakt] THE WP5E BCOV MEASURE QUESTION DECIDED AT PROBE LEVEL (2026-07-22; delta-1e + delta-1f consolidated in the v518 style, one module, 39 checks, ~20 s; ERFOLG-A on the preregistered decision … · `v520_celestial_wp5e_measure_decision.py`
- **CELEST.WP5E.WM.01** [exakt] THE CONSTRUCTIVE w_m DERIVATION EXECUTED (2026-07-23; the residual [O] of v516/v520, one module, 26 checks, ~1 s; verdict ERFOLG per the frozen preregistration): is the completion … · `v523_celestial_wp5e_wm_derived.py`

### WOIT (9)

- **WOIT.OS.TWISTOR.01** [offen] THE OSTERWALDER-SCHRADER TWISTOR BRIDGE (2026-07-22; the new central research contract after the external review of the celestial/twistor branch -- THIS IS THE ACTUAL BRIDGE FROM COMPILER … · `v519_woit_theta_rp_free.py`
- **WOIT.THETA.FREE.01** [exakt] THE WOIT-ALPHA MILESTONE (2026-07-22): the real structure exists, and free reflection positivity picks the SAME family -- the alpha stage of WOIT.OS.TWISTOR.01 executed as one module (Part … · `v519_woit_theta_rp_free.py`
- **WOIT.BETA1.GSO.01** [bedingt] THE WOIT-BETA1 MILESTONE EXECUTED (2026-07-22; one module, 16 checks, ~5 s; typed result, verdict UNDECIDED per the frozen preregistration -- 'the clock is time-like, GSO is the gauge … · `v522_woit_beta1_gso_gauge.py`
- **WOIT.BETA2.OS.01** [bedingt] THE WOIT-BETA2 MILESTONE EXECUTED (2026-07-23; one module, 20 checks, ~5 s; verdict SUCCESS per the frozen preregistration, [C]-typed per contract precision (iii)): contract clause beta_2 … · `v524_woit_beta2_os_quotient.py`
- **WOIT.BETA3.DUALITY.01** [exakt] THE WOIT-BETA3 MILESTONE EXECUTED (2026-07-30; one module, 22 checks, ~2 s; probe woit_os_beta3_pt_dual_probe.py 25/25, verdict SUCCESS per the frozen preregistration): 'the PT <-> PT* … · `v565_woit_beta3_pt_duality.py`
- **WOIT.GAMMA.FREE.01** [exakt, Teil bedingt oder offen] The free gamma slice of the WOIT roadmap (one module, 9 checks, ~1 s; discovery probe gamma_free_slice_probe.py 9/9, verdict GAMMA-FREE-SLICE-LANDED; the v524 machinery imported read-only). … · `v608_gamma_free_slice.py`
- **WOIT.GAMMA.TOY.01** [exakt, Teil bedingt oder offen] The interacting gamma slice at toy level (one module, 11 checks, ~2 s; discovery probe gamma_toy_interacting_probe.py 11/11, verdict GAMMA-TOY-LANDED). The gamma battery executed ON the … · `v615_gamma_toy_interacting.py`
- **WOIT.SU2.KIN.01** [exakt, Teil bedingt oder offen] The internal SU(2), kinematic separation at twistor level (one module, 17 checks, <1 s; discovery probe su2_internal_kinematic_probe.py 17/17, verdict SU2-KINEMATIC-SEPARATED; the first … · `v616_su2_internal_kinematic.py`
- **WOIT.GAMMA.SEMI.01** [exakt, Teil bedingt oder offen] The interacting Klein-Landau semigroup -- the clock survives the interaction as a positive transfer exactly on the mu4-commensurate steps (one module, 6 checks, ~2 s; discovery probe … · `v621_interacting_semigroup.py`

## Gravitation, Horizont und Geometrie

### QGEO (76)

- **QGEO.RIGID.01** [exakt] THE QGEO RIGIDITY THEOREM (the FINITE half of the Tier-1 hinge GATE.QGEO, hardened to [E]; v167 Tier 1): a genus-0 boundary with a faithful D4 action, one orbit of four parabolic marks, H^1 … · `v168_qgeo_rigidity.py`
- **QGEO.REALIZE.01** [bedingt] THE SEAM-COLLAR REALISATION (the physical half of GATE.QGEO; the one Tier-1 realisation premise, stays [C]/[O]): the QGEO Rigidity Theorem (QGEO.RIGID.01/v168) closes the PHYSICS only if … · `v168_qgeo_rigidity.py`
- **QGEO.REALIZE.02** [exakt, Teil bedingt oder offen] THE SEAM COLLAR REALISATION THEOREM, ASSEMBLED AS ONE CENTRAL STATEMENT + REDUCTION CERTIFICATE (the single structural proof obligation of the whole theory, formulated cleanly; NOT a … · `v176_seam_collar_realisation.py`
- **QGEO.TREE.01** [exakt, Teil bedingt oder offen] THE QGEO PROOF TREE (the v176 'one realisation premise' split into its honest constituents; the open residual sharpened from one diffuse premise into exactly TWO named obligations; external … · `v177_seam_marking_kernel.py`
- **QGEO.MARKS.01** [bedingt] THE SEAM MARKING OBLIGATION (genuine open obligation 1 of 2; the structural residual of the theory, NOT a finite computation): the RAW RP seam collar must canonically produce a genus-0, … · `v177_seam_marking_kernel.py`
- **QGEO.KERNEL.01** [bedingt] THE SEAM CALDERON-KERNEL OBLIGATION (genuine open obligation 2 of 2; the structural residual of the theory, NOT a finite computation): the RAW RP seam Calderon operator must equal the … · `v177_seam_marking_kernel.py`
- **QGEO.REDUCE.01** [exakt, Teil bedingt oder offen] DEEPER REDUCTION OF THE TWO OBLIGATIONS (an honest attempt at QGEO.MARKS.01 + QGEO.KERNEL.01; NEITHER is closed -- they are genuine constructive-geometry/AQFT statements, not finite … · `v178_marks_kernel_attempt.py`
- **QGEO.CONF.01** [exakt, Teil bedingt oder offen] THE SEAM CONFORMAL-REALISATION PREMISE -- the SINGLE structural residual of the whole theory (the two obligations QGEO.MARKS.01 + QGEO.KERNEL.01 collapse onto it; honest unification, the … · `v179_conformal_realisation.py`
- **QGEO.ISO.01** [exakt, Teil bedingt oder offen] THE SEAM-ISOMETRY PREMISE -- the milder, non-circular form of the single structural residual (QGEO.CONF.01 is REPLACED by it via three classical theorems; honest reduction, still OPEN): … · `v180_clock_is_mobius.py`
- **QGEO.SYM.01** [exakt, Teil bedingt oder offen] A COROLLARY OF SEAM.EQUIV.01 (no independent open premise; reclassified 2026-06-22, v335): per v323/QGEO.BW.01 it is downstream of SEAM.EQUIV.01 + a rotation-invariant vacuum, and the … · `v181_clock_is_conformal_symmetry.py`
- **QGEO.ENERGY.01** [exakt, Teil bedingt oder offen] THE ENERGY-CONSERVING-CLOCK REFORMULATION (external-review proposal 2026-06-14; HONESTLY typed: RELOCATES the bedrock, does NOT eliminate it): the proposed premise QGEO.ENERGY.01 (rho^* … · `v192_qgeo_energy_clock.py`
- **QGEO.ENERGY.02** [offen] THE ENERGY-COMMUTATOR CONTRACT (external-review proposal 2026-06-14; the NON-CIRCULAR sharpening of QGEO.ENERGY.01): the proposed premise is the operator identity [rho, Lambda_Sigma]=0 on … · `v193_qgeo_energy_commutator.py`
- **QGEO.DTN.01** [exakt, Teil bedingt oder offen] SUBCLAIM 2 OF QGEO.ENERGY.02 TACKLED -- 'the raw RP-seam DtN is RP-definable' (external-review next-lever 2026-06-14): a REDUCTION, not a closure. (1) FINITE CORE [E]: rho:z|->iz on … · `v194_raw_seam_dtn_rp.py`
- **QGEO.MARKS.02** [exakt, Teil bedingt oder offen] LEFSCHETZ/CHARACTER DERIVATION OF THE FOUR SEAM MARKS (external-review path 2/priority 2, 2026-06-14): instead of POSITING D4 boundary marks (the old QGEO.MARKS.01 finite core), the marks … · `v195_marks_lefschetz_character.py`
- **QGEO.VARI.01** [exakt, Teil bedingt oder offen] THE E_fail VARIATIONAL FUNCTIONAL (external-review path 1, 2026-06-14): a concrete, numerically-testable sharpening of the v194/QGEO.ENERGY.02 target. E_fail(Sigma,rho) = … · `v196_seam_energy_functional.py`
- **QGEO.MODULAR.01** [exakt, Teil bedingt oder offen] CRACKING THE BEDROCK TO A STATE-INVARIANCE (2026-06-14): the v194 residual -- the full-L^2 commutator [rho,Lambda_Sigma]=0, framed as Bisognano-Wichmann geometric modular flow with a … · `v198_modular_commutator_reduction.py`
- **QGEO.STATE.01** [exakt, Teil bedingt oder offen] ATTACKING omega o rho = omega DIRECTLY ON THE RAW SEAM (work package A, external-review priority 1, 2026-06-14): a further localising reduction, not a closure. (1) SETUP NON-CIRCULAR [E]: … · `v199_seam_state_invariance.py`
- **QGEO.SUBPRIN.01** [exakt, Teil bedingt oder offen] THE v199 BEDROCK RESIDUAL REDUCED ONE GENUINE STEP FURTHER (2026-06-15): block-diagonality of the sub-principal symbol is NOT an independent postulate. (1) [E] the seam DtN is … · `v201_seam_subprincipal_marks.py`
- **QGEO.PILLOW.01** [exakt, Teil bedingt oder offen] PILLOWCASE REDUCTION OF QGEO.SYM.01 (2026-06-15): the two separately-tracked open QGEO residuals -- QGEO.ISO.01 (v180; the carrier clock is an order-4 ISOMETRY of the seam metric) and … · `v214_seam_pillowcase.py`
- **QGEO.KILL.01** [offen] BEDROCK AS A FALSIFIABLE KILL-TEST (2026-06-15): QGEO.SYM.01 (omega o rho = omega) is recast as four independently-falsifiable predictions about the raw seam DtN, NOT a proof (a … · `v215_seam_deck_killtest.py`
- **QGEO.MARKS.03** [exakt, Teil bedingt oder offen] THE FOUR MARKS EMERGE FROM GAUSS-BONNET (2026-06-15; executes the K1 lever of QGEO.KILL.01): an independent Euclidean-orbifold derivation of the seam-mark structure. (1) [E] MARK COUNT = 2 … · `v216_marks_gauss_bonnet.py`
- **QGEO.EMERGE.01** [exakt, Teil bedingt oder offen] NUMERICAL EMERGENCE SCAN (2026-06-15; the numerical sibling of QGEO.MARKS.03, as v210 is to v201): with the NUMBER of marks n AND their POSITIONS free, the 4-equally-spaced (square) … · `v217_marks_emergence_scan.py`
- **QGEO.VARI.02** [exakt, Teil bedingt oder offen] THE DISCRETISED VARIATIONAL TEST (work package B, external-review, 2026-06-14): v196 showed E_fail=0 at the mu4 config on the fixed H^1 block; here the clock angle is a FREE parameter and … · `v200_seam_variational_scan.py`
- **QGEO.MARKLOCAL.01** [exakt, Teil bedingt oder offen] MARK-LOCAL DtN CERTIFIED ON REALISTIC STEKLOV PROFILES + AT THE STATE LEVEL (2026-06-15): a numerical SHARPENING of QGEO.SUBPRIN.01/v201, carried to the quasi-free STATE (the actual omega o … · `v210_mark_local_dtn.py`
- **QGEO.RR.01** [exakt, Teil bedingt oder offen] RIEMANN-ROCH INDEX GATE (2026-06-16): P2 as the mode space of a degree-4 seam divisor (bedrock LANGUAGE for P2, NOT a proof from nothing). deg D=4=|mu4| on P^1 => h0(O(D))=deg+1=5=g_car … · `v228_rr_index_gate.py`
- **QGEO.ENERGY.03** [offen] RAW-SEAM => MARK-LOCALITY REDUCTION (2026-06-18): composes v216+v201+v198 on the canonical pillowcase metric to make the converse implication explicit and isolate the bedrock. The flat … · `v264_raw_seam_marklocal.py`
- **QGEO.SYM.02** [exakt, Teil bedingt oder offen] QGEO RIGIDITY / MINIMAL-AXIOM FORM (2026-06-18): the bedrock QGEO.SYM.01 is sharpened to its minimal form -- the single bare statement 'the raw RP seam admits the carrier mu4 clock as a … · `v267_qgeo_rigidity_minimal_axiom.py`
- **QGEO.SYM.03** [exakt, Teil bedingt oder offen] FLAT PILLOWCASE CLOSES THE MODULAR COMMUTATOR TO ALL ORDERS (2026-06-18): collapses the bedrock QGEO.SYM.01 to ONE geometric postulate (does NOT prove the premise). For a quasi-free seam … · `v276_qgeo_flat_closes_commutator.py`
- **QGEO.OBLIG.01** [exakt, Teil bedingt oder offen] QGEO.SYM.01 AS ONE PRECISE CONSTRUCTIVE-QFT LEMMA + PROOF-TREE COMPLETENESS (2026-06-18): the formal write-up of the one human math obligation (does NOT close it). LEMMA (constructive-QFT … · `v279_qgeo_obligation_lemma.py`
- **QGEO.STEKLOV.01** [exakt, Teil bedingt oder offen] FLAT tau=i PILLOWCASE DtN EXPERIMENT (2026-06-18): a direct numerical experiment on the actual seam orbifold, realising the QGEO chain on the geometric side (does NOT prove the premise). … · `v280_pillowcase_steklov.py`
- **QGEO.ROUTEI.01** [exakt, Teil bedingt oder offen] ROUTE (i) (RP-STATE UNIQUENESS) LEMMA DECOMPOSITION (2026-06-18): turns Route (i) from one hard target into a 6-lemma chain with 5 lemmas discharged [E]/[F] and exactly one open lemma … · `v284_route_i_rp_uniqueness.py`
- **QGEO.MODHAM.01** [exakt, Teil bedingt oder offen] BEDROCK omega o rho = omega VIA AN EXPLICIT MODULAR HAMILTONIAN (2026-06-21, TOE attack point 2): the next concrete step on the open bedrock QGEO.SYM.01. [E] NON-CIRCULAR upstream: marks=2 … · `v309_modular_bedrock.py`
- **QGEO.BW.01** [exakt, Teil bedingt oder offen] THE KEYSTONE PUSHED ONE STEP -- BISOGNANO-WICHMANN MAKES THE SEAM MODULAR FLOW GEOMETRIC (2026-06-22): the mu4 deck invariance is NOT an independent axiom. [E] CLOCK = GEOMETRIC ROTATION: … · `v323_bw_geometric_modular.py`
- **QGEO.KEYSTONE.01** [Regel oder Synthese] THE SEAM KEYSTONE AS ONE CITABLE THEOREM (2026-06-22): the raw seam state is the flat tau=i pillowcase state, with EXACTLY ONE open lemma. Assembles the keystone pieces into one statement … · `v325_pillowcase_keystone.py`
- **QGEO.OSGAP.01** [exakt, Teil bedingt oder offen] CLOSING THE OS STEP: THE BULK GAP IS THE TRANSFER GAP (2026-06-22, second quantization) + the single sufficient premise. The v308 chain left ONE analytic input -- the OS identification … · `v329_os_gap_reduction.py`
- **QGEO.NECESS.01** [exakt, Teil bedingt oder offen] THE NECESSITY OF H REDUCED (2026-06-22): rotation-invariance of the seam DtN is EQUIVALENT to the four marks sitting at the square (tau=i), so the one open lemma is exactly the order-4 CM … · `v331_necessity_of_H.py`
- **QGEO.CMSELECT.01** [exakt, Teil bedingt oder offen] THE CM-SELECTION MECHANISM (2026-06-22): tau=i (the square) is the UNIQUE order-4 modulus and the attractor of mark-equilibration, so the last keystone residual is 'the deck acts … · `v333_cm_selection.py`
- **QGEO.MULTILOCAL.01** [exakt, Teil bedingt oder offen] THE FOUR mu4 MARKS ARE A FOUR-INTERVAL FREE-FERMION GEOMETRY (2026-07-14, bird's-eye round): the seam circle minus the four Gauss-Bonnet marks (v216) IS the n = 4 disjoint-interval region … · `v480_multilocal_four_interval.py`
- **QGEO.EMERGE.LIGHT.01** [exakt] THE v215 DEFERRED EMERGENCE TEST (LEVER 7) EXECUTED LIGHT + RESIDUAL CONVERGENCE QGEO-R1 == P2-R1 (2026-07-21; analysis/kill-test class, v347/v264 style; NO marker moves -- QGEO.REALIZE.01 … · `v503_qgeo_emergence_light.py`
- **QGEO.QRECON.01** [exakt] THE NO-CARRIER-RECONSTRUCTION AUDIT, STAGE A (2026-07-31; Priority 1 of the self-code review; one module, 9 checks, ~1 s, exact enumeration with integer prefilters + sympy verification; … · `v567_carrier_free_q_reconstruction.py`
- **QGEO.WINDING.SELF.01** [exakt, Teil bedingt oder offen] Stage B of the No-Carrier-Reconstruction audit, first narrowing (one module, 11 checks, ~6 s; discovery probe stageb_winding_probe.py 12/12, verdict FORCED): the involution TYPE is killed … · `v574_stageb_winding.py`
- **QGEO.INVOL.01** [exakt, Teil bedingt oder offen] The sheet-involution existence question of GATE.QGEO.01, narrowed to a two-demand selection with UNCONDITIONAL existence (one module, 7 checks, <1 s; discovery probe … · `v590_involution_existence.py`
- **QGEO.COVER.01** [exakt, Teil bedingt oder offen] Strategy I/II, first executed step -- the mu3-cover of P^1 branched at mu4 (the curve y^3 = x^4 - 1) carries FOUR independent compiler fingerprints (one module, 6 checks, <1 s; discovery … · `v597_mu4_cover_model.py`
- **QGEO.DICT.01** [exakt] The cover dictionary, round 2 -- THREE compiler objects located as canonical cusp-Gram combinations on the mu3-cover (one module, 4 checks, <1 s; discovery probe cover_dictionary_probe.py … · `v598_cover_dictionary.py`
- **QGEO.REAL.01** [exakt, Teil bedingt oder offen] The real structure on the mu3-cover -- the v598 suspect constructed (one module, 18 checks, ~3 s; discovery probe real_structure_probe.py 18/18, verdict REAL-STRUCTURE-LANDED). EXISTENCE … · `v599_real_structure.py`
- **QGEO.EMBED.01** [exakt, Teil bedingt oder offen] The joint embedding -- the compiler pair (U, V) is realized on the mu3-cover (one module, 15 checks, ~4 s; discovery probe joint_embedding_probe.py 15/15, verdict JOINT-EMBEDDING-LANDED). … · `v600_joint_embedding.py`
- **QGEO.CANON.01** [exakt, Teil bedingt oder offen] The canonical equivariant realization (one module, 14 checks, ~12 s; discovery probe equivariant_dual_probe.py 14/14, verdict EQUIVARIANT-DUAL-LANDED). THE REPAIR: the naive cusp-Gram … · `v601_equivariant_dual.py`
- **QGEO.DUALFORMS.01** [exakt, Teil bedingt oder offen] The duality-and-forms round (one module, 17 checks, ~8 s; discovery probe duality_forms_probe.py 17/17, verdict DUALITY-FORMS-LANDED). THE DUALITY MECHANISM [E]: U_1 = 3 GJ_1 is … · `v602_duality_forms.py`
- **QGEO.SEAMMARKS.01** [exakt, Teil bedingt oder offen] The seam-and-marks round (one module, 15 checks, ~7 s; discovery probe seam_marks_probe.py 15/15, verdict SEAM-REAL-LANDED). THE SEAM LINE IS REAL [E]: the line-stabilizer pair (U_1, V_1^J) … · `v603_seam_marks.py`
- **QGEO.EQORDER.01** [exakt, Teil bedingt oder offen] The equivariant order and the infinity bridge (one module, 12 checks, ~4 s; discovery probe equivariant_order_probe.py 12/12, verdict EQUIVARIANT-ORDER-LANDED). CORRECTION + THE INTEGRAL … · `v604_equivariant_order.py`
- **QGEO.TRANSLATOR.01** [exakt, Teil bedingt oder offen] The translator round (one module, 11 checks, ~11 s; discovery probe translator_probe.py 11/11, verdict TRANSLATOR-LANDED). THE TRANSLATOR EXISTS [E, centerpiece]: pulling the corpus RP … · `v605_translator.py`
- **QGEO.MODESEP.01** [exakt, Teil bedingt oder offen] The mode-separation round (one module, 7 checks, ~2 s; discovery probe mode_separation_probe.py 7/7, verdict MODES-SEPARATED). ALL THREE MODES ARE PURE TORSION AT INFINITY [E]: the three … · `v606_mode_separation.py`
- **QGEO.SELRULE.01** [exakt, Teil bedingt oder offen] The selection rule round (one module, 10 checks, ~2 s; discovery probe selection_rule_probe.py 10/10, verdict SELECTION-RULE-LANDED). THE F_3 STRUCTURE [E]: the torsion of coker(1 - T_inf) … · `v607_selection_rule.py`
- **QGEO.CFACTOR.01** [exakt, Teil bedingt oder offen] The translator factorization (one module, 8 checks, ~8 s; discovery probe translator_factorization_probe.py 8/8, verdict C-FACTORIZATION-LANDED). THE MIRROR INTERTWINER [E]: the … · `v609_translator_factorization.py`
- **QGEO.LANDSCAPE.01** [exakt, Teil bedingt oder offen] The curve landscape round -- what else the mu3-cover carries (one module, 9 checks, ~1 s; discovery probe curve_landscape_probe.py 9/9, verdict LANDSCAPE-MAPPED). THE JACOBIAN SPLITS [E … · `v610_curve_landscape.py`
- **QGEO.PERIODS.01** [exakt, Teil bedingt oder offen] The periods round -- the analytic layer of the cover is classical (one module, 6 checks, ~1 s; discovery probe periods_probe.py 6/6, verdict PERIODS-CLASSICAL). EXACT BETA IDENTITIES [E, … · `v611_periods.py`
- **QGEO.POLFRAME.01** [exakt, Teil bedingt oder offen] The polarization frame round (one module, 9 checks, ~2 s; discovery probe polarization_frame_probe.py 9/9, verdict POLARIZATION-FRAMED). THE ROTATION SPECTRUM [E]: the braid rotation r has … · `v612_polarization_frame.py`
- **QGEO.PERNORM.01** [exakt, Teil bedingt oder offen] The canonical period normalization -- the analytic bridge closes (one module, 19 checks, ~6 s; discovery probe canonical_periods_probe.py 19/19, verdict PERIOD-NORMALIZED; executes v612's … · `v613_canonical_periods.py`
- **QGEO.TRANSPORT.01** [exakt, Teil bedingt oder offen] The Gauss-Manin transport round -- the deck step IS the boundary monodromy (one module, 12 checks, ~7 s; discovery probe transport_deck_probe.py 12/12, verdict TRANSPORT-CANONICAL; closes … · `v614_transport_deck.py`
- **QGEO.SEAMFORCE.01** [exakt, Teil bedingt oder offen] The seam forcing round -- the conformal seam axioms force the mu3-cover y^3 = x^4 - 1 up to the conjugate sheet (one module, 14 checks, <1 s; discovery probe seam_cover_forcing_probe.py … · `v617_seam_cover_forcing.py`
- **QGEO.NCENSUS.01** [exakt, Teil bedingt oder offen] The cyclic-N census -- the seam admits N in {3, 5} as primitive-Lorentz orders and the compiler's integral constants pin N = 3 (one module, 9 checks, <1 s; discovery probe … · `v620_cyclic_n_census.py`
- **QGEO.SEAMID.01** [exakt, Teil bedingt oder offen] The seam identification round -- the physical seam IS the conformal seam (one module, 10 checks, ~1 s; discovery probe seam_identification_probe.py 10/10, verdict SEAM-IDENTIFIED; closes … · `v622_seam_identification.py`
- **QGEO.COVERSEAM.01** [exakt, Teil bedingt oder offen] The covered seam round -- the lattice skeleton of the covering-level identification (one module, 9 checks, ~4 s; discovery probe covered_seam_probe.py 9/9, verdict COVERED-SEAM-SKELETON; … · `v623_covered_seam.py`
- **QGEO.ORBCAS.01** [exakt] The orbifold Casimir round (one module, 5 checks, ~1 s; discovery probe orbifold_casimir_probe.py 5/5, verdict ORBIFOLD-CASIMIR-EXACT). The covered seam's three deck-mode classes {1/6, 1/2, … · `v628_orbifold_casimir.py`
- **QGEO.TWISTOPE.01** [exakt, Teil bedingt oder offen] The twist-field/OPE round (one module, 25 checks, ~1 s; discovery probe orbifold_twist_ope_probe.py 25/25, verdict TWIST-OPE-COHERENT-RP-OPEN). The interacting/CFT-level slice of the Z3 … · `v639_twist_ope.py`
- **QGEO.KLEINRP.01** [exakt, Teil bedingt oder offen] Parafermionic Klein-twist reflection positivity of the Z3 twist sector (one module, 17 checks, ~30 s; discovery probe parafermion_klein_rp_probe.py 17/17, verdict KLEIN-RP-POSITIVE). The … · `v645_klein_rp.py`
- **QGEO.ORBMOD.01** [exakt, Teil bedingt oder offen] The modular data of the Z3 orbifold of the seam on the lattice (one module, 17 checks, ~25 s; discovery probe orbifold_modular_probe.py 17/17, verdict ORB-MODULAR-COHERENT). SECTORS [E]: … · `v650_orbifold_modular.py`
- **QGEO.ORBASM.01** [exakt, Teil bedingt oder offen] The full Z3-orbifold assembly of the seam (one module, 17 checks, ~35 s; discovery probe orbifold_assembly_probe.py 17/17, verdict ORB-ASSEMBLY-SEAM). COHOMOLOGY [E]: H^2(Z3, U(1)) = 0 by … · `v651_orbifold_assembly.py`
- **QGEO.ORBARF.01** [exakt, Teil bedingt oder offen] Arf/GSO hardening of the B-vs-C3 discriminator (one module, 15 checks, ~30 s; discovery probe orbifold_arf_probe.py 15/15, verdict ARF-B-STRUCTURAL). THE LADDER [E]: the mu6 defect ladder … · `v652_orbifold_arf.py`
- **QGEO.BONDDEF.01** [exakt, Teil bedingt oder offen] The bond-defect premise upgraded to a lattice THEOREM (one module, 18 checks, ~40 s; discovery probe bond_defect_premise_probe.py 18/18, verdict BOND-DEFECT-CANONICAL). THE THEOREM [E]: the … · `v653_bond_defect.py`
- **QGEO.ORBOS.01** [exakt, Teil bedingt oder offen] The first continuum-OS slice of the seam orbifold B -- the Osterwalder-Schrader reconstruction data measured at the CONVERGENCE level (one module, 19 checks, ~3 s; discovery probe … · `v679_orbifold_continuum_os.py`
- **QGEO.DECKSECTOR.01** [exakt] TFPT scope of original v705: only D1 NS-circle clock/deck finite algebra, 48 modes with three 16-mode sectors. The prime/digamma application is maintained in Prime Front. · `v705_chain_deck_sector.py`
- **QGEO.CARLIMIT.01** [exakt] QGEO-CAR -- the quasi-free CAR route to the formal orbifold continuum limit (GATE.QGEO), first slice (one module, 22 checks, ~10 s; discovery probe qgeo_car_continuum_probe.py 22/22, … · `v715_qgeo_car_continuum.py`
- **QGEO.CARL1.01** [exakt] The (L1) sector lemma of the CAR continuum majorant, proven modulo elementary steps with explicit constants (one module, 11 checks, ~1 s; discovery probe qgeo_car_l1_sector_lemma_probe.py … · `v745_qgeo_car_l1_sector_lemma.py`
- **QGEO.FANODISC7.01** [exakt] The preregistered Fano -> discriminant(-7) bridge -- the three-cycle half is real, the parity half dies arithmetically: the seven-correlation is typed decorative for good (one module, 26 … · `v757_fano_disc7_bridge.py`
- **QGEO.GRIDSUP.01** [exakt] Remainder (a) of the (L1) majorant lemma FULLY closed in the corrected (1 + ln M) form -- one module from two probes (19 + 19 checks, ~40 s; discovery probes qgeo_grid_supremum_probe.py … · `v778_qgeo_gridsup_logsum.py`

### GRAV (28)

- **GRAV.FR.01** [exakt, Teil bedingt oder offen] R+R^2 covariant field equation CLOSED [I/P]: f(R)=R+R^2/(6M^2), M=c3^(7/2)Mbar=3.06e13 GeV=scalaron mass (m^2=1/(3 f_RR)=M^2); f(R) Einstein-side eq + 4D trace; N*=57 -> … · `v28_gravity_fR.py`
- **GRAV.SCAL.01** [exakt] scalaron exponent 7; M_scal^2/Mbar^2 = c3^7 (3.06e13 GeV) · `v7_gravity_cosmo.py`
- **GRAV.BACKREACT.01** [exakt, Teil bedingt oder offen] THE MATTER-GRAVITY BACKREACTION WITH THE EXPLICIT J3 FLUX (2026-06-23, Direction 5 of next-plan.md). With the parameter-free Einstein equation (v359, G_ab+Lambda g_ab=(1/c3)T_ab) and the … · `v361_grav_backreaction.py`
- **GRAV.GAPCORR.01** [exakt, Teil bedingt oder offen] THE GAP-INDUCED HIGHER-CURVATURE CORRECTION TO GR -- THE DISCIPLINED RESULT (2026-06-23, Direction 2 of next-plan.md). Question: does the seam transfer gap {1,(2/3)^6,(1/3)^6} yield a NEW, … · `v360_grav_gap_corrections.py`
- **GRAV.NONLINEAR.01** [exakt, Teil bedingt oder offen] THE FULL COVARIANT EINSTEIN EQUATION, PARAMETER-FREE (2026-06-23, Direction 1 of next-plan.md, extends GRAV.ENTROPY.EQUILIBRIUM.01/v358). v358 ran the entanglement first law at fixed radius … · `v359_grav_nonlinear_einstein.py`
- **GRAV.ENTROPY.EQUILIBRIUM.01** [exakt, Teil bedingt oder offen] THE ENTANGLEMENT-EQUILIBRIUM DERIVATION OF THE LINEARISED EINSTEIN EQUATION, PARAMETER-FREE (2026-06-23). Carries out the Jacobson-2015 / Faulkner-et-al. 'first law of entanglement => … · `v358_grav_entropy_equilibrium.py`
- **GRAV.SPIN2.UNITARITY.01** [exakt, Teil bedingt oder offen] Track 2 (2026-06-23): the Barnes-Rivers spin decomposition of the graviton propagator -- perturbative graviton unitarity SECTOR-BY-SECTOR. Extends v304 (scalar pole algebra) to the tensor … · `v370_grav_spin2_unitarity.py`
- **GRAV.KMS.HESSIAN.01** [exakt, Teil bedingt oder offen] The KMS Entire Hessian (2026-06-23): the Stelle ghost is EXACTLY the Seeley-DeWitt truncation, and resummation pushes it to infinity -- upgrading v304/v370 from the ASSUMPTION 'the resummed … · `v380_grav_kms_hessian.py`
- **GRAV.XI.01** [exakt, Teil bedingt oder offen] GRAVITATIONAL 3/4 OVER-DETERMINATION (2026-06-15, archive integration of the old Eliminating-K note): the torsion-compression ansatz kappa^2 = xi phi0/c3^2 with the Einstein limit kappa^2=8 … · `v205_xi_threequarter.py`
- **GRAV.H0.01** [exakt, Teil bedingt oder offen] COSMOLOGICAL-CONSTANT BRANCH H0 = 66.5..67.1 (2026-06-15, archive integration of the tfpt-31 CC note; sharpens the qualitative H0-vs-Lambda row): the seam number delta_Sigma = m_1 e^{-2 pi … · `v206_h0_lambda_branch.py`
- **GRAV.ASYMP.01** [offen] ASYMPTOTIC SAFETY -- A SECOND [P] ROUTE TO THE METRIC SECTOR (2026-06-15, archive integration of the old Five-Problems note; NOT load-bearing, does NOT close G_metric): the current theory … · `v207_asymptotic_safety.py`
- **GRAV.AMPLITUDE.01** [exakt, Teil bedingt oder offen] The entire-form-factor graviton-exchange AMPLITUDE is finite, UV-softened and tree-unitary (2026-06-24) -- perturbative gravity as an explicit scattering problem, extending v304/v370/v380 … · `v386_grav_amplitude.py`
- **GRAV.LOOP.FINITE.01** [exakt] The entire-form-factor graviton is UV-finite at the LOOP level by power counting (2026-06-24) -- the honest next step past v386 (which established the tree amplitude), extending … · `v389_grav_loop_finiteness.py`
- **GRAV.NONLOCAL.01** [exakt, Teil bedingt oder offen] THE NONLOCAL SPECTRAL-GRAVITY ACTION; LOCAL R+R^2 IS ITS IR PROJECTION (2026-06-24, Paper C / the 4D gravity-action contract): consolidates v304/v370/v380/v386 and adds the IR-matching … · `v400_grav_nonlocal_action.py`
- **GRAV.ENTROPIC.ACTION.01** [exakt, Teil bedingt oder offen] THE ENTROPIC-ACTION BRIDGE (2026-07-13) -- Bianconi's 'Gravity from entropy' (Phys. Rev. D 111, 066001 (2025); arXiv:2408.14391; DOI 10.1103/PhysRevD.111.066001) quantified as an EXTERNAL … · `v473_entropic_action_bridge.py`
- **GRAV.ENTROPIC.HODGE.01** [exakt, Teil bedingt oder offen] THE OPERATOR LEVEL OF THE ENTROPIC-ACTION BRIDGE (2026-07-13; work packages 1 + 4-algebra of GRAV.ENTROPIC.ACTION.01/v473) -- upgrades the v473 vector-space identity to the operator level … · `v474_entropic_hodge_carrier.py`
- **GRAV.ENTROPIC.SCALARON.01** [exakt, Teil bedingt oder offen] THE R^2 KILL TEST OF THE ENTROPIC-ACTION BRIDGE, EXECUTED (2026-07-13; work package 5 of GRAV.ENTROPIC.ACTION.01/v473) + THE EXPLICIT LORENTZIAN-POSITIVITY WITNESS -- decides v473's … · `v475_entropic_scalaron_gate.py`
- **GRAV.ENTROPIC.COMPRESS.01** [exakt, Teil bedingt oder offen] THE AP2 COMPRESSION CONJECTURE MADE WELL-POSED (2026-07-13; the last open work package of the v473 entropic-action-bridge round, exercised on an exactly solvable model) -- decides the … · `v476_entropic_compression_toy.py`
- **GRAV.ENTROPIC.SCALEFLOW.01** [exakt, Teil bedingt oder offen] THE SCALE-FLOW REPRESENTATION OF THE ENTROPIC ACTION (2026-07-13; the follow-up to GRAV.ENTROPIC.SCALARON.01/v475) -- the surviving interpretation 3 (KMS-spectral renormalisation) typed … · `v477_entropic_scaleflow.py`
- **GRAV.ENTROPIC.LEGS.01** [exakt, Teil bedingt oder offen] FIRST COMPUTABLE STEPS ON THE TWO REMAINING LEGS OF THE ENTROPIC-ACTION BRIDGE (2026-07-13; after v473-v477 the open items were the AP2 continuum statement and the derivation of the KMS … · `v478_entropic_continuum_legs.py`
- **GRAV.NONCIRCULAR.01** [Regel oder Synthese] THE BISOGNANO-WICHMANN DIRECTION FIREWALL (registered 2026-08-27, with audit census). BW is a theorem WITHIN an already-local Poincare-covariant QFT with suitable vacuum (wedge modular flow … · `-`
- **GRAV.SPIN2.EMERGENCE.01** [offen] THE SPIN-2 EMERGENCE CONTRACT [O] (registered 2026-08-28, review wave 3). DEMAND: from the SAME spectral determinant that carries the gauge fixpoints, a massless transversal spin-2 pole … · `v997_spin2_smg_witnesses.py`
- **GRAV.WEINBERG.WITTEN.01** [offen] THE WEINBERG-WITTEN FIREWALL [O] (registered 2026-08-30, review wave 6). BINDING BEFORE ANY spin-2 marker move. DEMAND: the microscopic Hamiltonian family must NOT possess a fundamental … · `-`
- **GRAV.PARENT.GAP.NO_GO.01** [exakt, Teil bedingt oder offen] GAPPED-PARENT SPIN-2 BOUNDARY (registered 2026-09-05, v1024; no parent marker move). A uniform strictly positive full-parent gap excludes a massless spin-2 pole. The exact … · `v1024_parent_internal_boundaries.py`
- **GRAV.SPIN2.LOCAL.CONSTRAINT.01** [exakt] Polynomial first-class constraints for a quadratic tensor target leave exactly two positive linear-dispersion modes at nonzero momentum. Staggered differences give a finite stencil and no … · `v1029_local_tensor_constraints.py`
- **GRAV.FREE.CURVATURE.FOCK.01** [exakt, Teil bedingt oder offen] Declared local tensor target after removal of the entire unstable global zero block: C=r^4 P_TT, positive pure quasi-free Weyl/Fock state with 2(V-1) real oscillators and thermodynamic then … · `v1031_quantum_tensor_curvature.py`
- **GRAV.FREE.CURVATURE.COVARIANT.01** [exakt, Teil bedingt oder offen] On the null cone the linearized Riemann sandwich of the de Donder form equals a rank-two positive TT polarization factorization. The polynomial kernel defines a positive-energy … · `v1032_covariant_curvature.py`
- **GRAV.LINEAR.SOURCE.WARD.01** [exakt, Teil bedingt oder offen] For the declared staggered tensor target and prescribed conserved sources use c0=-s^T q and ci=i Bp with the real paired Fourier Hamiltonian. The constraints propagate exactly; the explicit … · `v1035_matter_coupling.py`

### HOR (22)

- **HOR.01** [bedingt] horizon unit code: 1/(2pi)=4c3; 1920=|W(D5)|; S_dS; Page; beta_rad · `v8_horizon.py`
- **HOR.NARIAI.01** [exakt, Teil bedingt oder offen] THE MAXIMAL BLACK HOLE IS THE ANCHOR (carrier-in-the-bulk surface; Schwarzschild-de Sitter in seam units): (1) NARIAI = ANCHOR - at M_N = 1/(N_fam sqrt(L)) the horizon cubic is t^3-3t+2 = … · `v101_horizon_anchor.py`
- **HOR.ORIENT.01** [exakt, Teil bedingt oder offen] ONE ORIENTATION - the anchor is the stationary repeller in BOTH sectors (variational upgrade of the v101 temperature lemma + the v82/v99 flavor flow): FLAVOR [I] - the canonical flow dq/dt … · `v102_seam_orientation.py`
- **HOR.TRISECT.01** [exakt, Teil bedingt oder offen] TRISECTION NORMAL FORM (the canonical coordinate of the SdS moduli line; answers the HOR.ORIENT.01 curvature-matching question at the geometric level): the horizon cubic is uniformized by … · `v103_trisection_normal_form.py`
- **HOR.CLOCK.01** [exakt, Teil bedingt oder offen] THE CLASSICAL NARIAI CLOCK IS THE ANCHOR (third appearance) + honest (2/3)-test (executes the HOR.TRISECT.01 clock question at the classical level): STATIC PIN [I] - phi(rho) = rho solves … · `v104_nariai_clock.py`
- **HOR.CLOCK.02** [exakt, Teil bedingt oder offen] QUANTUM-CLOCK TARGET (R1 made quantitative; nothing fitted): the Ginsparg-Perry multi-mode clock has per-l equations lambda^2+lambda+(l(l+1)-2)=0 - l=0 anchor pair {1,-2}, l=1 zero modes … · `v107_quantum_clock_target.py`
- **HOR.CLOCK.03** [exakt, Teil bedingt oder offen] THE RESUMMED CLOCK (R1's bend acquires a closed form; the HOR.CLOCK.02 job description transforms): (1) THE SPECTRUM READS (1 - n/N)^{p2} [I]: the frozen transfer eigenvalues {1, (2/3)^6, … · `v124_resummed_clock.py`
- **HOR.CLOCK.04** [exakt, Teil bedingt oder offen] THE CLOCK AND THE WALL SHARE A SPECTRUM (the v124 resummed-clock weights ARE the Mehta-Seshadri parabolic weights; R1's first checkpoint passes classically): (1) THE WEIGHTS ARE THE … · `v126_clock_wall_bridge.py`
- **HOR.CLOCK.05** [exakt, Teil bedingt oder offen] THE CLOCK IS A LOG-DETERMINANT (the geometric tail of HOR.CLOCK.03/04 is the standard ring/RPA resummation; the 'O(1) resummation' demanded by HOR.CLOCK.02 is now TYPED): (1) LOG-DET … · `v127_ring_resummation.py`
- **HOR.CLOCK.06** [exakt, Teil bedingt oder offen] THE CLOCK IS AN ENTROPY POWER LAW (the transfer eigenvalues are the Nariai entropy fractions to the hexagon power; the ring coupling is the fractional entropy deficit): (1) THE THREE … · `v129_entropy_power_law.py`
- **HOR.CLOCK.07** [exakt, Teil bedingt oder offen] THE BORN SQUARE (the amplitude weight is h = N_fam = half the zero-mode count; the clock exponent p2 = 2h is its Born square - the HOR.CLOCK.05 multiplicity and the HOR.CLOCK.06 weight are … · `v130_born_square.py`
- **HOR.CLOCK.08** [exakt, Teil bedingt oder offen] THE MEASURE IS THE AREA (the per-mode S^{1/2} of the Born square is the zero-mode NORM, computed exactly; the R1 residue narrows to one determinant-ratio cancellation): (1) ZERO-MODE NORM = … · `v131_measure_is_area.py`
- **HOR.CLOCK.09** [exakt, Teil bedingt oder offen] THE DET-RATIO ANOMALY IS THE KOIDE CONSTANT (the HOR.CLOCK.08 'last step' acquires its first exact invariant): (1) SPECTRAL SETUP (re-pinned): the Ginsparg-Perry sphere operator on the unit … · `v132_detratio_anomaly.py`
- **HOR.CLOCK.10** [exakt, Teil bedingt oder offen] THE ZETA BUDGET, BOTH ROUTES COMPUTED EXACTLY (the reduced seam route carries pure anchor-ratio anomalies; the naive 4d route does not - the budget computation structurally selects the … · `v133_zeta_budget.py`
- **HOR.CLOCK.11** [exakt, Teil bedingt oder offen] THE VOLKOV-WIPF FIREWALL (R1 typed closure: the external full-4d graviton/ghost value is pinned and the edge part is exactly two copies of the reduced seam budget; external values verified … · `v138_vw_firewall.py`
- **HOR.CLOCK.12** [exakt, Teil bedingt oder offen] DET-RATIO CANCELLATION EXACT WITHIN THE SdS FAMILY (the v131 step closed in the within-family sense): e2-RIGIDITY - the traceless SdS cubic (t^3-3t+2c) forces r_b r_c = 1 - Delta^2/3 … · `v144_detratio_family.py`
- **HOR.CLOCK.13** [exakt, Teil bedingt oder offen] THE CLOCK IS A GAUSSIAN ZERO-MODE INTEGRAL AND THE BEND IS DET'-CLEAN (R1's derivation structured to its floor): (1) GAUSSIAN MODEL THEOREM - for p2 modes with per-mode variance ratio … · `v147_clock_gaussian.py`
- **HOR.NARIAI.02** [exakt] THE ONE-LINE NARIAI ENTROPY BOUND (external-review proposal 2026-06-14): the SdS entropy ratio S_tot/S_dS = (x^2+1)/Phi3(x) (x=r_b/r_c; Phi3=x^2+x+1 the N_fam=3 cyclotomic; HOR.NARIAI.01) … · `v190_nariai_entropy_bound.py`
- **HOR.BRANCHLINE.01** [bedingt] THE UNIVERSAL BRANCH LINE (external-review proposal 2026-06-14; HONESTLY typed [C], an exact affine RELABELING, NOT a theorem): the map q=7/2+(N_fam^2/2)m carries the SdS branch points … · `v191_universal_branch_line.py`
- **HOR.EHT.01** [exakt, Teil bedingt oder offen] EHT ACHROMATIC POLARIZATION-ROTATION INTERCEPT (2026-06-15, archive integration): beta_BH(r) = Q_e^eff Q_m^eff/(256 pi^4 r^2) = 16 c3^4 (Q_e Q_m / r^2) around a horizon-scale black hole. … · `v203_eht_achromatic.py`
- **HOR.BHTHERMO.01** [exakt, Teil bedingt oder offen] BLACK-HOLE THERMODYNAMICS FROM THE SEAM: modular Hawking temperature + scalaron-corrected Wald entropy (2026-06-15, archive integration of the dropped tfpt-42 … · `v208_bh_thermodynamics.py`
- **HOR.BHCORE.01** [offen] BLACK HOLE AS A SEAM FIXED POINT (topological defect), NOT a singularity (2026-06-15, archive integration of the old Five-Problems BH picture; the SUPERSEDED torsion/RN-vortex metric is … · `v209_bh_regular_core.py`

### QGAMB (13)

- **QGAMB.SADDLE.01** [exakt, Teil bedingt oder offen] Corrected 6 November 2026: Finite saddle diagnostics [E]/[C]/[O]. The scalar stiffness is 1/c3=8pi and the gap margin is 6 ln(3/2)-31/(4pi^2)>0. On the stated two positive OS modes, Tr' … · `v365_qg_oneloop_saddle.py`
- **QGAMB.REDUNDANCY.01** [offen] Track 2 (2026-06-23, next-plan v2): the AMBIENT REDUNDANCY statement -- the holographic route to Gravity-complete that SIDESTEPS building the general Euclidean QG measure (C7/QG.AMB.01). … · `v369_qgamb_redundancy.py`
- **QGAMB.ROADMAP.01** [Regel oder Synthese] OBLIGATION ROADMAP FOR QG.AMB.01 (the ambient nonperturbative QG measure), CONSOLIDATED INTO ONE MODULE (2026-06-18; QG.AMB.01 since reclassified as a DECOUPLED general QG problem, NOT a … · `v275_qgamb_roadmap.py`
- **QGAMB.TIERB.01** [exakt, Teil bedingt oder offen] SEAM-CALDERON -> (E8)_1 MATCHING CERTIFICATE (2026-06-18): the concrete Tier-B step of QG.AMB.01 (v275) -- computes the full (E8)_1 matching target and pins the residual to ONE bit (does … · `v277_seam_calderon_e8_match.py`
- **QGAMB.MODULAR.01** [exakt, Teil bedingt oder offen] MODULAR-DATA / ANYON-CONDENSATION LAYER OF THE Tier-B HOLOMORPHY BIT (2026-06-18): complements v277 (det+character) with the anyon/modular-data view. [E] #anyons = |det Gram|: E8 -> 1 … · `v281_holomorphy_modular_data.py`
- **QGAMB.UNIFY.01** [exakt, Teil bedingt oder offen] THE E8-AT-tau=i UNIFICATION: TWO OBLIGATIONS -> ONE OBJECT (2026-06-18): the candidate simplification. The (E8)_1 character as a modular function of the seam torus modulus tau is … · `v282_e8_tau_i_unification.py`
- **QGAMB.ROUTEII.01** [exakt, Teil bedingt oder offen] ROUTE (ii) (SEAM-NET CONDENSATION) LEMMA DECOMPOSITION + ROUTE-COINCIDENCE (2026-06-18): a 4-lemma chain with 3 discharged and the open lemma shown to COINCIDE with Route (i)'s (the v282 … · `v285_route_ii_seam_condensation.py`
- **QGAMB.IDG.01** [offen] INFINITE-DERIVATIVE / NONLOCAL READING OF THE GRAVITY SECTOR (2026-06-19): a CONDITIONAL narrowing of the red-team rt_F / v278 Stelle-ghost attack -- NOT a closure of QG.AMB.01 and NO new … · `v304_idg_nonlocal_ghost.py`
- **QGAMB.CLUSTER.01** [exakt, Teil bedingt oder offen] THE PHYSICAL QG MEASURE AS A GAP-DECOUPLED, CLUSTERING OBJECT (2026-06-21, TOE attack point 4): the honest statement in place of an over-claimed 4D measure for QG.AMB.01/G_metric. [E] the … · `v311_gap_decoupled_measure.py`
- **QGAMB.MEASURE.01** [exakt, Teil bedingt oder offen] CONSTRUCTING THE AMBIENT MEASURE ON THE GAP-DECOUPLED ADMISSIBLE SECTOR (2026-06-22) as a bona fide Osterwalder-Schrader object, with the honest fence (the constructive step v275's roadmap … · `v330_qgamb_admissible_measure.py`
- **QGAMB.METRIC.01** [exakt, Teil bedingt oder offen] THE METRIC SECTOR OF QG.AMB.01 CHARACTERIZED (2026-06-22): the obstruction is exactly the conformal-mode (Gibbons-Hawking-Perry) problem, the physical sector is gap-insulated, and the only … · `v332_qgamb_metric_sector.py`
- **QGAMB.CONFORMAL.01** [exakt, Teil bedingt oder offen] THE CONFORMAL-MODE RESOLUTION (2026-06-22): the Gibbons-Hawking-Perry contour rotation (IR sign) + the infinite-derivative (IDG) dressing (UV pole) give a CONDITIONALLY convergent … · `v334_conformal_resolution.py`
- **QGAMB.NODEP.01** [exakt, Teil bedingt oder offen] The No-Ambient-Dependency certificate (2026-06-25): no frozen TFPT readout is COMPUTED from the ambient QG measure QG.AMB.01 -- sharpening the v369/v384 prose redundancy ('QG.AMB.01 is a … · `v423_no_ambient_dependency.py`

### HYP (10)

- **HYP.0D.01** [exakt, Teil bedingt oder offen] HYPERGRAPH/WOLFRAM BRIDGE, HARDENED (2026-06-19): the 3-arm (2,3,r) branch type is the IRREDUCIBLE seed on which a LOCAL witness-based spectral growth dynamics REACHES E8 = T(2,3,5) as the … · `v299_hypergraph_0d_autonomous.py`
- **HYP.INJECT.01** [exakt, Teil bedingt oder offen] THE HYPERGRAPH SUBSTRATE SHARPENED -- WHAT A REWRITE MUST INJECT (2026-06-21, TOE attack point 5): turns the open rewrite-to-full-structure question into a precise list. [E] the affine-E8 … · `v312_hypergraph_substrate.py`
- **HYP.FIBER.01** [exakt, Teil bedingt oder offen] THE MINIMAL HYPERGRAPH SUBSTRATE THAT WORKS -- A (2,3,5)-NETWORK FIBERED BY A 3-FOLD CUSP CHANNEL (2026-06-22, the constructive v312 follow-up): reproduces BOTH the E8 skeleton AND the … · `v324_hypergraph_fiber.py`
- **HYP.REWRITE.01** [exakt, Teil bedingt oder offen] DERIVING THE CUSP WEIGHT 2/3 FROM A MINIMAL REWRITE (2026-06-22, sharpening HYP.FIBER.01/v324 + HYP.INJECT.01/v312). A minimal local rule -- a family channel with N_fam=3 slots, one … · `v327_hypergraph_rewrite.py`
- **HYP.COUPLED.01** [exakt, Teil bedingt oder offen] HYP.COUPLED.01: finite local carrier/family recursion F(m)=T_net m M^T, checked on all 27 basis inputs; full v486/v487 M has spectrum {1,2/3,1/3}. The joint Perron line is marks tensor e0. … · `v395_hypergraph_coupled.py`
- **HYP.PHI0.01** [exakt, Teil bedingt oder offen] THE phi0 LEADING TERM FROM ICOSAHEDRAL COMBINATORICS (2026-06-24, sharpening v312/v395): the tree-level retained seed phibase = 1/(6 pi) has a forced combinatorial reading on the … · `v396_phi0_icosahedral.py`
- **HYP.PHI0.PUNCTURE.01** [exakt, Teil bedingt oder offen] NAMES THE phi0 PUNCTURE TERM AS A TRACKED [O] HEAT-KERNEL TARGET (2026-06-24, v408): dtop = Omega_adm c3^4 = 48 c3^4 = 3/(256 pi^4) IS the order-4 (|mu4|) local boundary Seeley-DeWitt … · `v408_phi0_puncture_heatkernel.py`
- **HYP.PHI0.GEOM.01** [exakt, Teil bedingt oder offen] THE GEOMETRIC SIDE OF THE phi0 PUNCTURE TARGET IS EXACT (2026-07-14, deep-dive round; narrows HYP.PHI0.PUNCTURE.01/v408, does NOT close it): every twisted/orbifold heat trace of the flat … · `v483_pillowcase_exact_traces.py`
- **HYP.REWRITE.02** [exakt, Teil bedingt oder offen] THE FULL TRANSFER SPECTRUM FROM ONE FORCED LOCAL RULE (2026-07-15, generator-hunt round; sharpens HYP.REWRITE.01/v327): the verified seam transfer spectrum {1, (2/3)^6, (1/3)^6} (v54/v221) … · `v486_transfer_full_rule.py`
- **HYP.REWRITE.03** [exakt, Teil bedingt oder offen] THE LAZY SPLIT FORCED BY THE RESUMMED CLOCK (2026-07-15, closing the v486 [O] 'split selection'; sharpens HYP.REWRITE.02): the split (stay, hop) = (1/2, 1/6) is not a choice -- it is … · `v487_transfer_clock_rungs.py`

### FLATAWAY (7)

- **FLATAWAY.RP.01** [exakt, Teil bedingt oder offen] THE FLAT-AWAY LEMMA AS A NAMED MINI-THEOREM (2026-06-18). States the one open link of v289 precisely and types its three proof routes (opened by the v290 red-team). LEMMA: given RawRPSeam … · `v291_flataway_contract.py`
- **FLATAWAY.HEAT.01** [exakt, Teil bedingt oder offen] HEAT-KERNEL ROUTE FOR FLAT-AWAY (2026-06-18): reduces the open Flat-Away lemma to ONE named spectral invariant by showing the heat-trace deviation is a positive-definite quadratic form in … · `v292_flataway_heat_reduction.py`
- **FLATAWAY.SPEC.01** [exakt, Teil bedingt oder offen] SPECTRAL-RIGIDITY ROUTE FOR FLAT-AWAY (2026-06-18): upgrades the v291 one-direction scan to the FULL smooth-Z4 deformation space, showing the flat seam is a STRICT isolated minimum of the … · `v293_flataway_spectral_hessian.py`
- **FLATAWAY.ENERGY.01** [exakt, Teil bedingt oder offen] RP-ENERGY / TROYANOV ROUTE FOR FLAT-AWAY (2026-06-18): the third independent reduction -- the flat pillowcase is the UNIQUE constant-curvature (energy-minimising) metric in the conformal … · `v294_flataway_rp_energy.py`
- **FLATAWAY.A2.01** [exakt, Teil bedingt oder offen] EXACT ANALYTIC a_2 COEFFICIENT FOR FLAT-AWAY (2026-06-18): upgrades the v292 positive-definiteness of the heat-trace deviation from a numerical scan to a proof, and derives the exact … · `v295_flataway_a2_exact.py`
- **FLATAWAY.A2.CLOSED.01** [exakt] EXACT a_2 COEFFICIENT IN CLOSED FORM (2026-06-18): evaluates the v295 infinite mode sum analytically. W_k(t)=(1/2) sum_{n in Z} Phi(|n|,|n+4k|) with the divided-difference kernel of … · `v296_flataway_a2_closed.py`
- **FLATAWAY.RIGID.01** [exakt, Teil bedingt oder offen] THE FLAT-AWAY ATTACK -- HARDEN + DERIVE (2026-06-19): two moves on the single external fact the v291/v293/v295 routes left open ('pin the seam's Steklov spectrum / fix the a_2 heat … · `v300_flataway_rigid.py`

### QG (3)

- **QG.AMB.01** [offen] DECOUPLED GENERAL QG PROBLEM, NOT a TFPT-specific gate (reclassified 2026-06-22, v335): the ambient/metric-sector measure IS the general Euclidean-QG conformal-factor problem … · `-`
- **QG.G2.01** [exakt] G2 grounded: spectral-action Seeley-DeWitt for Dirac (E=-R/4, tr I=4) gives a2~-R/3 (Einstein-Hilbert) + a4~R^2/72 (R^2), i.e. R+R^2 STRUCTURE convention-independently; … · `v36_spectral_action_g2.py`
- **QG.GAPDEC.01** [bedingt] SIMPLER resolution: gap decoupling. 2||V_metric||=0.785<Delta=6log(3/2)=2.433 => Delta_eff=1.648>0; admissible sector keeps positive gap under metric dressing => decoupled from un-built … · `v36_spectral_action_g2.py`

### VGEO (2)

- **VGEO.SHARPEN.01** [exakt, Teil bedingt oder offen] v_geo IS THE SINGLE DIMENSIONFUL INPUT OF THE WHOLE THEORY (2026-06-23, Direction 6 of next-plan.md), now confirmed across matter AND gravity, and irreducibly primitive. The No-Unit Theorem … · `v364_vgeo_sharpen.py`
- **VGEO.TORSOR.01** [exakt, Teil bedingt oder offen] T1 CLOSURE -- the calibration/scale-torsor audit of the v_geo interface (one module, 21 checks, ~1 s; discovery probe phys_vgeo_torsor_audit_probe.py 21/21, verdict TORSOR-SCHLIESSUNG … · `v725_phys_vgeo_torsor_audit.py`

### BULK (1)

- **BULK.PRIME_COMPLETION.01** [offen] THE PRIME-COMPLETION CONTRACT [O] (registered 2026-09-04, v1020): replaces the "unique from seam data alone" reading of SEAM.BULK4D.RECON.01/B7 by uniqueness IN THE CATEGORY OF … · `v1020_spectator_sector_obstruction.py`

### CAUSAL (1)

- **CAUSAL.NOGO.01** [exakt, Teil bedingt oder offen] Causal Boundary Engineering - CONDITIONAL NO-GO for traversable seam shortcuts/CTCs (proof attempt that turned into a no-go). CORE [I] (machine-checked, v56-bound): a gapped positive … · `v64_causal_boundary_nogo.py`

### GRAVITY (1)

- **GRAVITY.COMPLETE.01** [exakt, Teil bedingt oder offen] GRAVITY-COMPLETE = BOUNDARY-COMPLETE + AMBIENT-REDUNDANCY (2026-06-24, Paper B / the gravity closure accounting): the explicit verdict combining the parameter-free LOCAL field equation … · `v399_gravity_complete.py`

### SCALE (1)

- **SCALE.OVERDET.01** [exakt] THE SINGLE MASS ANCHOR IS OVER-DETERMINED; NO DERIVABLE ABSOLUTE SCALE IS A THEOREM (2026-06-18): sharpens the No-Unit theorem (v153). [E] one anchor: every dimensionful output is a pure … · `v274_scale_overdetermination.py`

## Kosmologie, Vorhersagen und Daten

### FR (27)

- **FR.KOIDE.02** [bedingt] CONJECTURE: Q_l source=0.66446; Q_source+phi0/|W(A3)|=Q_source+phi0/24=0.666679 overshoots 2/3 by 1.3e-5 (source->pole correction conjecture, |W(A3)|=24) · `v25_frontier_conjectures.py`
- **FR.DM.02** [bedingt] determinant-line axion f_a=M_scal/(2 dim S+ |mu4|)=M_scal/128~2.39e11 GeV, m_a~23.8 ueV. POST-INFLATIONARY REFRAMING [C] (2026-07-10; isocurvature analysis … · `v25_frontier_conjectures.py`
- **FR.ETAB.01** [bedingt] baryon asymmetry eta_B = 6.1e-10 (downstream of Omega_b) · `-`
- **FR.MPME.01** [offen] proton/electron ratio m_p/m_e (cross-sector QCD) · `-`
- **FR.KOIDE.01** [bedingt] Koide Q = 0.664 (near 2/3; not exact) · `-`
- **FR.DM.01** [offen] dark matter: axion candidate fixed; relic scale (f_a; m_a) open · `-`
- **FR.KOIDE.03** [bedingt] sharper Koide source->pole conjecture: multiplicative seed transfer u->(1-1/(2 N_fam^3))u=(53/54)u (54=2*3^3=|Z2| N_fam^3, no new free number) gives Q=0.66666610, dev -5.7e-7 from 2/3 (10x … · `v37_plucker_anchor.py`
- **FR.POLE.SOLVE.01** [exakt, Teil bedingt oder offen] Track 3 (2026-06-23, F_transfer): the Koide source->pole transfer as a TYPED solver -- promotes experiments/ftransfer/koide_source_to_pole into the ledger with a kill test. Does NOT … · `v371_ftransfer_pole.py`
- **FR.BOLTZMANN.SOLVE.01** [bedingt] Track 3 (2026-06-23, F_transfer): eta_B from the INTEGRATED Buchmueller-Di Bari-Pluemacher Boltzmann network at the TFPT-frozen heavy scale -- promotes … · `v372_ftransfer_boltzmann.py`
- **FR.RELIC.SOLVE.01** [bedingt] Track 3 (2026-06-23, F_transfer): the axion relic abundance from the FINITE-T misalignment ODE -- promotes experiments/ftransfer/axion_relic and DECIDES the two competing TFPT angle … · `v373_ftransfer_relic.py`
- **FR.QCD.BUDGET.01** [exakt, Teil bedingt oder offen] Track 3 (2026-06-23, F_transfer): the m_p/m_e transfer as a TYPED uncertainty-budget + kill test -- the falsifiability companion to v262 (the precise 2-loop point), replacing the v326 toy … · `v374_ftransfer_qcd.py`
- **FR.KOIDE.04** [bedingt] Koide RG attractor is FORCED, not postulated [I] + ONE [P] identification (sharpens FR.KOIDE.01/02/03 via the v81 double cover): a branch-divisor-preserving Moebius map of the pencil line … · `v82_koide_attractor_splitting.py`
- **FR.KOIDE.05** [exakt, Teil bedingt oder offen] Koide relaxation toy (next.txt P2 NARROWED, not closed): (1) BASIN LEMMA [I] - under the canonical dictionary q=N_fam*Q the physical Koide range Q in [1/3,1] (Cauchy-Schwarz) maps to q in … · `v93_koide_relaxation_toy.py`
- **FR.KOIDE.06** [exakt, Teil bedingt oder offen] KOIDE FLOW TIME (P2 sharpened on both ends; corrects the ROBUSTNESS of v93 negative (i), not its arithmetic): [I] CANONICAL GENERATOR - the v82 flow is generated by dq/dt = … · `v99_koide_flow_time.py`
- **FR.KOIDE.07** [exakt, Teil bedingt oder offen] OPERATOR ORIGIN OF THE KOIDE SOURCE->POLE FACTOR 53/54 (external-review proposal 2026-06-14, validated; upgrades the FR.KOIDE.03 conjecture from a bare arithmetic guess 1-1/(2 N_fam^3) to … · `v183_koide_f_corner_transfer.py`
- **FR.ETAB.03** [bedingt] HONEST TEST OF THE ANCHORED-BOLTZMANN eta_B PROPOSAL (external review 2026-06-14; firewall verdict: SHARPER SCENARIO, not a closure): the proposal M_1 = M_R phi0^4 and m~1 = m_3/A_Lambda … · `v184_etaB_anchored_boltzmann.py`
- **FR.DM.03** [bedingt] AXION RELIC HONEST MISALIGNMENT ESTIMATE (external-review proposal 2026-06-14; the determinant-line axion can be all DM but stays a [C] scenario): (1) SCALES [I]: f_a = M_scal/(2 dim S+ … · `v185_axion_relic_solver.py`
- **FR.TRANSFER.01** [exakt] THE F_TRANSFER DISCIPLINE GUARD (a CI-style firewall check, the machine-readable form of the Frontier rule; external-review proposal 2026-06-14): the four continuous-transfer frontier … · `v187_ftransfer_laws.py`
- **FR.RAREKAON.01** [bedingt] RARE KAON SECTOR K -> pi nu nubar FROM THE TFPT CKM POINT (2026-06-15, archive integration): the kernel fixes the full CKM matrix (s12=lambda_C, s23=phi0/(1+lambda_C)=|V_cb|, … · `v202_rare_kaon.py`
- **FR.ETAB.02** [bedingt] ETA_B LEPTOGENESIS TRANSFER OPERATIONALISED (the Tier-3 F_transfer Boltzmann path made falsifiable; NOT a new structural class, NOT a compiler theorem): type-I thermal leptogenesis eta_B ~ … · `v169_etaB_boltzmann_interface.py`
- **FR.MUONG2.01** [bedingt] MUON g-2 SEAM VERTEX (2026-06-15, archive integration of the old muon-g2 note, re-typed for the current kernel): the carrier fixes a second-order topological defect delta_2 = … · `v204_muon_seam_g2.py`
- **FR.ETAB.04** [bedingt] LEPTOGENESIS SCALARON-DECUPLE BRANCH (2026-06-15): BOTH Boltzmann inputs share the decuple A_Lambda = 10 = |E(K5)| -- M_1 = M_scal phi0^2/A_Lambda ~ 8.65e9 GeV and m~1 = m_3/A_Lambda ~ 5 … · `v212_leptogenesis_decuple.py`
- **FR.MPME.02** [exakt, Teil bedingt oder offen] F_QCD TRANSFER REPRODUCTION OF m_p/m_e (2026-06-18): the contract-only F_transfer/qcd_matching_mp_me pipeline implemented as a real solver. 2-loop alpha_s(M_Z)=0.1179 run with n_f … · `v262_fqcd_mp_me.py`
- **FR.DYNAMICS.01** [Regel oder Synthese] THE DISCRETE->DYNAMIC LENS ON F_transfer (2026-06-19): asks whether the main-branch mechanism -- a gapped, positivity-preserving relaxation to a UNIQUE attractor (the local-averaging update … · `v303_ftransfer_dynamics.py`
- **FR.SUITE.01** [exakt, Teil bedingt oder offen] THE F_transfer SOLVER SUITE (2026-06-22): the four transfer readouts as ONE runnable harness, each a gapped relaxation to a unique attractor, with values + uncertainties + honest typing … · `v326_ftransfer_suite.py`
- **FR.CLOCK.PROBE.01** [exakt, Teil bedingt oder offen] A FALSIFIABLE probe of how external the 'external physics' really is (2026-06-24): do the frontier scheme-numbers land on the {2,3,5} clock (= the anchor atoms |Z2|,N_fam,g_car, v394) or … · `v397_external_clock_probe.py`
- **FR.POLE.QED.01** [exakt, Teil bedingt oder offen] F_pole, the external side (2026-06-25): the Koide source->pole transfer is a ONE-LOOP ELECTROMAGNETIC effect, not an algebraic step -- the honest 'do the external physics' follow-up to v93 … · `v420_fpole_qed_running.py`

### COSMO (6)

- **COSMO.LAM.01** [exakt] rho_Lambda/Mbar^4 = (3/4pi^2) e^{-2 a^-1} = 7.125e-121 (dark-energy scale 2.2375 meV) · `v7_gravity_cosmo.py`
- **COSMO.INF.01** [exakt] Starobinsky A_s; n_s=1-2/N*; r=12/N*^2 · `v7_gravity_cosmo.py`
- **COSMO.OMB.01** [exakt] Omega_b = (1-1/4pi) phi0 = 0.04894 · `v7_gravity_cosmo.py`
- **COSMO.NSTAR.01** [bedingt] N_star band -> POINT [P] (the announced 'future derivation from a reheating model' executed; demotion of N_star as a compiler output STANDS): the only compiler input is … · `v86_nstar_reheating.py`
- **COSMO.RUNNING.01** [exakt, Teil bedingt oder offen] Spectral-index running of the Starobinsky branch, PROMOTED from the 2026-07-20 scorecard prediction-of-record row (same N_star_reheating group as n_s/r/A_s; never an independent hit). [E] … · `v494_cosmo_running_mu_record.py`
- **COSMO.MUDIST.01** [bedingt] CMB mu-distortion as the A_s/N_star BRANCH DISCRIMINATOR, PROMOTED from the 2026-07-20 scorecard prediction-of-record row (data-limited; shares the Nstar_branch alternative group with the … · `v494_cosmo_running_mu_record.py`

### PRED (3)

- **PRED.LAYER.01** [exakt] Prediction/falsification layer with EXPLICIT failure thresholds (Alessandro residual gate #6) [N], data-confrontation Python-only. Each prediction carries an explicit falsification … · `v65_falsification_layer.py`
- **PRED.LIVE.01** [bedingt] LIVE KILL-TEST SCORECARD OF THE RECENT ROUND'S COMPUTABLE PREDICTIONS (2026-06-18): a data-confrontation consolidation (confronts, does not add physics), each row tagged with the decisive … · `v283_live_killtest_scorecard.py`
- **PRED.JOINTLIKELIHOOD.01** [offen] THE JOINT LIKELIHOOD CONTRACT [O] (registered 2026-08-27; from hit-counting to overdetermined relations). The 27+ individually matched values descend from the same few atoms (phi_0, c_3, … · `v992_joint_likelihood_v1.py`

### CCC (2)

- **CCC.SEAM.KINEMATICS.01** [exakt] THE CCC CROSSOVER KINEMATICS DICTIONARY [E] (rounds DL-DLV, notes 2026-08-24; promoted 2026-08-24). EXACT (v957; six probes 28/28 + 29/29 + 18/18 + 21/21 + 14/14 + 7/7 embedded byte-exact, … · `v957_ccc_crossover_kinematics.py`
- **CCC.SEAM.CROSSOVER.01** [offen] THE CCC CROSSOVER CONTRACT [O] (registered 2026-08-24 from the external TFPT-vs-CCC adjudication of the six gates; the kinematic base is CCC.SEAM.KINEMATICS.01 [E]/v957). DEMAND (the … · `v957_ccc_crossover_kinematics.py`

### DATA (2)

- **DATA.SCORECARD.01** [exakt] Confrontation of TFPT predictions with current research data (2024/25) - honest scorecard incl. tensions [N], data-confrontation (Python-only, external data). STRONG MATCHES: … · `v62_data_scorecard.py`
- **DATA.WATCHDOG.01** [exakt] DATA WATCHDOG / DECISION PIPELINE (2026-06-21): the operational decision layer over the frozen registry (REG.FREEZE.01/v84) [N] (Python-only, external data). [N] reads … · `v307_data_watchdog.py`

### DM (2)

- **DM.AXION.SPINE.01** [bedingt] AXION DM SPINE-ANGLE BRANCH theta_i = pi N_fam/g_car = 3 pi/5 = 108 deg (2026-06-15): an ALTERNATIVE, more ROBUST misalignment angle than the determinant-line hilltop theta_i ~ 170 deg … · `v211_axion_spine_angle.py`
- **DM.AXION.PENTAGON.01** [exakt, Teil bedingt oder offen] THE 'UNMAPPED' ICOSAHEDRAL/GOLDEN E8 STRUCTURE SURFACES IN THE ONE EXTERNAL COSMOLOGICAL INPUT theta_i (2026-06-26, v429): a geometric re-reading answering the parsimony question 'why does … · `v429_axion_pentagon_phi.py`

### AS (1)

- **AS.RECONCILE.01** [exakt, Teil bedingt oder offen] A_s is predicted (not missing), is DIMENSIONLESS (so it does NOT define an absolute mass), and the apparent -11 sigma tension is a reheating-channel question -- reconciled with the archived … · `v340_As_reconcile.py`

### KILL (1)

- **KILL.BOARD.01** [exakt] THE FORWARD KILL-TEST BOARD: THE DECISIVE UPCOMING MEASUREMENTS, RANKED (2026-06-22, forward companion to v307's current-data board). [N] assembles the decisive NEXT measurements, each … · `v321_killtest_board.py`

### LAMBDA (1)

- **LAMBDA.BRANCH.01** [exakt, Teil bedingt oder offen] Lambda-metrology branch selection => G_N as horizon-metrology OUTPUT, the 123-orders double overlap, and the ACT DR6 birefringence comparison. (1) The form … · `v60_lambda_metrology_branch.py`

### METROLOGY (1)

- **METROLOGY.CLOSURE.01** [exakt, Teil bedingt oder offen] THE FINAL INPUT SET {a, pi, v_geo} CERTIFIED (2026-06-24, Paper D / the scale-modulus closure): 0 dimensionless dials + 1 dimensionful unit, consolidating v153 (No-Unit), v274 … · `v401_metrology_closure.py`

### OBS (1)

- **OBS.TRANSDUCTION.01** [bedingt] THE TRANSDUCTION PRINCIPLE CONTRACT (registered 2026-08-27; elevates the S15 typing rule to a named principle). PRINCIPLE: observable signature = internal eigenvalue [compiler] + … · `v990_wj_transduction_mechanism.py`

### OBSERVATORY (1)

- **OBSERVATORY.REGISTRY.01** [exakt] Track 4 (2026-06-23): the prediction OBSERVATORY -- a status-typed CI over the frozen prediction registry (verification/freeze_file.csv) that makes the falsifiability surface … · `v375_observatory_registry.py`

### UVBRANCH (1)

- **UVBRANCH.KILLTEST.01** [exakt, Teil bedingt oder offen] The optional carrier-Pati-Salam UV branch on the ALL-ORDER footing -- a single citable kill-test surface with the proton-decay SAFETY hierarchy computed (2026-06-24). A consolidation: does … · `v385_uvbranch_killtest.py`

## Methodik, Pruefung und Sonstiges

### REDTEAM (6)

- **REDTEAM.A.01** [bedingt] Red Team / Stress Test Target A - adversarial attack on the (E8)_1 boundary-net identification (GATE.METRIC.03). CONFIRMS the conservative wording: central charge is NECESSARY not … · `redteam/rt_A_e8net.py`
- **REDTEAM.B.01** [exakt] Red Team Target B - attack on the UPSTREAM physical premise of Theorem A (the Lean arithmetic 2^g=g^2+g+2 <=> g=5 is NOT attacked; it is verified [F]). Finding: the carrier rank is selected … · `redteam/rt_B_pascal.py`
- **REDTEAM.C.01** [exakt] Red Team Target C - dimensional firewall on k=c3/2. The reduced coefficient k_red=c3/2=1/(16pi) is dimensionless; the only legal chain is k_phys=k_red/G, S=4pi k_phys A_phys=A_phys/(4G). … · `redteam/rt_C_kc3.py`
- **REDTEAM.D.01** [exakt] Red Team Target D - stress test of the U_point->v_geo bijection (FLAV.UPOINT.01). For positive reals (ratios,product) <=> (individuals) is a bijection and the TFPT sectors satisfy it … · `redteam/rt_D_upoint.py`
- **REDTEAM.E.01** [exakt] Red Team Target E - single-scale audit of v_geo (ANCHOR.VGEO.01). The closed compiler + protected IR readouts all reduce to (pure number) x v_geo: M_scal=c3^(7/2) Mbar, f_a=(c3/|mu4|)^(7/2) … · `redteam/rt_E_vgeo.py`
- **REDTEAM.F.01** [bedingt] RED TEAM Target F -- adversarial audit of the perturbative 4D-QFT + scale round (v269-v275), published in release 5.2. TWO attacks land and force a narrowing (now folded into v269/v274); … · `redteam/rt_F_qft4d.py`

### CONTRACT (4)

- **CONTRACT.G.01** [offen] Research Contract 2 (G_metric): lemma chain G0-G8 + Theorem G; priority 2 (deep analytic programme). Certified: G5 gap-dominance 2||V||=0.785<Delta=2.433, Delta_eff>1.647 · `v29_research_contract_certs.py`
- **CONTRACT.F.01** [bedingt] THE F_transfer FUNCTOR CONTRACT (2026-06-15): the four frontier transfers are ONE typed functor F_transfer = F_observable o F_threshold o F_RG satisfying four structural axioms, each … · `v213_ftransfer_functor.py`
- **CONTRACT.QFT4D.01** [exakt, Teil bedingt oder offen] 4D LAGRANGIAN VIA THE SPECTRAL ACTION (2026-06-17): the honest CONTRACT for the Phase-3 structural gap (the 2d chiral seam net v156/v175 + its dynamics v238-v243 -> the physical 4d theory). … · `v244_spectral_action_contract.py`
- **CONTRACT.QFT4D.DIRAC.01** [exakt, Teil bedingt oder offen] THE SPECTRAL-TRIPLE CONTRACT (2026-06-17): the object that would turn 'carrier is gauged' from a fork into a theorem (Reviewer-A Theorem 3 + the sigma=scalaron identification). Target: a … · `v250_dirac_triple_contract.py`

### CORRECTIONS (3)

- **CORRECTIONS.GAP.01** [Regel oder Synthese] The practical harvest of v383 (2026-06-24): the same spectral gap that makes each sector parameter-free ALSO sets the SIZE of its first correction -- for a gapped operator the deviation … · `v387_corrections_gap.py`
- **CORRECTIONS.BUDGET.01** [exakt, Teil bedingt oder offen] The gap correction (v387) is a TYPED budget, NOT a uniform band on every prediction (2026-06-24): the (lambda_2/lambda_1)^n size is a property of the distance from a gapped operator's … · `v388_corrections_budget.py`
- **CORRECTIONS.NUMERIC.01** [exakt, Teil bedingt oder offen] The typed correction budget (v388) turned into ACTUAL first-correction magnitudes -- a computed theoretical-uncertainty surface (2026-06-24), the short-term testable harvest of v388: … · `v393_corrections_numeric.py`

### AUDIT (2)

- **AUDIT.WORDING.01** [exakt] THE FRONTIER-WORDING GUARD (a prose sentinel, no new physics; external-review proposal 2026-06-14): v187 protects the ledger TYPING, this guard protects the PROSE. It forbids stale … · `v188_frontier_wording_guard.py`
- **AUDIT.NULLMC.01** [exakt] NUMEROLOGY NULL TEST (look-elsewhere corrected; executes the legacy 'numerology grammar test' TODO from the CMB paper): (A) EXACT CENSUS - for each of the 13 scored frozen-registry … · `v100_numerology_null_mc.py`

### REVIEW (2)

- **REVIEW.REPRO.02** [Regel oder Synthese] EXTERNAL REPRODUCTION of the v83 package (Alessandro, 2026-06): TeX build.sh notes = 9 ok 0 failed; Python run_all v1-v83 = ALL CHECKS PASSED; Wolfram tfpt_readouts.wl = 116/116 ALL WOLFRAM … · `verification/make_manifest.py (--check)`
- **REVIEW.MAP.01** [exakt, Teil bedingt oder offen] THE FOUR EXTERNAL-REVIEW CONCERNS REDUCE TO THE TWO ALREADY-NAMED RESIDUALS + ONE DATA-ONLY ITEM (a careful 2026-06-14 review raised A mapping 2D->4D, B UV-gravity abstractness, C … · `v182_reviewer_residual_map.py`

### RH (2)

- **RH.SPLIT.01** [exakt] RH route progress: exponent collapse sum_k U^k A0 U^-k = 4 diag(A0) exact (U=diag(1,i,-i)); integer-exponent options = perms of (2,1,1),(2,2,0); SPLITTING type made ALGEBRAIC: O(-2)+O(-1)^2 … · `v32_rh_splitting.py`
- **RH.BUNDLE.01** [exakt] EXPLICIT valid (U_wall) flat bundle constructed via numerical RH solve: A0 (eig {0,1/3,2/3}, diag (1/2,1/4,1/4)) realises cusp class + splitting O(-2)+O(-1)^2 (SigmaA_k=diag(2,1,1)) + … · `v33_explicit_flat_bundle.py`

### TOPO (2)

- **TOPO.E8.01** [exakt, Teil bedingt oder offen] THE SEAM AS THE E8 KLEINIAN SINGULARITY (2026-06-17): a canonical algebraic-geometric MODEL for the open seam-realisation premise QGEO.REALIZE.01, the du Val side of the McKay … · `v232_e8_kleinian_seam.py`
- **TOPO.BRIESKORN.01** [exakt] BRIESKORN CAPSTONE (2026-06-17): the (2,3,5) Brieskorn singularity x^2+y^3+z^5 is the ONE generator of the discrete skeleton -- the capstone of the icosahedral round … · `v236_brieskorn_capstone.py`

### EXTREV (1)

- **EXTREV.LATTICE.01** [exakt, Teil bedingt oder offen] The external-review lattice audit (one module, 14 checks, <1 s; discovery probe external_lattice_audit_probe.py 14/14, verdict LATTICE-AUDIT-PASSED; external review 2026-08-02, three … · `v624_external_lattice_audit.py`

### MESSAGE (1)

- **MESSAGE.LADDER.01** [exakt] v844 -- MESSAGE.LADDER.01 + DOILY.PASCAL.RANK.01: the E8 addressing formula and the Doily-Pascal rank theorem -- two exact corpus compressions, ONE module from two probes (23/23 + 16/16 … · `v844_message_doily_rank.py`

### REG (1)

- **REG.FREEZE.01** [exakt] Blind-prediction registry FROZEN 2026-06-09, machine-enforced (review point: pre-registration before JUNO/CMB-S4/LiteBIRD/DESI). predictions_frozen.json freezes every dimensionless … · `v84_frozen_registry.py`

### TRANSLATE (1)

- **TRANSLATE.CLOCK.01** [Regel oder Synthese] THE TRANSLATION CLOCK: THE DISCRETE<->DYNAMIC BRIDGE IS THE ORDER-30 COXETER ELEMENT, A CLOCK WITH TWO COPRIME HANDS 5 x 6 (2026-06-22, synthesis of the v313-v318 arc + v124/v55, answering … · `v319_translation_clock.py`

### WITNESS (1)

- **WITNESS.ECON.01** [exakt] WITNESS-INDEPENDENCE / GENERATOR-ECONOMY AUDIT (2026-06-21): the STRUCTURAL anti-numerology firewall complementing the v100 phenomenological null test. [E] GENERATOR ECONOMY: 13 headline … · `v305_witness_independence.py`
