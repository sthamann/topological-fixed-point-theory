import json,sys
from pathlib import Path
repo=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4');sys.path.insert(0,str(repo/'semantic_explorer'));import extract
corp=extract.Corpus(repo);path=repo/'semantic_explorer/data/static-candidate-resolutions.json';doc=json.loads(path.read_text())
notes={
('PRIME.LORENTZ.SPINOR.01','verification/v807_lorentz_nullselector.py',302,'CR_CONST_CV'):{
 'semantic_role':'author_chosen_diagnostic_threshold','purpose':'Coefficient-of-variation gate used only by the LORENTZ-COORDS-REVEAL classification when the monotonicity gate fails. It is a diagnostic threshold chosen in this script, not a TFPT-derived constant or independently calibrated physical/statistical constant.', 'selection_basis':'author_chosen_test_threshold','selection_rationale':'The v807 source labels the threshold but does not justify its numerical value; it is also absent from the explicit frozen-bar list. Treat it as an author-selected diagnostic bar.', 'limits':'The reported enum depends on this threshold; no physical conclusion follows from its passing.'},
('PRIME.LORENTZ.SPINOR.01','verification/v807_lorentz_nullselector.py',1036,'a_anchor'):{
 'semantic_role':'conditional_integer_split_anchor','purpose':'The chosen anchor tuple (1,1,2) has elementary symmetric values (4,5,2), matching (|mu4|, g_car, N_fam) in the stated positive-integer rank-3 split. v807 checks consequences of that anchor; it does not derive the physical/module identification.', 'selection_basis':'conditional_on_integer_split_and_typing','selection_rationale':'origin_theory.tex states the positive 3-multiset and its elementary-symmetric typing; v499 derives bundle consequences while retaining module identification/(U) as residual. Thus the tuple is forced only within that typing, not unconditionally selected by physical origin.', 'limits':'Physical identification of the integer split remains open; v807 is a consequences check, not an anchor-origin proof.'},
('PRIME.GAUGE.PPR.01','verification/v556_gauge_ppr_identities.py',1020,'m_c'):{
 'semantic_role':'finite_matrix_test_size','purpose':'Sets one matrix/test size for the parity-sine eigenbasis and Toeplitz identity computation.', 'selection_basis':'finite_numerical_test_instance','selection_rationale':'The companion checks v557 use m=64,128,256; 128 is one finite test size, not the separate TFPT 128-spinor extension and not a physical parameter.', 'limits':'The selected finite size only identifies the tested matrix instance.'},
('SEAM.BIT.FREEDOM.01','verification/v510_seam_bit_freedom.py',664,'(L, W)'):{
 'semantic_role':'finite_cylinder_instance','purpose':'Defines the discrete cylinder Z_(2L) x {0,...,W-1} used in this finite deck-involution test.', 'selection_basis':'one_finite_model_instance','selection_rationale':'The source comment identifies the cylinder; odd W supplies the unique central seam row needed by the construction, but this particular L=12,W=5 instance is not derived from TFPT and W=5 is not tied to g_car in this code.', 'limits':'The code supports the construction for L>0 and odd W; this run only checks one instance.'},
('PRIME.RANK3ZERO.01','verification/v684_rank3_zeroside.py',153,'S_T_BAR'):{
 'semantic_role':'declared_scan_tolerance','purpose':'Tolerance for the scan-count and Riemann-von-Mangoldt checkpoint deviations.', 'selection_basis':'declared_diagnostic_budget','selection_rationale':'The script labels this a declared budget and uses it to permit at most three-count deviations; it is not evidence that no zero was missed.', 'limits':'A finite checkpoint tolerance does not certify scan completeness.'},
('PRIME.RANK3ZERO.01','verification/v684_rank3_zeroside.py',156,'DELTA_POS'):{
 'semantic_role':'declared_zero_position_error_budget','purpose':'Declared position-error budget used in the explicit-formula zero terms.', 'selection_basis':'declared_numerical_budget_without_per_root_certificate','selection_rationale':'The source labels it a declared budget. The implemented bisection width certifies roots of the Riemann-Siegel approximation, but the script does not provide per-root true-zeta interval/error checks that establish this 1e-3 bound.', 'limits':'This positional-error premise is not established by the reported approximate-root checks; downstream zero-side bounds remain conditional on it.'},
('PRIME.RANK3ZERO.01','verification/v684_rank3_zeroside.py',157,'GABCKE_C'):{
 'semantic_role':'unused_cited_constant_declaration','purpose':'Declares the cited Gabcke C0 coefficient 0.127 for the classical Riemann-Siegel remainder.', 'selection_basis':'classical_source_cited_but_unused_here','selection_rationale':'The script comment states the bound, and the TeX cites the Gabcke remainder; inspection of v684 shows this symbol is not read after assignment. Its presence does not support the numerical zero scan in this module.', 'limits':'Unused declaration; no result in this script is grounded by this constant.'},
('PRIME.RANK3ZERO.01','verification/v684_rank3_zeroside.py',161,'RH_VERIFIED_T'):{
 'semantic_role':'external_verification_height','purpose':'Upper height through which the cited Platt–Trudgian verification supports RH on the critical line; the script uses it to classify heights in its declared surface.', 'selection_basis':'external_cited_verification_input','selection_rationale':'The script explicitly cites Platt–Trudgian (2021); the paper text also identifies beta=1/2 up to 3e12. This supports the on-line classification through that height, not completeness or location of the script’s numerically scanned zero list.', 'limits':'External verified-height input only; does not certify this script’s zero scan.'},
('PRIME.RANK3ZERO.01','verification/v684_rank3_zeroside.py',162,'TAIL_SLACK'):{
 'semantic_role':'declared_tail_margin','purpose':'Multiplies the leading Riemann-von-Mangoldt density integral used in the analytic tail estimate.', 'selection_basis':'declared_numerical_margin','selection_rationale':'The module calls this a declared numeric estimate; no explicit global upper envelope for the zero-count remainder is applied here to prove that factor 1.10 is a rigorous bound.', 'limits':'Tail estimate is numerical/conditional, not a certified rigorous remainder bound.'},
('PRIME.RANK3ZERO.01','verification/v684_rank3_zeroside.py',163,'N_QUAD_CELL'):{
 'semantic_role':'quadrature_order','purpose':'Sets four Gauss–Legendre points per lag cell in the RTERM calculation.', 'selection_basis':'numerical_discretization_setting','selection_rationale':'The source comment identifies the order. No interval quadrature error is charged in the script, so the computed term is numerical rather than exact.', 'limits':'Numerical quadrature error remains unbounded by this setting alone.'}
}
extra={
 ('PRIME.LORENTZ.SPINOR.01',302): [('verification/v807_lorentz_nullselector.py',59,90),('verification/v807_lorentz_nullselector.py',658,662),('verification/v807_lorentz_nullselector.py',801,805)],
 ('PRIME.LORENTZ.SPINOR.01',1036): [('origin_theory.tex',121,129),('verification/v499_p2_weights_deligne_bg.py',368,382),('verification/v807_lorentz_nullselector.py',1036,1061)],
 ('PRIME.GAUGE.PPR.01',1020): [('verification/v556_gauge_ppr_identities.py',1016,1035),('verification/v557_cascade_vector_identities.py',842,862)],
 ('SEAM.BIT.FREEDOM.01',664): [('verification/v510_seam_bit_freedom.py',662,689)],
 ('PRIME.RANK3ZERO.01',153): [('verification/v684_rank3_zeroside.py',151,164),('tfpt_3_e8_audit_bootstrap.tex',4325,4338)],
 ('PRIME.RANK3ZERO.01',156): [('verification/v684_rank3_zeroside.py',151,164),('verification/v684_rank3_zeroside.py',227,241),('verification/v684_rank3_zeroside.py',448,463)],
 ('PRIME.RANK3ZERO.01',157): [('verification/v684_rank3_zeroside.py',151,164),('tfpt_3_e8_audit_bootstrap.tex',4325,4338)],
 ('PRIME.RANK3ZERO.01',161): [('verification/v684_rank3_zeroside.py',151,164),('tfpt_3_e8_audit_bootstrap.tex',4325,4338)],
 ('PRIME.RANK3ZERO.01',162): [('verification/v684_rank3_zeroside.py',151,164),('verification/v684_rank3_zeroside.py',403,412)],
 ('PRIME.RANK3ZERO.01',163): [('verification/v684_rank3_zeroside.py',151,164),('verification/v684_rank3_zeroside.py',289,297)]
}
for case in doc['cases']:
 key=(case['claim_id'],case['source']['path'],case['source']['line'],case['candidate_name'])
 if key in notes:
  case.update(notes[key])
  case['method']='source_reviewed_mapping'
  case['review_scope']='exact assignment, enclosing computation, downstream use, and specified original source windows; no script execution'
  case['supporting_sources']=[]
  for p,a,b in extra[(case['claim_id'],case['source']['line'])]:
   sp=repo/p;source=corp.source(sp,a,b)
   case['supporting_sources'].append({'id':source['id'],'path':source['path'],'sha256':source['sha256'],'line_start':a,'line_end':b,'excerpt':source['text']})
# Add explicit rule for every remaining case: role and selection metadata must be present.
assert len(doc['cases'])==195 and all(c.get('semantic_role') and c.get('purpose') and c.get('selection_rationale') for c in doc['cases'])
path.write_text(json.dumps(doc,ensure_ascii=False,indent=2)+'\n')
print('annotated',len(notes),'cases')
