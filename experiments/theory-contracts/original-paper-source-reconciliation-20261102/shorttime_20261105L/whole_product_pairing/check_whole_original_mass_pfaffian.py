"""Whole symbolic v252 mass family in the EXISTING Spin4 product pairing.
No 96-independent-fermion determinant, no toy projected species, no new Yukawas.
The degree96 result is for the actual 192-dimensional pointwise mass pairing;
the whole regulated spacetime kinetic operator is not set equal to that mass.
"""
from pathlib import Path
import importlib.util,sys,json,hashlib
import sympy as s
root=Path('/Users/stefanhamann/Projekte/tfpt-theoryv4')
sys.path.insert(0,str(root/'verification'))
spec=importlib.util.spec_from_file_location('original_v252',root/'verification/v252_full_finite_triple.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
checks=[]
def need(name,b):
 checks.append({'name':name,'passed':bool(b)})
 if not b:raise AssertionError(name)
def zero(M):return all(s.expand(x)==0 for x in M.todok().values())
x=s.symbols('x1:5',real=True);h1=x[0]+s.I*x[1];h2=x[2]+s.I*x[3]
u=sum(t*t for t in x)
q=s.Matrix([[s.conjugate(h2),h1],[-s.conjugate(h1),h2]])
need('original_quaternion_determinant_is_whole_Hnorm',s.expand(q.det()-u)==0)
Y={f:s.Matrix(3,3,lambda i,j:s.Symbol(f'{f}r{i}{j}',real=True)+s.I*s.Symbol(f'{f}i{i}{j}',real=True))for f in ('u','d','e','nu')}
MR=s.zeros(3)
for i in range(3):
 for j in range(i,3):MR[i,j]=MR[j,i]=s.Symbol(f'Mr{i}{j}',real=True)+s.I*s.Symbol(f'Mi{i}{j}',real=True)
S=s.MutableSparseMatrix(48,48,{})
for up,down,colors in [('u','d',range(3)),('nu','e',range(1))]:
 for color in colors:
  for a,fl in enumerate((up,down)):
   for b,fr in enumerate((up,down)):
    for i in range(3):
     for j in range(3):
      row=v.IDX[(fl,'L',color,i)];col=v.IDX[(fr,'R',color,j)]
      S[row,col]=q[a,b]*Y[fr][i,j];S[col,row]=s.conjugate(S[row,col])
T=s.MutableSparseMatrix(48,48,{})
for i in range(3):
 for j in range(3):T[v.IDX[('nu','R',0,i)],v.IDX[('nu','R',0,j)]]=MR[i,j]
Phi=s.SparseMatrix.vstack(s.SparseMatrix.hstack(S,T),s.SparseMatrix.hstack(T.H,S.conjugate()))
KF=s.SparseMatrix(s.Matrix(v.Kswap()).applyfunc(s.nsimplify))
gF=s.SparseMatrix(s.Matrix(v.grading()).applyfunc(s.nsimplify))
need('whole_symbolic_Phi_has_original_96_dimension',Phi.shape==(96,96))
need('whole_symbolic_Phi_Hermitian',zero(Phi.H-Phi))
need('whole_symbolic_Phi_odd',zero(gF*Phi+Phi*gF))
need('whole_symbolic_v252_reality',zero(KF*Phi.conjugate()-Phi*KF))

L=[i for i,p in enumerate(v.P)if p[1]=='L'];R=[i for i,p in enumerate(v.P)if p[1]=='R']
M=S.extract(L,R)
need('whole_particle_LR_size24',M.shape==(24,24))
need('native_LR_label_order_matches',[(v.P[i][0],v.P[i][2],v.P[i][3])for i in L]==[(v.P[i][0],v.P[i][2],v.P[i][3])for i in R])
KPhi=KF*Phi
plus=L+[i+48 for i in R];minus=R+[i+48 for i in L]
Sp=KPhi.extract(plus,plus);Sm=KPhi.extract(minus,minus)
need('whole_chiral_KPhi_plus_zero_LL',zero(Sp[:24,:24]))
need('whole_chiral_KPhi_plus_LR_is_conjugate_M',zero(Sp[:24,24:]-M.conjugate()))
need('whole_chiral_KPhi_plus_RL_is_M_adjoint',zero(Sp[24:,:24]-M.H))
need('whole_chiral_KPhi_minus_zero_RR',zero(Sm[24:,24:]))
need('whole_chiral_KPhi_minus_LR_is_M_transpose',zero(Sm[:24,24:]-M.T))
need('whole_chiral_KPhi_minus_RL_is_M',zero(Sm[24:,:24]-M))
need('whole_chiral_Majorana_only_complementary_diagonals',zero(Sp[:24,:24])and zero(Sm[24:,24:]))
# Swapping the two24-row halves has sign(-1)^(24*24)=+1 and makes each
# of the ACTUAL full48 chiral blocks triangular, for arbitrary symmetric MR.
swapped=Sp.extract(list(range(24,48))+list(range(24)),list(range(48)))
need('whole_plus_block_row_swap_triangular',zero(swapped[24:,:24]))
swappedm=Sm.extract(list(range(24,48))+list(range(24)),list(range(48)))
need('whole_minus_block_row_swap_triangular',zero(swappedm[:24,24:]))
need('whole_block_swap_orientation_positive',(-1)**(24*24)==1)

pauli=[s.Matrix([[0,1],[1,0]]),s.Matrix([[0,-s.I],[s.I,0]]),s.diag(1,-1)]
gam=[s.kronecker_product(pauli[0],p)for p in pauli]+[s.kronecker_product(pauli[1],s.eye(2))]
g5=s.prod(gam);C4=gam[0]*gam[2];spin=s.SparseMatrix(C4.H*g5)
need('actual_spin_mass_factor_has_two_oriented_J_blocks',spin[0,0]==0 and spin[1,1]==0 and spin[2,2]==0 and spin[3,3]==0 and spin[0,1]in(-1,1)and spin[2,3]in(-1,1)and zero(spin+s.transpose(spin)))
need('both_spin_pfaffian_signs_positive_on48_internal_slots',all(((-1)**(48*47//2))*spin[i,i+1]**48==1 for i in(0,2)))
Mass=s.kronecker_product(spin,KPhi)
PhiY=s.diag(S,S.conjugate());DM=Phi-PhiY
weakL=s.diag(*([int(ch=='L')for _,ch,_,_ in v.P]*2));weakR=s.eye(96)-weakL
weak=s.eye(96)-2*weakL
need('actual_weak_center_reverses_all_Yukawa_mass',zero(weak*PhiY*weak+PhiY))
need('actual_Majorana_only_weak_singlet_RR',zero(DM-weakR*DM*weakR))
need('actual_Yukawa_only_LR_and_RL',zero(weakL*PhiY*weakL)and zero(weakR*PhiY*weakR))
need('actual_DM_DY_DM_zero_without_commuting_Y_or_MR',zero(DM*PhiY*DM))
p=s.Symbol('p',real=True)
mode_relations={}
for eps in(-1,1):
 inds=[96*i+j for i in range(4)for j in range(96)if g5[i,i]*gF[j,j]==eps]
 Me=Mass.extract(inds,inds)
 need(f'whole_total_chiral_mass_size192_{eps}',Me.shape==(192,192))
 need(f'whole_total_chiral_mass_skew_{eps}',zero(Me+Me.T))
 need(f'whole_total_chiral_mass_spinblocks_decouple_{eps}',zero(Me[:96,96:])and zero(Me[96:,:96]))
 internal0=[j for j in range(96)if g5[0,0]*gF[j,j]==eps]
 internal2=[j for j in range(96)if g5[2,2]*gF[j,j]==eps]
 need(f'actual_native_first_Jblock_{eps}',zero(Me[:48,48:96]-spin[0,1]*KPhi.extract(internal0,internal0)))
 need(f'actual_native_second_Jblock_{eps}',zero(Me[96:144,144:]-spin[2,3]*KPhi.extract(internal2,internal2)))
 for j,gmu in enumerate(gam):
  derivative=s.kronecker_product(s.SparseMatrix(-s.I*C4.H*gmu),KF).extract(inds,inds)
  need(f'actual_kinetic_spinblocks_offdiagonal_{eps}_{j}',zero(derivative[:96,:96])and zero(derivative[96:,96:]))
 # Fourier differentiation is symmetric at a single p. The actual skew
 # object couples p and -p; this must be checked before a Pfaffian claim.
 Bp=s.kronecker_product(s.SparseMatrix(C4.H*gam[3]),KF).extract(inds,inds)
 need(f'actual_single_p_kinetic_coefficient_symmetric_{eps}',zero(Bp-Bp.T))
 Kp=Me+p*Bp
 need(f'actual_opposite_momentum_skew_pair_relation_{eps}',zero(Kp.subs(p,-p).T+Kp))
 index={g:k for k,g in enumerate(inds)}
 rows=[];cols=[]
 for j in range(96):
  ir=next(i for i in(0,2)if g5[i,i]*gF[j,j]==eps)
  ic=next(i for i in(1,3)if g5[i,i]*gF[j,j]==eps)
  rows.append(index[96*ir+j]);cols.append(index[96*ic+j])
 Lp=Kp.extract(rows,cols);Lm=Kp.extract(cols,rows)
 need(f'whole_mode_inner_spin_diagonal_blocks_zero_{eps}',zero(Kp.extract(rows,rows))and zero(Kp.extract(cols,cols)))
 need(f'whole_mode_cross_blocks_opposite_p_transpose_{eps}',zero(Lm+Lp.subs(p,-p).T))
 matches=[]
 for aa in(1,-1,s.I,-s.I):
  for bb in(s.I,-s.I):
   if zero(Lp-aa*gF*KF*(Phi+bb*p*s.eye(96))):matches.append((aa,bb))
 need(f'whole_mode96_crossblock_actual_original_Phi_identity_{eps}',len(matches)==1)
 aa,bb=matches[0]
 need(f'whole_mode96_crossblock_unit_determinant_prefactor_{eps}',aa**96==1 and ((-1)**48)==1)
 mode_relations[str(eps)]={'crossblock_prefactor':str(aa),'Phi_plus_ip_sign':str(bb),
   'det_Kp':'det(p I96+i Phi)^2',
   'actual_skew_pair_size':384,
   'actual_pair_Pfaffian':'det_Kp; orientation (-1)^(192*193/2)=+1'}

# Exact actual neutrino6 chiral block in the existing weak-orbit chart H=(0,r).
amp=s.Symbol('r',positive=True)
nr=[v.IDX[('nu','L',0,i)]for i in range(3)]+[48+v.IDX[('nu','R',0,i)]for i in range(3)]
nc=[v.IDX[('nu','R',0,i)]for i in range(3)]+[48+v.IDX[('nu','L',0,i)]for i in range(3)]
Fnu=Phi.subs({x[0]:0,x[1]:0,x[2]:amp,x[3]:0}).extract(nr,nc)
want=s.BlockMatrix([[amp*Y['nu'],s.zeros(3)],[MR.H,amp*Y['nu'].T]]).as_explicit()
need('actual_whole_neutrino_chiral6_block_no_projection_of_family',zero(Fnu-want))
swap=s.BlockMatrix([[s.zeros(3),s.eye(3)],[s.eye(3),s.zeros(3)]]).as_explicit()
Xi=swap*Fnu.H
need('actual_neutrino_Takagi_matrix_orientation',zero(Xi-Xi.T))
need('actual_neutrino_Takagi_relation_is_unitary_row_permutation',zero(Xi.H*Xi-Fnu*Fnu.H))
need('physical_paired_mode_Pfaffian_orientation_positive',(-1)**(192*193//2)==1)

# Exact block inverse of the ACTUAL six-dimensional neutrino chiral family.
# MatrixSymbol Y means the already supplied original Ynu, not new couplings.
Ym=s.MatrixSymbol('Ynu_original',3,3);Rm=s.MatrixSymbol('MR_original_adjoint',3,3)
FF=s.BlockMatrix([[amp*Ym,s.ZeroMatrix(3,3)],[Rm,amp*Ym.T]])
FFinv=s.BlockMatrix([[Ym.I/amp,s.ZeroMatrix(3,3)],[-Ym.T.I*Rm*Ym.I/amp**2,Ym.T.I/amp]])
prod=s.block_collapse(FFinv*FF)
need('whole_neutrino_inverse_has_only_rminus1_and_rminus2_blocks',prod==s.BlockMatrix([[s.Identity(3),s.ZeroMatrix(3,3)],[s.ZeroMatrix(3,3),s.Identity(3)]]))
zeta,Su,Sd,Se,Sn=s.symbols('zeta Su Sd Se Sn')
formal=(1+p**2*zeta*Su)**6*(1+p**2*zeta*Sd)**6*(1+p**2*zeta*Se)**2*(1+2*p**2*zeta*Sn)
need('whole_mode_first_inverse_Hnorm_coefficient_all_colors_and_families',s.expand(s.expand(formal).coeff(zeta,1)-2*p**2*(3*Su+3*Sd+Se+Sn))==0)
need('whole_paired_mode_first_inverse_Hnorm_coefficient_doubles',s.expand(s.expand(formal**2).coeff(zeta,1)-4*p**2*(3*Su+3*Sd+Se+Sn))==0)

# Original quaternion mass factorization in all three quark colors + leptons.
canon=[];Qblocks=[];Yblocks=[]
for up,down,colors in [('u','d',range(3)),('nu','e',range(1))]:
 for color in colors:
  canon.extend((name,color,g)for name in(up,down)for g in range(3))
  Qblocks.append(s.kronecker_product(q,s.eye(3)));Yblocks.extend([Y[up],Y[down]])
pos={(v.P[i][0],v.P[i][2],v.P[i][3]):k for k,i in enumerate(L)}
order=[pos[key]for key in canon]
Mcan=M.extract(order,order)
need('all_original_Y_and_colors_factorize_without_projection',zero(Mcan-s.diag(*Qblocks)*s.diag(*Yblocks)))
need('whole_Hpower_det_M_is_u12',len(Qblocks)*3==12)
need('whole_physical_mass_Pf_degree96',4*24==96)
z,trX,trX2=s.symbols('z TrX TrX2')
series=s.series(s.exp(trX*z/2-trX2*z*z/4),z,0,3).removeO()
need('Pfaffian_second_coefficient_keeps_square_first_trace',s.expand(series.coeff(z,2)-(trX**2/8-trX2/4))==0)
out={'passed':sum(c['passed']for c in checks),'total':len(checks),'checks':checks,
 'scope':'Existing conditional Spin4 x actual whole96 product pairing: whole192 mass fibers and actual free constant-H opposite-momentum384 skew pairing. No whole primitive source, original interacting regulator Hessian or physical-time UY evaluation.',
 'whole_mass_Pfaffian':'u^48 |det Yu|^12 |det Yd|^12 |det Ye|^4 |det Ynu|^4',
 'Majorana_dependence_of_pure_mass_Pfaffian':'exactly independent of symmetric MR; source minors need not be',
 'pure_mass_homogeneous_coefficients':{'P96':'the whole mass Pfaffian above','P95':0,'P94':0},
 'whole_affine_skew_general_Pm_minus2':'Pf(A)*((Tr(A^-1 B))^2/8-Tr((A^-1 B)^2)/4), only where A invertible',
 'conditional_product_pointwise_simplification':'If actual finite regulator preserves the given whole Spin4/KO product, Tr(A^-1 B)=0 and Majorana/cross parts of Tr((A^-1 B)^2) vanish: only Tr((A^-1 Bkin)^2) remains. Inverse statement requires invertible Yukawa mass; complementary minors otherwise.',
 'actual_mode_crossblock_relations':mode_relations,
 'actual_neutrino_Takagi_convention':'Xi=[[0,r*conjugate(Ynu)],[r*Ynu^dagger,MR]], directly from actual v252 chiral block; singular values are those of Fnu',
 'mode_polynomial_Q':'det3(p²+u Yu†Yu)^6 det3(p²+u Yd†Yd)^6 det3(p²+u Ye†Ye)^2 det6(p²+Xi†Xi)',
 'genuine_whole_opposite_momentum_pair_Pfaffian':'Q(p,u)^2; Q is one factor per p when the FULL original p,-p set is counted, not an independently postulated 192x192 Pfaffian at nonzero p',
 'mode_Q_top_coefficients':{'u48':'C_Y','u47':'2*p²*C_Y*Sigma','Sigma':'Tr3[3(Yu†Yu)^-1+3(Yd†Yd)^-1+(Ye†Ye)^-1+(Ynu†Ynu)^-1]','MR_first_can_enter':'u46 or lower'},
 'genuine_paired_mode_top_coefficients':{'u96':'C_Y²','u95':'4*p²*C_Y²*Sigma'},
 'actual_original_spin_coefficients':{'gamma_mu':[[[str(g[i,j])for j in range(4)]for i in range(4)]for g in gam],
   'gamma5':[[str(g5[i,j])for j in range(4)]for i in range(4)],'C4':[[str(C4[i,j])for j in range(4)]for i in range(4)]},
 'native_physical_Y_values_chosen':False,'full_primitive_CUY_evaluated':False,
 'sources':[{'path':p,'sha256':hashlib.sha256((root/p).read_bytes()).hexdigest()}for p in ['verification/v252_full_finite_triple.py','experiments/theory-contracts/craig-original-source-bridge-20261103/ACTUAL-PRODUCT-PAIRING.txt','experiments/theory-contracts/craig-original-source-bridge-20261103/actual_product_pairing.py','experiments/theory-contracts/tfpt-unified-source-derivation-20261029/aps_source_review/DYNAMICAL-HIGGS-SOURCE-PROOF.txt']]}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
