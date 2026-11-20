"""LEPTONTRANSFER20261119: the map's sharpest untried empirical point.  Source ratios
m_mu/m_tau = (8/7) phi0, m_e/m_mu = (12/7) phi0^2 against pole masses (+2.20 %, +0.21 %).
Does one-loop QED running to a common scale close the gap?  Masses quoted from memory."""
from common import *
from mpmath import sqrt
c = Checks("lepton_transfer")
phi0 = PHI_TREE + DTOP
me, mmu = mpf("0.51099895"), mpf("105.6583755")
al = 1 / CODATA
out = {}
for mtau, dm, tag in ((mpf("1776.86"), mpf("0.12"), "PDG 2022"), (mpf("1776.93"), mpf("0.09"), "newer average")):
    a_src, b_src = mpf(8) / 7 * phi0, mpf(12) / 7 * phi0**2
    a_pol, b_pol = mmu / mtau, me / mmu
    ga, gb = a_src / a_pol - 1, b_src / b_pol - 1
    qa, qb = -3 * al / (2 * pi) * log(mtau / mmu), -3 * al / (2 * pi) * log(mmu / me)
    koide = lambda a, b: (1 + a + a * b) / (1 + sqrt(a) + sqrt(a * b))**2
    print(f"  [{tag}] gap mu/tau {mp.nstr(ga * 100, 5)} %  e/mu {mp.nstr(gb * 100, 5)} %   one-loop QED common scale: {mp.nstr(qa * 100, 4)} %, {mp.nstr(qb * 100, 4)} %")
    print(f"           Koide Q: source {mp.nstr(koide(a_src, b_src), 8)}  pole {mp.nstr(koide(a_pol, b_pol), 8)}  (2/3 = 0.66666667)")
    t = (ga - 3 * al) / (dm / mtau)
    print(f"           gap mu/tau minus 3 alpha = {mp.nstr(ga - 3 * al, 4)}  = {mp.nstr(t, 3)} sigma of m_tau")
    out[tag] = {"gap_mu_tau": mp.nstr(ga, 8), "gap_e_mu": mp.nstr(gb, 8), "qed_mu_tau": mp.nstr(qa, 6), "qed_e_mu": mp.nstr(qb, 6),
                "koide_source": mp.nstr(koide(a_src, b_src), 10), "three_alpha_pull": mp.nstr(t, 4)}
c.ok("gaps reproduce the map: +2.20 % and +0.21 %", abs(ga - mpf("0.0220")) < 3e-4 and abs(gb - mpf("0.0021")) < 2e-4)
c.ok("one-loop QED running to a common scale has the WRONG sign (-0.98 %, -1.86 %): it widens both gaps", qa < 0 and qb < 0 and abs(qa + mpf("0.0098")) < 2e-4 and abs(qb + mpf("0.0186")) < 2e-4)
c.ok("source ratios do not satisfy Koide (Q differs from 2/3 by more than 1e-3), pole masses do", abs(koide(a_src, b_src) - mpf(2) / 3) > 1e-3 and abs(koide(a_pol, b_pol) - mpf(2) / 3) < 1e-4)
c.ok("numerical remark only: the mu/tau gap equals 3 alpha within 2 sigma of m_tau; the e/mu gap (0.29 alpha) shows no such form", abs(mpf(out["newer average"]["three_alpha_pull"])) < 2 and abs(gb / al - mpf("0.29")) < 0.02)
from mpmath import findroot
a_k = findroot(lambda a: koide(a, b_src) - mpf(2) / 3, a_pol)
pred = mmu / a_k
print(f"  Koide + source e/mu alone: m_mu/m_tau = {mp.nstr(a_k, 9)} ({mp.nstr((a_k / a_pol - 1) * 100, 3)} %),  m_tau = {mp.nstr(pred, 7)} MeV ({mp.nstr((pred - mtau) / dm, 3)} sigma)")
c.ok("if Koide is taken as the interface, the source e/mu ratio alone gives m_mu/m_tau to 0.02 % (instead of 2.2 %): two gaps reduce to one", abs(a_k / a_pol - 1) < 3e-4)
c.ok("but not to measurement precision: predicted m_tau = 1777.30 MeV is about 4 sigma high, driven by the 0.21 % e/mu gap", 3 < (pred - mtau) / dm < 5, mp.nstr((pred - mtau) / dm, 4))
out["koide_route"] = {"mmu_over_mtau": mp.nstr(a_k, 10), "mtau_pred_MeV": mp.nstr(pred, 8)}
c.done({"data": out, "verdict": "QED running to a common scale EXCLUDED as the transfer (wrong sign). '3 alpha' is an unexplained coincidence with look-elsewhere, not a result"})
