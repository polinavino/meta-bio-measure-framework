"""
TP53 domain, protocol step 3: the constitutive/artifactual split via cross-PLATFORM reproducibility.

Two experimentally INDEPENDENT platforms measure the same latent concept (p53 loss of function):
  Cohort A = 8 Kato/Ishioka YEAST transactivation promoters
  Cohort B = MAMMALIAN cDNA proliferation DMS (Giacomelli 2018, three conditions, on 2,308 core variants)
Different organism, readout, and lab -> genuinely independent, unlike two cohorts of one assay.

We report cross-platform consensus agreement and separation-dependent direction agreement.
These summaries neither compare the full incomparable sets nor establish a concept-level cause.
Directions are fixed by the ProteinGym fitness convention, not fitted to yeast agreement.
"""
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

# --- persist this script's console output to a tracked file (single source of truth for paper numbers) ---
import sys as _sys
_OUT = "/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/analysis/outputs/reproducibility.txt"
class _Tee:
    def __init__(self, p): self._f = open(p, "w"); self._o = _sys.stdout
    def write(self, s): self._o.write(s); self._f.write(s)
    def flush(self): self._o.flush(); self._f.flush()
_sys.stdout = _Tee(_OUT)

M = pd.read_csv("/Users/polina/Documents/BioInfStuff/meta-bio-measure-framework/tp53/data/tp53_matrix.csv")
prom = ["waf1","mdm2","bax","h1433s","aip1","gadd45","noxa","p53r2"]
giac = ["giac_wtnut","giac_nullnut","giac_nulletop"]
keep = M[prom+giac].notna().all(axis=1).to_numpy()
D = M[keep].reset_index(drop=True); N=len(D)
print(f"=== Cross-platform reproducibility (yeast promoters vs mammalian DMS), N={N} ===\n")

# concept-positive (higher = more deleterious). A: negate transactivation.
A = np.column_stack([-D[p].to_numpy(float) for p in prom])
consensus_A = np.column_stack([rankdata(A[:,k]) for k in range(A.shape[1])]).mean(1)
# ProteinGym fitness direction is opposite to deleteriousness in these three conditions.
B_raw = np.column_stack([D[c].to_numpy(float) for c in giac])
B = -B_raw
consensus_B = np.column_stack([rankdata(B[:,k]) for k in range(B.shape[1])]).mean(1)
print("  Giacomelli signs fixed in advance: [-1, -1, -1] (fitness -> deleteriousness)")

# (1) consensus agreement across platforms, overall + extremes vs middle
rho = spearmanr(consensus_A, consensus_B).correlation
cpA = (rankdata(consensus_A)-1)/(N-1)
q1,q2 = np.quantile(cpA,[1/3,2/3]); ext=(cpA<=q1)|(cpA>=q2); mid=~ext
rho_e = spearmanr(consensus_A[ext],consensus_B[ext]).correlation
rho_m = spearmanr(consensus_A[mid],consensus_B[mid]).correlation
print(f"\n(1) Spearman(consensus_yeast, consensus_mammalian):")
print(f"    overall {rho:+.3f} | extremes {rho_e:+.3f} | middle {rho_m:+.3f}")

# Pooling the bottom and top tertiles spans the whole axis while the middle spans one
# third, and Spearman is attenuated by range restriction, so the pooled contrast above
# is partly an artifact of the split. Each tertile on its own carries the same
# restriction and is the fair comparison.
bot = cpA <= q1
top = cpA >= q2
print("    per-tertile (same range restriction in each, so comparable):")
for _nm, _m in [("bottom (functional)", bot), ("middle", mid), ("top (deleterious)", top)]:
    print(f"      {_nm:22s} {spearmanr(consensus_A[_m], consensus_B[_m]).correlation:+.3f}  (n={_m.sum()})")
_rng = np.random.default_rng(0)
_ib, _im, _it = np.flatnonzero(bot), np.flatnonzero(mid), np.flatnonzero(top)
_dp, _df = [], []
for _ in range(2000):
    _ra = _rng.choice(np.flatnonzero(ext), ext.sum(), True)
    _rm = _rng.choice(_im, _im.size, True)
    _rb = _rng.choice(_ib, _ib.size, True)
    _rt = _rng.choice(_it, _it.size, True)
    _mm = spearmanr(consensus_A[_rm], consensus_B[_rm]).correlation
    _dp.append(spearmanr(consensus_A[_ra], consensus_B[_ra]).correlation - _mm)
    _df.append(0.5*(spearmanr(consensus_A[_rb], consensus_B[_rb]).correlation
                    + spearmanr(consensus_A[_rt], consensus_B[_rt]).correlation) - _mm)
for _nm, _d in [("pooled extremes - middle", _dp), ("mean(functional,deleterious) - middle", _df)]:
    _lo, _hi = np.percentile(_d, [2.5, 97.5])
    print(f"    {_nm}: {np.mean(_d):+.3f}  bootstrap 95% CI [{_lo:+.3f}, {_hi:+.3f}]")
print("    (The agreement at the extremes is carried by the deleterious end; the functional")
print("     end agrees no better than the middle, so this is a floor effect shared by both")
print("     platforms, not a symmetric extremes-versus-middle pattern.)")

# (2) near-tie structure reproduces across platforms
i,j = np.triu_indices(N,1)
sepA = np.abs(cpA[i]-cpA[j])
dirA = np.sign(consensus_A[i]-consensus_A[j])
dirB = np.sign(consensus_B[i]-consensus_B[j])
agree = (dirA==dirB)
print("\n(2) Direction agreement between platforms, by yeast-platform separation:")
qs=np.quantile(sepA,[0,.2,.4,.6,.8,1.0])
for a,b in zip(qs[:-1],qs[1:]):
    m=(sepA>=a)&(sepA<b if b!=qs[-1] else sepA<=b)
    print(f"    |Δconsensus_yeast| in [{a:.3f},{b:.3f}]: cross-platform agreement={agree[m].mean():.3f}  (n={m.sum():,})")
print("    These are direction-agreement summaries of two scalar consensuses, not a direct")
print("    replication test of their incomparable sets. Near-chance agreement cannot distinguish")
print("    biological ambiguity from assay resolution and measurement error.")
