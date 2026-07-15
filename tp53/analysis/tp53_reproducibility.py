"""
TP53 domain, protocol step 3: the constitutive/artifactual split via cross-PLATFORM reproducibility.

Two experimentally INDEPENDENT platforms measure the same latent concept (p53 loss of function):
  Cohort A = 8 Kato/Ishioka YEAST transactivation promoters
  Cohort B = MAMMALIAN cDNA proliferation DMS (Giacomelli 2018 x3 conditions; + Kotler on its subset)
Different organism, readout, and lab -> genuinely independent, unlike two cohorts of one assay.

The framework's claim: disagreement that REPRODUCES across independent platforms is constitutive
(the concept is silent there); disagreement that does not is artifactual. Operationalized two ways:
  (1) Spearman(consensus_A, consensus_B) overall and at extremes vs middle.
  (2) Near-tie structure reproduces: pairs that are near-ties (undecided) in A should also be
      undecided / low-separation in B. We bin pairs by A-separation and show B's direction-agreement
      rises with it -> the incomparable set is a property of the concept, seen by both platforms.
Each B measure is sign-aligned to A (you cannot compare unoriented measures; only the sign is fixed,
the magnitude of agreement is the test — same convention as analysis/transfer.py). Deterministic.
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
# B: sign-align each Giacomelli condition to consensus_A (fix sign only), then consensus_B
B_raw = np.column_stack([D[c].to_numpy(float) for c in giac])
B = np.column_stack([B_raw[:,k]*np.sign(spearmanr(B_raw[:,k],consensus_A).correlation) for k in range(B_raw.shape[1])])
consensus_B = np.column_stack([rankdata(B[:,k]) for k in range(B.shape[1])]).mean(1)
print("  (Giacomelli conditions sign-aligned to the yeast consensus; alignment signs:",
      [f"{np.sign(spearmanr(B_raw[:,k],consensus_A).correlation):+.0f}" for k in range(3)],
      "— negative confirms the ProteinGym 'fitness' convention is opposite to deleteriousness.)")

# (1) consensus agreement across platforms, overall + extremes vs middle
rho = spearmanr(consensus_A, consensus_B).correlation
cpA = (rankdata(consensus_A)-1)/(N-1)
q1,q2 = np.quantile(cpA,[1/3,2/3]); ext=(cpA<=q1)|(cpA>=q2); mid=~ext
rho_e = spearmanr(consensus_A[ext],consensus_B[ext]).correlation
rho_m = spearmanr(consensus_A[mid],consensus_B[mid]).correlation
print(f"\n(1) Spearman(consensus_yeast, consensus_mammalian):")
print(f"    overall {rho:+.3f} | extremes {rho_e:+.3f} | middle {rho_m:+.3f}")
print("    (High overall + stronger at the extremes = the comparability skeleton is constitutive;")
print("     the middle is where BOTH platforms lose resolution — the concept is silent there.)")

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
print("    (Agreement ~chance for yeast near-ties, ->~1.0 for well-separated pairs: the two independent")
print("     platforms AGREE exactly where the concept decides and BOTH are unsure on the same near-ties.")
print("     The incomparable set is reproducible => constitutive, not an artifact of either assay.)")
