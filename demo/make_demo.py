#!/usr/bin/env python3
"""
Generate a fully synthetic GPBLUP demonstration dataset (no real-animal data).

Creates, deterministically (fixed seed):
  ped.txt   : pedigree            "animal sire dam"          (0 = unknown parent)
  geno.txt  : genotypes           "animal dosage-string"     (one 0/1/2 char per SNP)
  phen.txt  : phenotypes (header) "animal dummy hy sex ADG BF age litter"

A small pig-like nucleus is simulated: 60 founders + 4 generations of 160, with
genotypes gene-dropped through the pedigree (so relatives share genotypes and the
genomic relationship matrix reflects the pedigree). Two correlated traits (ADG, BF)
are simulated from marker effects plus fixed (herd-year, sex), covariate (age),
litter (IID), and residual terms, matching model_demo.par.
"""
import numpy as np

rng = np.random.default_rng(2026)

# ---- dimensions ----
N_FOUND, N_GEN, PER_GEN = 60, 4, 160
NSNP = 800
G0 = np.array([[12.0, 3.0], [3.0, 9.0]])     # genetic (co)variance
R0 = np.array([[30.0, 5.0], [5.0, 25.0]])    # residual (co)variance
VAR_LITTER = 5.0                              # litter (IID) variance, trait 1
MU = np.array([120.0, 14.0])                  # trait means (ADG days-ish, BF mm-ish)
B_AGE = np.array([-0.20, 0.05])               # age covariate slope per trait

# ---- pedigree (gene-dropped) ----
ids, sires, dams, sexes, gens = [], [], [], [], []
hap1 = {}; hap2 = {}
p = rng.uniform(0.15, 0.85, NSNP)             # founder allele frequencies

def draw_founder():
    return (rng.random(NSNP) < p).astype(np.int8), (rng.random(NSNP) < p).astype(np.int8)

nid = 0
gen_members = []
# generation 0: founders
g0 = []
for _ in range(N_FOUND):
    nid += 1
    ids.append(nid); sires.append(0); dams.append(0)
    sexes.append(int(rng.integers(1, 3))); gens.append(0)
    hap1[nid], hap2[nid] = draw_founder()
    g0.append(nid)
gen_members.append(g0)

# generations 1..N_GEN
for g in range(1, N_GEN + 1):
    prev = gen_members[g - 1]
    males = [i for i in prev if sexes[i - 1] == 1] or prev
    females = [i for i in prev if sexes[i - 1] == 2] or prev
    cur = []
    for _ in range(PER_GEN):
        nid += 1
        s = int(rng.choice(males)); d = int(rng.choice(females))
        ids.append(nid); sires.append(s); dams.append(d)
        sexes.append(int(rng.integers(1, 3))); gens.append(g)
        # Mendelian sampling: one allele from each parent haplotype per locus
        pick_s = rng.integers(0, 2, NSNP)
        pick_d = rng.integers(0, 2, NSNP)
        h_s = np.where(pick_s == 0, hap1[s], hap2[s])
        h_d = np.where(pick_d == 0, hap1[d], hap2[d])
        hap1[nid], hap2[nid] = h_s.astype(np.int8), h_d.astype(np.int8)
        cur.append(nid)
    gen_members.append(cur)

N = nid
gens = np.array(gens); sexes = np.array(sexes)
Z = np.vstack([hap1[i] + hap2[i] for i in ids]).astype(np.int8)   # [N, NSNP] dosage 0/1/2

# ---- true breeding values from marker effects ----
sum2pq = float(np.sum(2.0 * p * (1.0 - p)))
Umark = rng.multivariate_normal(np.zeros(2), G0 / sum2pq, size=NSNP)   # [NSNP, 2]
Zc = Z - 2.0 * p                                                       # centered
A = Zc @ Umark                                                          # [N, 2] breeding values

# ---- fixed / covariate / litter / residual ----
n_hy = 6
hy = rng.integers(1, n_hy + 1, N)
hy_eff = rng.normal(0, np.sqrt(5.0), size=(n_hy + 1, 2))
sex_eff = np.array([[0, 0], [0.0, 0.0], [2.0, -1.0]])                  # index by sex (1,2)
age = np.clip(rng.normal(170, 15, N), 120, 220).round().astype(int)
n_litter = 130
litter = rng.integers(1, n_litter + 1, N)
litter_eff = rng.normal(0, np.sqrt(VAR_LITTER), size=n_litter + 1)
E = rng.multivariate_normal(np.zeros(2), R0, size=N)

Y = (MU + hy_eff[hy] + sex_eff[sexes] + np.outer(age - 170, B_AGE)
     + A + np.column_stack([litter_eff[litter], np.zeros(N)]) + E)

# ---- write files ----
# genotyped = generations >= 2 (later cohorts); phenotyped = generations >= 1
geno_mask = gens >= 2
phen_mask = gens >= 1

with open("ped.txt", "w") as f:
    for i in range(N):
        f.write(f"{ids[i]} {sires[i]} {dams[i]}\n")

with open("geno.txt", "w") as f:
    for i in range(N):
        if geno_mask[i]:
            f.write(f"{ids[i]} " + "".join(map(str, Z[i].tolist())) + "\n")

with open("phen.txt", "w") as f:
    f.write("animal dummy herd_year sex ADG BF age litter\n")
    for i in range(N):
        if phen_mask[i]:
            f.write(f"{ids[i]} 0 {hy[i]} {sexes[i]} {Y[i,0]:.2f} {Y[i,1]:.2f} "
                    f"{age[i]} {litter[i]}\n")

with open("demo_truth.txt", "w") as f:                 # hidden true breeding values
    f.write("animal geno true_BV_ADG true_BV_BF\n")
    for i in range(N):
        f.write(f"{ids[i]} {int(geno_mask[i])} {A[i,0]:.4f} {A[i,1]:.4f}\n")

print(f"animals={N}  genotyped={int(geno_mask.sum())}  phenotyped={int(phen_mask.sum())}  "
      f"SNP={NSNP}")
print(f"var(BV) ADG={A[:,0].var():.1f} BF={A[:,1].var():.1f}  "
      f"(target {G0[0,0]:.0f}/{G0[1,1]:.0f})")
print("wrote ped.txt, geno.txt, phen.txt")
