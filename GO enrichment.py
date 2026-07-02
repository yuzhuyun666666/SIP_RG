"""
GO Biological Process (BP) enrichment with Fisher's exact test + BH FDR
Inputs (exploded mapping):
  - /mnt/data/C1_expanded.csv
  - /mnt/data/universe_expanded.csv
Each file must contain:
  - A gene identifier column (e.g., Gene/SYMBOL/id/ENSEMBL...) -> auto-detected
  - A "BP" column with GO Biological Process terms (one row = one (Gene, BP) pair)
Output:
  - /mnt/data/go_enrichment_results.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact

# -------- paths (修改为你的路径也可以，比如 Windows 本地) --------
c1_path = Path("C:/Users/YZY/Desktop/C12_expanded.csv")
universe_path = Path("C:/Users/YZY/Desktop/universe_expanded.csv")
out_path = Path("C:/Users/YZY/Desktop/go_enrichment_C12_results.csv")

# -------- helpers --------
def guess_gene_col(df):
    """
    Try to find a reasonable gene column.
    Falls back to the first non-'BP' column.
    """
    candidates = [
        "gene", "Gene", "GENE",
        "symbol", "Symbol", "SYMBOL",
        "gene_id", "GeneID",
        "ENSEMBL", "Ensembl",
        "id", "ID"
    ]
    for c in candidates:
        if c in df.columns:
            return c
    for c in df.columns:
        if c != "BP":
            return c
    return df.columns[0]

def prep_df(path):
    """
    Load CSV, keep only (Gene, BP), strip spaces, drop empties.
    """
    df = pd.read_csv(path)
    if "BP" not in df.columns:
        raise ValueError(f"{path.name} 缺少 'BP' 列。请先确保文件已按 BP 展开（每行一个 (Gene, BP) 对）.")
    gene_col = guess_gene_col(df)
    out = df[[gene_col, "BP"]].rename(columns={gene_col: "Gene"}).copy()
    out["Gene"] = out["Gene"].astype(str).str.strip()
    out["BP"] = out["BP"].astype(str).str.strip()
    out = out[(out["Gene"] != "") & (out["BP"] != "") & (out["BP"].str.lower() != "nan")]
    return out

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini–Hochberg FDR correction (monotone).
    Returns an array of q-values aligned with p-values.
    """
    m = len(pvals)
    order = np.argsort(pvals)
    ranked = np.empty(m, dtype=float)
    ranked[order] = np.arange(1, m + 1)
    bh = pvals * m / ranked
    # enforce monotonicity from the end
    bh_sorted = bh[order]
    bh_sorted = np.minimum.accumulate(bh_sorted[::-1])[::-1]
    q = np.empty(m, dtype=float)
    q[order] = np.clip(bh_sorted, 0, 1)
    return q

# -------- load & clean --------
c1_df = prep_df(c1_path)
univ_df = prep_df(universe_path)

# universe gene set
U_genes = set(univ_df["Gene"].unique())

# warn/ignore C1 genes not in universe
extra_genes = set(c1_df["Gene"].unique()) - U_genes
if extra_genes:
    print(f"注意：C1 中有 {len(extra_genes)} 个基因不在背景集中，将被忽略。")

# filter C1 to universe
c1_df = c1_df[c1_df["Gene"].isin(U_genes)].copy()

# mappings
bp2genes_univ = univ_df.groupby("BP")["Gene"].apply(set).to_dict()
bp2genes_c1 = c1_df.groupby("BP")["Gene"].apply(set).to_dict()

# sizes
M = len(U_genes)                              # background size
N = len(set(c1_df["Gene"].unique()))          # C1 size

# -------- enrichment per BP --------
records = []
for bp, genes_in_univ in bp2genes_univ.items():
    K = len(genes_in_univ)                    # bg genes annotated to this BP
    genes_in_c1 = bp2genes_c1.get(bp, set())  # C1 genes annotated to this BP
    x = len(genes_in_c1)                      # overlap

    # 2x2 table
    #            BP       not BP
    # C1         x        N - x
    # not C1     K - x    M - K - (N - x)
    a = x
    b = N - x
    c = K - x
    d = M - K - b

    # skip inconsistent rows
    if min(a, b, c, d) < 0:
        continue

    oddsratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")

    records.append({
        "term": bp,
        "overlap": x,
        "c1_size": N,
        "bg_with_term": K,
        "bg_size": M,
        "odds_ratio": oddsratio,
        "p_value": p,
        "gene_ratio": f"{x}/{N}" if N > 0 else "0/0",
        "bg_ratio": f"{K}/{M}" if M > 0 else "0/0",
        "genes_in_overlap": ";".join(sorted(genes_in_c1)),
    })

res = pd.DataFrame.from_records(records)

# -------- multiple testing correction --------
if not res.empty:
    # BH FDR
    res = res.sort_values("p_value", ascending=True).reset_index(drop=True)
    res["fdr_bh"] = bh_fdr(res["p_value"].to_numpy())

    # convenience columns
    res["neg_log10_p"] = -np.log10(res["p_value"].replace(0, np.nextafter(0, 1)))
    res["neg_log10_fdr"] = -np.log10(res["fdr_bh"].replace(0, np.nextafter(0, 1)))

    # final sort
    res = res.sort_values(["fdr_bh", "p_value", "overlap"], ascending=[True, True, False])

# -------- save --------
out_path.parent.mkdir(parents=True, exist_ok=True)
res.to_csv(out_path, index=False)

# -------- report --------
print(f"C1 基因数（去重且位于背景中）：{N}")
print(f"背景基因数（去重）：{M}")
print(f"测试到的 BP 术语数：{len(res)}")
print(f"结果已保存：{out_path}")
if extra_genes:
    print(f"提示：忽略了 {len(extra_genes)} 个不在背景中的 C1 基因。")
