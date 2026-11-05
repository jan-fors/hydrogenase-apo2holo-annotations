"""
python -m apo2holo.database.analysis.analyze_fingerprint_tsv benchmarking/complete_benchmark/2026_09_01_16_16_experiment_04/0/db/fingerprintDB/dd_experiment_04_fpr7.0.tsv
"""

import argparse
from pathlib import Path
import pandas as pd
from collections import Counter
from pprint import pprint
import numpy as np
import seaborn as sns
from scipy.stats import entropy
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
from apo2holo.utils.constants import COFACTOR_WHITELIST
import matplotlib.pyplot as plt
import umap
from scipy.spatial.distance import cdist


from itertools import combinations
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

def analyze_fingerprint_tsv(fingerprint_tsv : Path):
    """
    """
    df = pd.read_csv(fingerprint_tsv, sep="\t")

    n_entries = df.shape[0]
    print(f"n_entries: {n_entries}")

    
    try:
        print("\nres_name distribution:")
        res_name_list = df["res_name"].to_list()
        counter = dict(Counter(res_name_list))

        negatives = 0
        for key in counter.keys():
            if key in COFACTOR_WHITELIST:
                print(f"{key} : {counter[key]}")
            else:
                negatives += 1
        print(f"Negatives: {negatives}")
    except:
        print("No Resname")

    print("\nformula distribution")
    formula_list = df["formula"].to_list()
    counter = dict(Counter(formula_list))

    pprint(counter)

    print("\ntype distribution:")
    type_list = df["type"].to_list()
    counter = dict(Counter(type_list))

    pprint(counter)
   
    #print(df.columns)


    outdir = Path("plots2")
    outdir.mkdir(exist_ok=True)

    X = df.drop(columns=["formula","aug_dist","structure","type","fingerprint"]).apply(pd.to_numeric, errors="coerce").fillna(0)

    # 1. Alpha diversity (richness and Shannon)
    div = pd.DataFrame({
        "class": df["formula"],
        "richness": (X > 0).sum(axis=1),
        "shannon": X.apply(lambda r: entropy(r[r > 0]) if r.sum() > 0 else 0, axis=1),
    })
    for metric in ["richness", "shannon"]:
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.violinplot(data=div, x="class", y=metric, inner="box", ax=ax)
        fig.savefig(outdir / f"alpha_{metric}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 2. Prevalence clustermap
    prev = (X > 0).groupby(df["formula"]).mean()
    prev = prev.loc[:, prev.max() > 0.05]
    g = sns.clustermap(prev, cmap="viridis", figsize=(14, 5))
    g.savefig(outdir / "prevalence_clustermap.png", dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    # 3. Jaccard embedding
    emb = umap.UMAP(metric="jaccard", random_state=0).fit_transform((X > 0).values)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(x=emb[:, 0], y=emb[:, 1], hue=df["formula"].values,
                    s=10, alpha=0.6, linewidth=0, ax=ax)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(title="class", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.savefig(outdir / "jaccard_umap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    

    B = (X > 0).values
    labels = df["formula"].values
    classes = np.unique(labels)

    rng = np.random.default_rng(0)
    max_n = 500
    idx = {c: rng.choice(np.flatnonzero(labels == c),
                        size=min(max_n, (labels == c).sum()), replace=False)
        for c in classes}

    M = pd.DataFrame(np.nan, index=classes, columns=classes)
    for i, a in enumerate(classes):
        for b in classes[i:]:
            d = cdist(B[idx[a]], B[idx[b]], metric="jaccard")
            if a == b:
                d = d[np.triu_indices_from(d, k=1)]   # skip self-comparisons
            M.loc[a, b] = M.loc[b, a] = d.mean()

    g = sns.clustermap(M, cmap="viridis", annot=True, fmt=".2f",
                    vmin=0, vmax=1, figsize=(8, 7),
                    cbar_kws={"label": "mean Jaccard distance"})
    g.savefig(outdir / "class_jaccard_distances.png", dpi=300, bbox_inches="tight")
    plt.close(g.fig)

    P = (X > 0).groupby(labels).mean()            # prevalence per class
    score = (P.max() - P.min()).sort_values(ascending=False)
    top = score.head(20).index

    fig, ax = plt.subplots(figsize=(1.2 * len(P) + 4, 8))
    sns.heatmap(P[top].T * 100, annot=True, fmt=".0f", cmap="Blues",
                vmin=0, vmax=100, cbar_kws={"label": "% of samples with feature"}, ax=ax)
    ax.set_yticklabels([f"{f}  ({score[f]:.2f})" for f in top], rotation=0)
    ax.set_xlabel("class")
    ax.set_ylabel("feature  (score)")
    ax.set_title("Features that differ most between classes")
    fig.savefig(outdir / "top_discriminating_features.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    score.to_csv(outdir / "feature_scores.tsv", sep="\t", header=["score"])

    from sklearn.tree import DecisionTreeClassifier, plot_tree, export_text
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import balanced_accuracy_score

    Xtr, Xte, ytr, yte = train_test_split(X, labels, stratify=labels,
                                        test_size=0.3, random_state=0)
    tree = DecisionTreeClassifier(max_depth=3, min_samples_leaf=10,
                                class_weight="balanced", random_state=0).fit(Xtr, ytr)
    acc = balanced_accuracy_score(yte, tree.predict(Xte))

    fig, ax = plt.subplots(figsize=(22, 10))
    plot_tree(tree, feature_names=X.columns.astype(str), class_names=tree.classes_.astype(str),
            filled=True, rounded=True, impurity=False, proportion=True, fontsize=9, ax=ax)
    ax.set_title(f"Decision tree (balanced accuracy on held-out data: {acc:.2f})")
    fig.savefig(outdir / "decision_tree.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    (outdir / "decision_tree_rules.txt").write_text(
        export_text(tree, feature_names=list(X.columns.astype(str)), show_weights=True))


    Pb = (X > 0).values
    classes = np.unique(labels)
    rows = []
    for a, b in combinations(classes, 2):
        ia, ib = labels == a, labels == b
        na, nb = ia.sum(), ib.sum()
        ca, cb = Pb[ia].sum(axis=0), Pb[ib].sum(axis=0)
        for j, f in enumerate(X.columns):
            if ca[j] + cb[j] == 0:          # feature absent in both classes
                continue
            _, p = fisher_exact([[ca[j], na - ca[j]], [cb[j], nb - cb[j]]])
            rows.append((f"{a} vs {b}", a, b, f, ca[j] / na, cb[j] / nb, p))

    res = pd.DataFrame(rows, columns=["pair", "class_a", "class_b", "feature", "prev_a", "prev_b", "p"])
    res["diff"] = res["prev_a"] - res["prev_b"]
    res["q"] = multipletests(res["p"], method="fdr_bh")[1]
    res["significant"] = (res["q"] < 0.05) & (res["diff"].abs() >= 0.1)
    res.sort_values("q").to_csv(outdir / "pairwise_tests.tsv", sep="\t", index=False)

    C = pd.DataFrame(0, index=classes, columns=classes)
    for (a, b), n in res[res["significant"]].groupby(["class_a", "class_b"]).size().items():
        C.loc[a, b] = C.loc[b, a] = n

    fig, ax = plt.subplots(figsize=(1.0 * len(classes) + 3, 1.0 * len(classes) + 2))
    sns.heatmap(C, annot=True, fmt="d", cmap="Blues", mask=np.eye(len(classes), dtype=bool),
                cbar_kws={"label": "number of significantly different features"}, ax=ax)
    ax.set_title("Significantly different features per class pair")
    fig.savefig(outdir / "pairwise_significant_counts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    sig = res[res["significant"]]
    top = sig.groupby("feature")["diff"].apply(lambda d: d.abs().max()) \
            .sort_values(ascending=False).head(20).index

    sub = res[res["feature"].isin(top)]
    D = sub.pivot(index="feature", columns="pair", values="diff").reindex(top).fillna(0) * 100
    S = sub.pivot(index="feature", columns="pair", values="significant").reindex(top).fillna(False)
    annot = D.round().astype(int).astype(str) + np.where(S.astype(bool), "*", "")

    fig, ax = plt.subplots(figsize=(1.0 * D.shape[1] + 4, 8))
    sns.heatmap(D, annot=annot, fmt="", cmap="RdBu_r", center=0, vmin=-100, vmax=100,
                cbar_kws={"label": "difference in % of samples (first − second class)"}, ax=ax)
    ax.set_xlabel("class pair")
    ax.set_ylabel("feature")
    ax.set_title("Largest significant differences between classes  (* = q < 0.05)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.savefig(outdir / "pairwise_feature_differences.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fingerprint_tsv",type=Path)
    args = parser.parse_args()

    analyze_fingerprint_tsv(args.fingerprint_tsv)