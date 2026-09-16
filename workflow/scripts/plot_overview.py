"""Four-panel composition figure for the finished database.

Standalone CLI: usable outside Snakemake.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TIER1, TIER2 = "tier1_host_and_related", "tier2_environmental_diet_microbiome"
C1, C2 = "#1f5fa9", "#b8b2a7"
STATUS_ORDER = ["complete genome", "coding-complete", "near-complete",
                "substantial fragment", "gene fragment", "unscored"]


def italicise(ax, labels):
    ax.set_yticklabels(labels)
    for tick, lab in zip(ax.get_yticklabels(), labels):
        if lab[:1].isupper() and "(" not in lab:
            tick.set_fontstyle("italic")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reps", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dataset", required=True)
    args = ap.parse_args()

    plt.rcParams.update({
        "savefig.dpi": 300, "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlelocation": "left", "axes.titlepad": 6})

    reps = pd.read_csv(args.reps, sep="\t")
    records = pd.read_csv(args.records, sep="\t")
    for t in (TIER1, TIER2):
        if t not in set(reps.tier):
            reps.loc[len(reps)] = pd.NA      # keeps both tiers present in group-bys
            reps.loc[len(reps) - 1, ["tier", "length_bp"]] = [t, np.nan]

    hosts = (records.host_scientific_name.fillna("not resolved to species")
             .replace("", "not resolved to species").value_counts().sort_values())
    fam = (reps[reps.virus_family.notna() & (reps.virus_family != "unclassified")]
           .groupby(["virus_family", "tier"]).size().unstack(fill_value=0))
    for t in (TIER1, TIER2):
        if t not in fam:
            fam[t] = 0
    fam["total"] = fam[[TIER1, TIER2]].sum(axis=1)
    fam = fam.nlargest(12, "total").sort_values("total")
    status = (reps.groupby(["genome_status", "tier"]).size().unstack(fill_value=0))
    for t in (TIER1, TIER2):
        if t not in status:
            status[t] = 0
    status = status.reindex([s for s in STATUS_ORDER if s in status.index]).fillna(0)

    fig, ((a, b), (c, d)) = plt.subplots(2, 2, figsize=(11, 8.2))

    a.hlines(range(len(hosts)), 0.85, hosts.values, color="0.75", lw=0.7)
    a.plot(hosts.values, range(len(hosts)), "o", color=C1, ms=5.5)
    a.set_yticks(range(len(hosts)))
    italicise(a, list(hosts.index))
    a.set_xscale("log")
    a.set_xlim(0.85, max(hosts.values) * 4)
    a.set_xlabel("GenBank records (log scale)")
    a.set_title("Records per host species")
    for i, v in enumerate(hosts.values):
        a.text(v * 1.35, i, str(v), va="center", fontsize=6)
    a.margins(y=0.03)

    if len(fam):
        b.barh(range(len(fam)), fam[TIER1].values, color=C1, height=0.72,
               label="Host related")
        b.barh(range(len(fam)), fam[TIER2].values, left=fam[TIER1].values, color=C2, height=0.72,
               label="Environmental, diet or microbiome")
        b.set_yticks(range(len(fam)))
        italicise(b, list(fam.index))
        b.legend(frameon=False, loc="lower right")
    b.set_xlabel("Representative genomes (n)")
    b.set_title("Viral families in the reference set")
    b.margins(x=0.06, y=0.03)

    y = np.arange(len(status))
    c.barh(y - 0.19, status[TIER1].values, height=0.34, color=C1, label="tier 1")
    c.barh(y + 0.19, status[TIER2].values, height=0.34, color=C2, label="tier 2")
    c.set_yticks(y)
    c.set_yticklabels([s.replace(" ", "\n") if len(s) > 16 else s for s in status.index])
    c.invert_yaxis()
    c.set_xlabel("Representative genomes")
    t1 = reps[reps.tier == TIER1]
    usable = ["complete genome", "coding-complete", "near-complete"]
    c.set_title("%d of %d tier-1 entries are genome-scale"
                % (int(t1.genome_status.isin(usable).sum()), len(t1)))
    c.legend(frameon=False, loc="lower right")
    c.margins(x=0.06)

    lengths = reps.length_bp.dropna()
    lo, hi = max(lengths.min(), 100), lengths.max()
    bins = np.logspace(np.log10(lo * 0.9), np.log10(hi * 1.1), 34)
    d.hist([reps.loc[reps.tier == TIER1, "length_bp"].dropna(),
            reps.loc[reps.tier == TIER2, "length_bp"].dropna()],
           bins=bins, stacked=True, color=[C1, C2], label=["Host related", "Environmental, diet or microbiome"])
    d.set_xscale("log")
    d.set_xlabel("Sequence length (bp, log scale)")
    d.set_ylabel("Representative genomes", labelpad=2)
    d.set_title("Sequence length spans %.1f kb to %.0f kb" % (lo / 1000, hi / 1000))
    d.legend(frameon=False, loc="upper left")
    d.margins(x=0.03)

    for ax, letter in zip((a, b, c, d), "abcd"):
        ax.text(-0.08, 1.06, letter, transform=ax.transAxes, fontweight="bold", fontsize=10, va="top")
    fig.suptitle("%s reference database -- composition (NCBI GenBank)" % args.dataset,
                 x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(args.out)
    print("wrote " + args.out)


if __name__ == "__main__":
    main()
