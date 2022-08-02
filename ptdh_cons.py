import os
from collections import Counter
from typing import Dict, List

import numpy as np
import pandas as pd
from prody import *
from scipy.stats import norm, zscore

POS_MAP = {
    "Alpha1": [154, 155, 156, 157],
    "Arch": [101, 104, 100, 97],
    "Base loop": [207, 208, 209, 210],
    "Pocket": [235, 236, 237, 261, 262, 292, 293, 294, 295],
    "Seatbelt loop": [76, 77, 78, 79],
}


def clean_aln(fasta_path: str, mafft_path: str, label: str, refined_path: str) -> None:

    # run mafft
    cmd = f"mafft {fasta_path} > {mafft_path}"
    os.system(cmd)

    # degap
    msa = parseMSA(mafft_path)
    refined_msa = refineMSA(msa, label=label, seqid=0.98)
    writeMSA(refined_path, refined_msa, format="fasta")


def get_entropy_df(msa: MSA, window: int = 6) -> pd.DataFrame:

    # per residue stats
    entropy = calcShannonEntropy(msa)
    entropy_df = pd.DataFrame({"se": entropy})
    entropy_df["z"] = zscore(entropy_df["se"])
    entropy_df["ppf"] = norm.cdf(entropy_df["z"])
    entropy_df["conservation_score"] = entropy_df["ppf"].max() - entropy_df["ppf"]
    entropy_df["rolling"] = entropy_df["conservation_score"].rolling(window).mean()

    pos = np.arange(1, len(entropy_df) + 1)
    entropy_df.insert(0, "pos", pos)

    return entropy_df


def get_pos_cons(
    entropy_df: pd.DataFrame, pos_map: Dict[str, List[str]]
) -> pd.DataFrame:

    # select regions neighboring cofactor
    save = []
    for region in pos_map.keys():
        positions = pos_map[region]
        for position in positions:
            pos_data = entropy_df.query("pos == @position")
            pos_data = pos_data.assign(region=region)
            save.append(pos_data)

    pos_cons_df = pd.concat(save)
    return pos_cons_df.sort_values(by="pos")


def group_cons(pos_cons_df: pd.DataFrame) -> pd.DataFrame:

    # aggregate per-residue stats for each region
    pos_grp_df = (
        pos_cons_df.groupby("region").agg({"rolling": ["mean", "std"]}).reset_index()
    )
    pos_grp_df.columns = ["_".join(col).rstrip("_") for col in pos_grp_df.columns]

    return pos_grp_df


def get_res_counts(msa: MSA, resi: int) -> Counter:

    counter = Counter()
    residues = msa[:, resi - 1 : resi]
    save = []
    for residue in residues:
        counter.update(str(residue))

    return counter


def get_all_counts(msa: MSA, pos_map: Dict[str, List[int]]) -> pd.DataFrame:

    # selected residue positions
    positions = []
    for k, v in pos_map.items():
        positions.extend(v)
    positions = sorted(positions)

    # count freq of amino acids at each position
    counts_list = []
    for position in positions:
        counts = get_res_counts(msa, position)
        counts_list.append(pd.Series(counts))

    counts_df = pd.concat(counts_list, axis=1)
    counts_df.columns = positions

    new_index = [
        "A", "V", "L", "I", "M", "C",
        "F", "Y", "W", "H", "S", "T",
        "N", "Q", "D", "E", "K", "R",
        "G", "P", "X", "-",
    ]

    counts_df = counts_df.reindex(new_index)
    counts_df = counts_df.sort_index(axis=1)
    counts_df = counts_df.fillna(0)

    # drop gap (-) and unknown (X)
    counts_df = counts_df.drop(["-", "X"], axis=0)

    # normalize
    total_sum = counts_df.sum(axis=0)

    return (counts_df / total_sum) * 100


if __name__ == "__main__":

    # degapped msa from blast hits
    clean_aln(
        "./4e5n_blast.fasta", "./4e5n_mafft.fasta", "4E5N", "./4e5n_refined.fasta"
    )
    msa = parseMSA("./4e5n_refined.fasta")

    # s2b
    entropy_df = get_entropy_df(msa)
    entropy_df.to_csv("./ptdh_cons_s2b.csv", index=False)

    # s2c
    pos_cons_df = get_pos_cons(entropy_df, POS_MAP)
    pos_cons_df.to_csv("./ptdh_cons_s2c.csv", index=False)
    pos_grp_df = group_cons(pos_cons_df)
    pos_grp_df.to_csv("./ptdh_cons_s2c2.csv", index=False)

    # s3
    counts_df = get_all_counts(msa, POS_MAP)
    counts_df.to_csv("./ptdh_cons_s3.csv")
