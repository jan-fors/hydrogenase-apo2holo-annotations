"""

"""
import argparse
import os
from pathlib import Path
import json
import pandas as pd

def score_model_results(ground_truth : Path, gt_index_col : str, prediction : Path, prediction_index_col : str):
    """
    """
    gt = pd.read_csv(ground_truth, sep="\t")
    pred = pd.read_csv(prediction, sep="\t")
    sum = 0

    wrong = []

    for index, row in pred.iterrows():    
        score = 0.0
        gt_row = gt[gt[gt_index_col] == row[prediction_index_col]]

        gt_prox = gt_row["proximal_formula"].values[0]
        if gt_prox == row["proximal"]:
            score += 1

        gt_med = gt_row["medial_formula"].values[0]
        if gt_med == row["medial"]:
            score += 1

        gt_dist = gt_row["distal_formula"].values[0]
        if gt_dist == row["distal"]:
            score += 1

        score = score / 3

        if score != 1.0:
            wrong_string = row[prediction_index_col]
            if gt_prox != row["proximal"]:
                wrong_string += f" p:{row['proximal']}|{gt_prox}"

            gt_med = gt_row["medial_formula"].values[0]
            if gt_med != row["medial"]:
                wrong_string += f" m:{row['medial']}|{gt_med}"

            gt_dist = gt_row["distal_formula"].values[0]
            if gt_dist != row["distal"]:
                wrong_string += f" d:{row['distal']}|{gt_dist}"
            wrong.append(wrong_string)

        print(f"{row[prediction_index_col]} : {score}")
        sum += score

    print("="*20, "SCORE", "="*20)
    print(sum/pred.shape[0])
    print("="*20, "WRONG", "="*20)
    print("<position>:<prediction>|<ground_truth>")
    for w in wrong:
        print("-",w)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ground_truth", type=Path)
    parser.add_argument("gt_index_col", type=str)
    parser.add_argument("prediction", type=Path)
    parser.add_argument("prediction_index_col", type=str)
    args = parser.parse_args()

    score_model_results(args.ground_truth, args.gt_index_col, args.prediction, args.prediction_index_col)