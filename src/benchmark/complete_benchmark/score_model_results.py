"""

"""
import argparse
import os
from pathlib import Path
import json
import pandas as pd
from sklearn.metrics import f1_score, accuracy_score, classification_report

CODING = {
    "3FE4S": (1,0,0),
    "4FE4S": (0,1,0),
    "4FE3S": (0,0,1),
    "n": (0,0,0)
}

def _keep_only_predictions(per_fes_gt, per_fes_pred):
    res_gt = []
    res_pred = []
    for i in range(len(per_fes_pred)):
        if per_fes_pred[i] == (0,0,0):
            continue
        else:
            res_gt.append(per_fes_gt[i])
            res_pred.append(per_fes_pred[i])

    return res_gt, res_pred

def score_model_results(ground_truth : Path, gt_index_col : str, prediction : Path, prediction_index_col : str):
    """
    """
    gt = pd.read_csv(ground_truth, sep="\t")
    pred = pd.read_csv(prediction, sep="\t")
    sum = 0

    wrong = []
    unpredicted = 0

    per_fes_gt = []
    per_fes_pred = []

    proximal_gt = []
    proximal_pred = []

    medial_gt = []
    medial_pred = []
    
    distal_gt = []
    distal_pred = []    

    for index, row in pred.iterrows():    
        score = 0.0
        gt_row = gt[gt[gt_index_col] == row[prediction_index_col]]

        prox_true = CODING[gt_row["proximal_formula"].values[0]]
        prox_pred = CODING[row["proximal"]]

        if prox_true == prox_pred:
            score += 1

        proximal_gt.append(prox_true)
        proximal_pred.append(prox_pred)
        per_fes_gt.append(prox_true)
        per_fes_pred.append(prox_pred)

        med_true = CODING[gt_row["medial_formula"].values[0]]
        med_pred = CODING[row["medial"]]

        if med_true == med_pred:
            score += 1

        medial_gt.append(med_true)
        medial_pred.append(med_pred)
        per_fes_gt.append(med_true)
        per_fes_pred.append(med_pred)

        dist_true = CODING[gt_row["distal_formula"].values[0]]
        dist_pred = CODING[row["distal"]]

        if dist_true == dist_pred:
            score += 1

        distal_gt.append(dist_true)
        distal_pred.append(dist_pred)
        per_fes_gt.append(dist_true)
        per_fes_pred.append(dist_pred)
        

        score = score / 3

        if score != 1.0:

            wrong_string = row[prediction_index_col]
            if gt_row["proximal_formula"].values[0] != row["proximal"]:
                wrong_string += f" p:{row['proximal']}|{gt_row['proximal_formula'].values[0]}"

            if gt_row["medial_formula"].values[0] != row["medial"]:
                wrong_string += f" m:{row['medial']}|{gt_row['medial_formula'].values[0]}"

            if gt_row["distal_formula"].values[0] != row["distal"]:
                wrong_string += f" d:{row['distal']}|{gt_row['distal_formula'].values[0]}"
            wrong.append(wrong_string)

        if row["proximal"] == "n" and row["medial"] == "n" and row["distal"] == "n":
            unpredicted += 1

        print(f"{row[prediction_index_col]} : {score}")
        sum += score

    print("="*20, "SCORE", "="*20)
    print("PER SAMPLE")
    print("Overall Acc.:\t", sum/pred.shape[0])
    print("Predicted Acc.:\t", sum/(pred.shape[0] - unpredicted))
    print("Coverage:\t", (pred.shape[0]-unpredicted)/pred.shape[0])
    print()

    print("PER POCKET")
    print(f"All pockets:\n- acc {accuracy_score(per_fes_gt, per_fes_pred)}\n- ", classification_report(per_fes_gt, per_fes_pred, zero_division=0), "\n")

    per_fes_gt_only_predictions, per_fes_pred_only_predictions = _keep_only_predictions(per_fes_gt, per_fes_pred)

    print(f"All pockets (only predicted):\n- acc {accuracy_score(per_fes_gt_only_predictions, per_fes_pred_only_predictions)}")


    print(f"Proximal:\n- acc {accuracy_score(proximal_gt, proximal_pred)}\n- ", classification_report(proximal_gt, proximal_pred, zero_division=0.0))
    print(f"Medial:\n- acc {accuracy_score(medial_gt, medial_pred)}\n- ", classification_report(medial_gt, medial_pred, zero_division=0.0))
    print(f"Distal:\n- acc {accuracy_score(distal_gt, distal_pred)}\n- ", classification_report(distal_gt, distal_pred, zero_division=0.0))

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