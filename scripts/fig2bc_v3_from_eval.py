#!/usr/bin/env python3
"""Fig2b/c driven directly by benchmark_v3_mettd_eval.csv.

MetTD's discriminant is taken from that file's Target_rank column: how well the
record's LABELED target ranked among the 1,126 enzymes (rank percentile, smaller =
better). Records whose labelled target is not an enzyme in the network (no rank) are
dropped, so all five methods are evaluated on the same records.
ML methods: leave-one-study-out CV on the metabolite delta matrix, out-of-fold
probabilities pooled.
"""
import sys, json, warnings
sys.stdout.reconfigure(encoding='utf-8'); warnings.filterwarnings('ignore')
import numpy as np, pandas as pd

src = open('D:/MetTD_recompute/fig2bc_v3_losoc.py', encoding='utf-8').read()
ns = {}
exec(compile(src.split('from sklearn.preprocessing import StandardScaler')[0],
             'fig2bc_v3_losoc.py(data)', 'exec'), ns)
globals().update(ns)

ev = pd.read_csv('D:/MetTD_recompute/benchmark_v3_mettd_eval.csv', encoding='utf-8-sig')
assert len(ev) == len(v3) and (ev.Program.values == v3.Program.values).all()
rank = pd.to_numeric(ev.Target_rank, errors='coerce').to_numpy(float)
keep = np.isfinite(rank)
pct = np.where(keep, rank / nE, np.nan)

y = y[keep]; METk = MET[keep]; grpk = np.array(grp)[keep]; pctk = pct[keep]
rankk = rank[keep]
print(f'records kept {keep.sum()}/{len(keep)} (dropped {int((~keep).sum())} whose labelled '
      f'target is not in the enzyme network) | IDH1 {int(y.sum())} vs rest {int((1-y).sum())}')

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, \
    recall_score, f1_score, confusion_matrix
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from xgboost import XGBClassifier

logo = list(LeaveOneGroupOut().split(METk, y, grpk))
MK = {
 'XGBoost': lambda: XGBClassifier(n_estimators=35, learning_rate=0.1, max_depth=6,
                                  subsample=0.8, colsample_bytree=0.8,
                                  eval_metric='logloss', random_state=123),
 'SVM': lambda: SVC(kernel='sigmoid', gamma=5, probability=True, random_state=1),
 'BP': lambda: MLPClassifier(hidden_layer_sizes=(64, 32), max_iter=1000, random_state=1),
 'RF': lambda: RandomForestClassifier(n_estimators=100, max_features='sqrt', random_state=123)}

def metrics(lab):
    tn, fp, fn, tp = confusion_matrix(y, lab, labels=[0, 1]).ravel()
    return dict(acc=float(accuracy_score(y, lab)),
                prec=float(precision_score(y, lab, zero_division=0)),
                rec=float(recall_score(y, lab)), f1=float(f1_score(y, lab)),
                spec=float(tn / (tn + fp)))

P = {'MetTD': -pctk}
for k, f in MK.items():
    p = np.zeros(len(y))
    for tr, te in logo:
        sc = StandardScaler().fit(METk[tr])
        m = f(); m.fit(sc.transform(METk[tr]), y[tr])
        p[te] = m.predict_proba(sc.transform(METk[te]))[:, 1]
    P[k] = p

out = {'roc': {}, 'metrics': {}, 'protocol': {
    'mettd': 'Target_rank from benchmark_v3_mettd_eval.csv; no training, no split; '
             'Fig2c call = labelled target within rank 20 of the 1,126 enzymes',
    'ml': f'leave-one-study-out CV, {len(logo)} folds, out-of-fold probabilities',
    'n_records': int(len(y)), 'n_pos_IDH1': int(y.sum()), 'n_studies': int(len(logo))}}
TOPK = 20
LAB = (rankk <= TOPK).astype(int)
print(f'MetTD operating point (Fig2c bars): labelled target rank <= {TOPK} | '
      f'tP={int(((LAB==1)&(y==1)).sum())} fp={int(((LAB==1)&(y==0)).sum())} '
      f'fn={int(((LAB==0)&(y==1)).sum())} tn={int(((LAB==0)&(y==0)).sum())}')
for k, p in P.items():
    fpr, tpr, _ = roc_curve(y, p)
    lab = LAB if k == 'MetTD' else (p >= 0.5).astype(int)
    mm = metrics(lab)
    out['roc'][k] = {'fpr': fpr.tolist(), 'tpr': tpr.tolist()}
    out['metrics'][k] = {'auc': float(roc_auc_score(y, p)), **mm}
    print(f"{k:8s} AUC={roc_auc_score(y,p):.3f} Acc={mm['acc']:.2f} Prec={mm['prec']:.2f} "
          f"Rec={mm['rec']:.2f} F1={mm['f1']:.2f} Spec={mm['spec']:.2f}")

json.dump(out, open('D:/MetTD_recompute/fig2bc_v3_evalcsv.json', 'w'))
print('wrote fig2bc_v3_evalcsv.json')
