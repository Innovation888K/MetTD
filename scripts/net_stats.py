"""B5: network construction statistics (numpy version)."""
import csv, sys, os
import numpy as np
from collections import Counter, defaultdict

os.chdir(r'D:/MetTD_recompute')
sys.stdout.reconfigure(encoding='utf-8')
N = 9492
nodes = [None] * (N + 1)
for r in csv.DictReader(open('nodes.csv')):
    nodes[int(float(r['idx']))] = r['node']
enz_idx = set(int(r['idx']) for r in csv.DictReader(open('enzyme.csv')))

M = np.memmap('model_dist.bin', dtype=np.float64, mode='r').reshape(N, N)  # col-major: [i,j] = distance(i+1,j+1)
# note: file was written column-major via as.vector (column-stacking), so reshape(N,N) gives element [row, col] directly
A = (M == 2.0)
np.fill_diagonal(A, False)
deg = A.sum(axis=0) + A.sum(axis=1) - (A.sum(axis=0) & A.sum(axis=1))  # symmetric; just use one side
deg = A.sum(axis=0).astype(np.int64)  # column sums = degree (symmetric matrix)
edges = int(A.sum()) // 2

parent = list(range(N))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x

ii, jj = np.nonzero(np.triu(A, 1))
for a, b in zip(ii.tolist(), jj.tolist()):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb
comp = Counter(find(i) for i in range(N))
sizes = sorted(comp.values(), reverse=True)

enz_mask = np.zeros(N, bool)
for i in enz_idx: enz_mask[i - 1] = True
print('nodes:', N, ' metabolites:', N - len(enz_idx), ' enzymes:', len(enz_idx))
print('edges (enzyme-metabolite):', edges)
print('connected components:', len(comp), ' largest sizes:', sizes[:10])
print('isolated nodes:', sum(1 for s in sizes if s == 1))
dm = deg[~enz_mask]; de = deg[enz_mask]
print('metabolite degree mean/median/max: %.1f / %d / %d' % (dm.mean(), np.median(dm), dm.max()))
print('enzyme degree mean/median/max: %.1f / %d / %d' % (de.mean(), np.median(de), de.max()))
order = np.argsort(-deg)[:15]
print('top-15 hubs:')
for o in order:
    print('   deg=%d %s %s' % (deg[o], nodes[o + 1], '(enzyme)' if (o + 1) in enz_idx else ''))

CURRENCY = {'atp','adp','amp','h2o','h+','nad+','nadh','nadp+','nadph','nadp','pi','coa','co2','nh3','nh4+','gtp','gdp','udp','utp','fad','fadh2','nad','nadh'}
cur_nodes = [i for i in range(1, N + 1) if i not in enz_idx and nodes[i].rsplit('.', 1)[0].lower() in CURRENCY]
cur_deg = int(sum(deg[i - 1] for i in cur_nodes))
print('currency metabolite nodes:', len(cur_nodes), ' degree share of all edge ends: %.1f%%' % (100 * cur_deg / (2 * edges)))
print('examples:', [nodes[i] for i in cur_nodes[:10]])

cc = Counter(nodes[i].rsplit('.', 1)[-1] for i in range(1, N + 1) if i not in enz_idx)
print('metabolite compartments:', cc.most_common(12))
namec = Counter(nodes[i].rsplit('.', 1)[0] for i in range(1, N + 1) if i not in enz_idx)
print('names in >1 compartment:', sum(1 for c in namec.values() if c > 1))

# degree distribution buckets for SI
hist, edgesb = np.histogram(dm, bins=[1,2,3,5,10,20,50,100,200,500,10**6])
print('metabolite degree histogram [1,2,3,5,10,20,50,100,200,500+]:', hist.tolist())
hist2, _ = np.histogram(de, bins=[1,2,3,5,10,20,50,100,200,500,10**6])
print('enzyme degree histogram:', hist2.tolist())
