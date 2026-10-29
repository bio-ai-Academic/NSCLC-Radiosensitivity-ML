from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from workflow import *

x,y=load(); out=output('candidates')
l=pd.read_csv(ROOT/'analysis_outputs/lasso_SF2_5000/gene_scores.csv').set_index('gene')
m=pd.read_csv(ROOT/'analysis_outputs/mlp_SF2_5000/gene_scores.csv').set_index('gene')
scores=l.join(m,how='outer').fillna(0)
def z(v):
    return (v-v.mean())/v.std(ddof=0) if v.std(ddof=0)>0 else v*0
scores['fusion_score']=.5*z(scores.selection_frequency)+.5*z(scores.mean_abs_shap)
scores=scores.sort_values('fusion_score',ascending=False)
scores.to_csv(out/'all_gene_scores.csv')
candidate=scores.head(100).copy(); candidate.to_csv(out/'ranked_candidates.csv')
correlation=x[candidate.index].corr().fillna(0)
d=1-correlation.abs(); np.fill_diagonal(d.values,0)
candidate['cluster']=fcluster(linkage(squareform(d.values,checks=False),method='complete'),.1,criterion='distance')
candidate.groupby('cluster',sort=False).head(1).to_csv(out/'representative_subset.csv')
correlation.to_csv(out/'candidate_correlations.csv')
def pairs(labels):
    return {(i,j) for i in range(len(labels)) for j in range(i+1,len(labels)) if labels[i]==labels[j]}
base=pairs(candidate.cluster.to_numpy()); rng=np.random.default_rng(SEED); stability=[]
for iteration in range(100):
    xx=x[candidate.index].iloc[rng.integers(0,len(x),len(x))]
    dd=1-xx.corr().abs().fillna(0); np.fill_diagonal(dd.values,0)
    labels=fcluster(linkage(squareform(dd.values,checks=False),method='complete'),.1,criterion='distance')
    current=pairs(labels); union=base|current
    stability.append(dict(iteration=iteration,pairwise_cluster_jaccard=len(base&current)/len(union) if union else 1.))
pd.DataFrame(stability).to_csv(out/'cluster_bootstrap.csv',index=False)
a=set(l.nlargest(100,'selection_frequency').index); b=set(m.nlargest(100,'mean_abs_shap').index)
observed=len(a&b); universe=np.array(scores.index); null=[]
for _ in range(10000): null.append(len(a & set(rng.choice(universe,len(b),replace=False))))
(out/'model_overlap.json').write_text(json.dumps(dict(intersection=observed,jaccard=observed/len(a|b),
    random_set_p=(1+sum(v>=observed for v in null))/10001),indent=2))
