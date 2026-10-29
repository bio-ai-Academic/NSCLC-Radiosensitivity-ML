import argparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from workflow import *

def evaluate(x,y,k):
    pred=np.zeros(len(y))
    for fold,(train,test) in enumerate(folds(len(y))):
        fit=lasso_search(x[train],y[train],k,SEED+fold)
        pred[test]=fit.predict(x[test])
    return metrics(y,pred)

p=argparse.ArgumentParser(); p.add_argument('--permutations',type=int,default=100); a=p.parse_args()
out=output('robustness'); x,y=load(); rows=[]
for target in ['SF2','AUC_recomputed']:
    xx,yy=load(target)
    for k in [3000,5000,10000]: rows.append(dict(target=target,genes=k,**evaluate(xx.to_numpy(),yy,k)))
pd.DataFrame(rows).to_csv(out/'lasso_dimension_and_endpoint.csv',index=False)
rng=np.random.default_rng(SEED); observed=evaluate(x.to_numpy(),y,5000); null=[]
for iteration in range(a.permutations):
    permuted=rng.permutation(y)
    null.append(dict(iteration=iteration,**evaluate(x.to_numpy(),permuted,5000)))
pd.DataFrame(null).to_csv(out/'label_permutation.csv',index=False)
(out/'permutation_summary.json').write_text(json.dumps(dict(observed_R2=observed['R2'],
    p=(1+sum(v['R2']>=observed['R2'] for v in null))/(a.permutations+1)),indent=2))
external=ROOT/'analysis_inputs/external_expression_log2.csv'
if external.exists():
    xx=pd.read_csv(external,index_col=0).T
    yy=pd.read_csv(ROOT/'analysis_inputs/external_response.csv',index_col='sample_id')
    ids=xx.index.intersection(yy.index); xx=xx.loc[ids]; yy=yy.loc[ids,'SF2']; rows=[]
    for label,file in [('ranked','ranked_candidates.csv'),('representative','representative_subset.csv')]:
        genes=pd.read_csv(ROOT/'analysis_outputs/candidates'/file)['gene']
        for gene in genes:
            if gene not in xx.columns: continue
            valid=np.isfinite(xx[gene]) & np.isfinite(yy)
            if valid.sum()<5 or xx.loc[valid,gene].std()==0: continue
            r,pv=spearmanr(xx.loc[valid,gene],yy.loc[valid])
            rows.append(dict(set=label,gene=gene,n=int(valid.sum()),rho=r,p=pv))
    frame=pd.DataFrame(rows)
    if len(frame):
        frame['FDR']=np.nan
        for label,part in frame.groupby('set'): frame.loc[part.index,'FDR']=multipletests(part.p,method='fdr_bh')[1]
        frame.to_csv(out/'external_gene_associations.csv',index=False)
else:
    (out/'external_status.txt').write_text('Independent measured response inputs absent; external validation not performed.')
