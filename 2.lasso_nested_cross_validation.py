import argparse
from sklearn.model_selection import LeaveOneOut
from workflow import *

def run(target='SF2',k=5000,stability=100,loocv=False):
    x,y=load(target); values=x.to_numpy(); out=output(f'lasso_{target}_{k}')
    predictions=np.zeros(len(y)); rows=[]
    for fold,(train,test) in enumerate(folds(len(y))):
        fit=lasso_search(values[train],y[train],k,SEED+fold)
        pred=fit.predict(values[test]); predictions[test]=pred
        rows.append(dict(fold=fold,alpha=fit.best_params_['model__alpha'],**metrics(y[test],pred)))
        joblib.dump(fit.best_estimator_,out/f'fold_{fold}.joblib')
        pd.DataFrame(fit.cv_results_).to_csv(out/f'fold_{fold}_tuning.csv',index=False)
    report(y,predictions,x.index,out,rows)
    final=lasso_search(values,y,k).best_estimator_
    joblib.dump(final,out/'final_model.joblib')
    scores=pd.DataFrame(dict(gene=x.columns,coefficient=0.,selection_frequency=0.))
    scores.loc[final['select'].indices_,'coefficient']=final['model'].coef_
    rng=np.random.default_rng(SEED)
    for iteration in range(stability):
        ii=rng.choice(len(y),int(.8*len(y)),replace=False)
        fit=lasso_search(values[ii],y[ii],k,SEED+iteration).best_estimator_
        chosen=fit['select'].indices_[np.abs(fit['model'].coef_)>1e-10]
        scores.loc[chosen,'selection_frequency']+=1/stability
    scores['stable_nonzero']=(scores.selection_frequency>=.6)&(scores.coefficient.abs()>1e-10)
    scores.to_csv(out/'gene_scores.csv',index=False)
    if loocv:
        pred=np.zeros(len(y))
        for train,test in LeaveOneOut().split(x):
            pred[test]=lasso_search(values[train],y[train],k).predict(values[test])
        pd.DataFrame(dict(sample_id=x.index,observed=y,predicted=pred)).to_csv(out/'loocv_predictions.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--target',default='SF2'); p.add_argument('--genes',type=int,default=5000)
    p.add_argument('--stability',type=int,default=100); p.add_argument('--loocv',action='store_true')
    a=p.parse_args(); run(a.target,a.genes,a.stability,a.loocv)
