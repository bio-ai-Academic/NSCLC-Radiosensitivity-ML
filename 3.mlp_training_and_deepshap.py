import argparse
import tensorflow as tf
import shap
from hyperopt import fmin,tpe,hp,Trials,STATUS_OK
from sklearn.model_selection import train_test_split
from workflow import *

def run(target='SF2',k=5000,trials=20):
    tf.keras.utils.set_random_seed(SEED)
    x,y=load(target); values=x.to_numpy(); out=output(f'mlp_{target}_{k}')
    predictions=np.zeros(len(y)); rows=[]; importance=[]; frequencies=[]
    for fold,(train,test) in enumerate(folds(len(y))):
        core,validation=train_test_split(train,test_size=.2,random_state=SEED+fold)
        select=TopVariance(k).fit(values[core]); scale=StandardScaler().fit(select.transform(values[core]))
        a=scale.transform(select.transform(values[core])); b=scale.transform(select.transform(values[validation]))
        ym,ys=y[core].mean(),y[core].std()
        if ys==0: raise ValueError('Constant training outcome')
        def objective(params):
            tf.keras.backend.clear_session(); tf.keras.utils.set_random_seed(SEED+fold)
            net=build_mlp(a.shape[1],**params)
            history=net.fit(a,(y[core]-ym)/ys,validation_data=(b,(y[validation]-ym)/ys),
                epochs=300,batch_size=16,callbacks=callbacks(),verbose=0)
            return dict(loss=float(min(history.history['val_loss'])),status=STATUS_OK,
                        epochs=int(np.argmin(history.history['val_loss'])+1))
        space=dict(lr=hp.loguniform('lr',np.log(1e-4),np.log(3e-3)),
            dropout=hp.uniform('dropout',.2,.5),l2=hp.loguniform('l2',np.log(1e-5),np.log(1e-3)))
        history=Trials()
        best=fmin(objective,space,algo=tpe.suggest,max_evals=trials,trials=history,rstate=np.random.default_rng(SEED+fold))
        epochs=history.best_trial['result']['epochs']
        select=TopVariance(k).fit(values[train]); scale=StandardScaler().fit(select.transform(values[train]))
        a=scale.transform(select.transform(values[train])).astype('float32')
        b=scale.transform(select.transform(values[test])).astype('float32')
        ym,ys=y[train].mean(),y[train].std()
        tf.keras.backend.clear_session(); tf.keras.utils.set_random_seed(SEED+fold)
        net=build_mlp(a.shape[1],**best)
        net.fit(a,(y[train]-ym)/ys,epochs=epochs,batch_size=16,verbose=0)
        pred=net.predict(b,verbose=0).ravel()*ys+ym; predictions[test]=pred
        rows.append(dict(fold=fold,**metrics(y[test],pred)))
        net.save(out/f'fold_{fold}.keras')
        joblib.dump(dict(select=select,scale=scale,y_mean=ym,y_sd=ys,train=train,test=test),out/f'fold_{fold}_preprocessing.joblib')
        (out/f'fold_{fold}_parameters.json').write_text(json.dumps(dict(**best,epochs=epochs),indent=2))
        explainer=shap.DeepExplainer(net,a[:min(30,len(a))])
        sv=explainer.shap_values(b)
        if isinstance(sv,list): sv=sv[0]
        sv=np.asarray(sv).reshape(len(test),a.shape[1])*ys
        pd.DataFrame(sv,index=x.index[test],columns=x.columns[select.indices_]).to_csv(out/f'fold_{fold}_shap.csv')
        full=np.zeros(x.shape[1]); full[select.indices_]=np.abs(sv).mean(axis=0)
        importance.append(full); top=np.zeros(x.shape[1])
        top[select.indices_[np.argsort(np.abs(sv).mean(axis=0))[::-1][:100]]]=1
        frequencies.append(top)
    report(y,predictions,x.index,out,rows)
    pd.DataFrame(dict(gene=x.columns,mean_abs_shap=np.mean(importance,axis=0),
                     top100_frequency=np.mean(frequencies,axis=0))).to_csv(out/'gene_scores.csv',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--target',default='SF2'); p.add_argument('--genes',type=int,default=5000)
    p.add_argument('--trials',type=int,default=20); a=p.parse_args(); run(a.target,a.genes,a.trials)
