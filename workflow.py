from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

ROOT = Path(__file__).resolve().parent
PREP = ROOT / 'analysis_outputs' / 'prepared'
SEED = 123

class TopVariance(BaseEstimator, TransformerMixin):
    def __init__(self, k=5000):
        self.k = k
    def fit(self, X, y=None):
        v = np.var(X, axis=0)
        valid = np.flatnonzero(v > 1e-12)
        self.indices_ = valid[np.argsort(v[valid])[::-1][:self.k]]
        if not len(self.indices_):
            raise ValueError('No variable genes in training partition')
        return self
    def transform(self, X):
        return np.asarray(X)[:, self.indices_]

def load(target='SF2'):
    x = pd.read_csv(PREP / 'expression_log2.csv', index_col=0).T
    response = pd.read_csv(PREP / 'response.csv', index_col='sample_id').reindex(x.index)
    y = response[target].to_numpy(float)
    valid = np.isfinite(y)
    return x.loc[valid], y[valid]

def output(name):
    path = ROOT / 'analysis_outputs' / name
    path.mkdir(parents=True, exist_ok=True)
    return path

def folds(n):
    return list(KFold(5, shuffle=True, random_state=SEED).split(np.arange(n)))

def metrics(y, pred):
    return dict(R2=float(r2_score(y,pred)), RMSE=float(np.sqrt(mean_squared_error(y,pred))),
                MAE=float(mean_absolute_error(y,pred)), Pearson_r=float(pearsonr(y,pred).statistic))

def report(y, pred, ids, out, rows):
    pd.DataFrame(dict(sample_id=ids,observed=y,predicted=pred,residual=y-pred)).to_csv(out/'oof_predictions.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'fold_metrics.csv',index=False)
    (out/'metrics.json').write_text(json.dumps(metrics(y,pred),indent=2))
    rng=np.random.default_rng(SEED); bootstrap=[]
    for _ in range(1000):
        ii=rng.integers(0,len(y),len(y))
        if np.std(y[ii])>0 and np.std(pred[ii])>0:
            bootstrap.append(metrics(y[ii],pred[ii]))
    pd.DataFrame(bootstrap).quantile([.025,.975]).to_csv(out/'oof_bootstrap_intervals.csv')
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    ax[0].scatter(y,pred); bounds=[min(y.min(),pred.min()),max(y.max(),pred.max())]
    ax[0].plot(bounds,bounds,'k--'); ax[0].set(xlabel='Observed response',ylabel='Out-of-fold prediction')
    ax[1].hist(y-pred,bins='auto'); ax[1].set(xlabel='Residual',ylabel='Count')
    fig.tight_layout(); fig.savefig(out/'diagnostics.pdf'); plt.close(fig)

def lasso_search(x,y,k=5000,seed=SEED):
    pipe=Pipeline([('select',TopVariance(k)),('scale',StandardScaler()),('model',Lasso(max_iter=50000))])
    search=GridSearchCV(pipe,{'model__alpha':np.logspace(-5,0,40)},
        cv=KFold(5,shuffle=True,random_state=seed),scoring='neg_mean_squared_error',n_jobs=-1)
    return search.fit(x,y)

def build_mlp(d,lr=.001,dropout=.3,l2=.0001):
    import tensorflow as tf
    reg=tf.keras.regularizers.l2(l2)
    model=tf.keras.Sequential([tf.keras.Input(shape=(d,)),
        tf.keras.layers.Dense(128,activation='relu',kernel_initializer='he_normal',kernel_regularizer=reg),
        tf.keras.layers.BatchNormalization(),tf.keras.layers.Dropout(dropout),
        tf.keras.layers.Dense(64,activation='relu',kernel_initializer='he_normal',kernel_regularizer=reg),
        tf.keras.layers.BatchNormalization(),tf.keras.layers.Dropout(dropout),tf.keras.layers.Dense(1)])
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr),loss='mse')
    return model

def callbacks():
    import tensorflow as tf
    return [tf.keras.callbacks.EarlyStopping(patience=30,restore_best_weights=True),
            tf.keras.callbacks.ReduceLROnPlateau(factor=.1,patience=10)]
