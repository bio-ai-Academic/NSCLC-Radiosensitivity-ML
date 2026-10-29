import tensorflow as tf
import shap
from sklearn.decomposition import PCA
from workflow import *

x,y=load(); out=output('pca_kernelshap'); source=ROOT/'analysis_outputs/mlp_SF2_5000'
for fold in range(5):
    state=joblib.load(source/f'fold_{fold}_preprocessing.joblib')
    model=tf.keras.models.load_model(source/f'fold_{fold}.keras',compile=False)
    a=state['scale'].transform(state['select'].transform(x.to_numpy()[state['train']]))
    b=state['scale'].transform(state['select'].transform(x.to_numpy()[state['test']]))
    pca=PCA(n_components=.95,svd_solver='full').fit(a)
    def predict(v): return model.predict(v,verbose=0).ravel()*state['y_sd']+state['y_mean']
    def reconstructed(v): return predict(pca.inverse_transform(pca.transform(v)))
    pd.DataFrame(dict(sample_id=x.index[state['test']],original=predict(b),reconstructed=reconstructed(b))).to_csv(out/f'fold_{fold}_prediction_comparison.csv',index=False)
    explainer=shap.KernelExplainer(reconstructed,a[:min(20,len(a))])
    values=np.asarray(explainer.shap_values(b,nsamples=2048,l1_reg='num_features(100)'))
    pd.DataFrame(values,index=x.index[state['test']],columns=x.columns[state['select'].indices_]).to_csv(out/f'fold_{fold}_gene_kernelshap.csv')
