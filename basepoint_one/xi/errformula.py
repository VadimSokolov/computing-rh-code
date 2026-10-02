import numpy as np, json
g=np.concatenate([np.load(f) for f in ['g_1_400.npy','g_401_1000.npy','g_1001_1700.npy']])
B=json.load(open('below.json'))
rho=np.concatenate([0.5+1j*g,0.5-1j*g])
def E(X,al,th):
    a=rho-al
    t=-X**a/a+0.5*X**(a+1j*th)/(a+1j*th)+0.5*X**(a-1j*th)/(a-1j*th)
    return float(np.real(np.sum(t)))
rows=[]
for al in [1.0,0.9,0.75,0.6]:
    for th in [5.0,14.134725]:
        for r in [x for x in B if x['alpha']==al and abs(x['theta']-th)<1e-6]:
            X=r['X']; pred=E(X,al,th)
            corrected=r['total']-pred
            rows.append(dict(alpha=al,theta=th,X=X,err=r['err'],pred=pred,corr_err=corrected-r['exact']))
            print('alpha %.2f theta %6.3f X %.0e observed err %10.3e predicted (1700 zeros) %10.3e  error after correction %10.3e'%(al,th,X,r['err'],pred,corrected-r['exact']))
json.dump(rows,open('errformula.json','w'))
