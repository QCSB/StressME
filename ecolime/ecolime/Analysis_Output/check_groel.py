import pandas as pd

mf=pd.read_csv('propensity_scaling_mf.csv')
mf=mf.sort_values(['strain'], ascending=[1])
df=mf[['strain','b0014','b4143','b0439']]
#df = df.loc[df['strain'].contains('AliKeff')] 
#print (df)
for i in range(0,len(df)):
    row = df.iloc[i]
    if 'AliKeff' not in row.strain and 'pd_025_pg_005' in row.strain:
    
       print(row['strain'], row['b0014'], row['b4143'], row['b0439'])
#print(mf['b0014'].tolist())

#import cobrame
#import ecolime
#import pickle

#me = pickle.load(open('../Models/FoldME_py36_ps_0.1_ETCfix.pickle','rb'))
#print(me.gam,me.ngam)


#rxn_names=['NADH16pp','NADH17pp','NADH18pp','ATPS4rpp']
#for n in rxn_names:
#    rxn=me.reactions.query(n)
#    for r in rxn:
#        print(r.id,r.reaction)

