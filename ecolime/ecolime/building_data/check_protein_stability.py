import pandas as pd
import argparse


# input 
parser = argparse.ArgumentParser(description='Check Folding Parameters')
parser.add_argument('temperature', type=int, nargs=1, help='temperature')
parser.add_argument('gene',type=str, nargs=1,help='gene')
args = parser.parse_args()

temperature = args.temperature
T = temperature[0]
gene = args.gene
gene = gene[0]

df=pd.read_csv('Folding_Rates_matrix_slope22000.csv')
cols = ['protein_length','k_f_'+str(T)+'C']
df = df[df.genome_region==gene]
print (df[cols])


df=pd.read_csv('Oobatake_Keq_matrix.csv')
cols=['Oobatake_Keq_'+str(T)+'C']
df = df[df.genome_region==gene]
print (df[cols])

df=pd.read_csv('Dill_dG_matrix.csv')
cols=['Dill_Keq_'+str(T)+'C']
df = df[df.genome_region==gene]
print (df[cols])
