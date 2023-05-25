## Example ##
## python simulation.py 37##

import cobrame
import ecolime
from cobra import Reaction
from cobra import Metabolite
from qminospy.me1 import ME_NLP1

import numpy as np
from numpy.random import normal, uniform
import scipy
import math
import pickle
import pandas as pd
import argparse


# input 
parser = argparse.ArgumentParser(description='Run temperature simulations')
parser.add_argument('temperature', type=int, nargs=1, help='temperature')
parser.add_argument('ML_keff', type=int, nargs=1, help='1:use David keff; 0: use SASA-scaled keff')
parser.add_argument('propensity_scaling', type=float, nargs=1, help='global parameter: propensity scaling')
args = parser.parse_args()

temperature = args.temperature
temperature = temperature[0]
ML_keff = args.ML_keff
ML_keff = ML_keff[0]
ps = args.propensity_scaling
ps = ps[0]

MU_PREC = 1e-3
MU_MIN  = 0.1
MU_MAX  = 1.1

# load FoldME model
modelfile = 'Models/FoldME_py36_ps_'+str(ps)+'_pd_0.50_pg_0.05_ETCfix.pickle'
with open(modelfile,'rb') as iofile:
    me = pickle.load(iofile)

int_temp = int(temperature)
str_temp = str(temperature)
me.global_info['temperature'] = int_temp

# use David's ML keffs, should be used before change temperature 
if ML_keff:
    ecolime.chaperones.update_ML_keffs(me)
# change temperature
ecolime.chaperones.change_temperature(me,temperature)

# anaerobic 
#me.reactions.EX_o2_e.upper_bound = 0.0
#me.reactions.EX_o2_e.lower_bound = 0.0

# ECOM4
#me.reactions.get_by_id('translation_b3029').upper_bound = 0 # ygiN
#me.reactions.get_by_id('translation_b3029').lower_bound = 0

# solve model
me_nlp = ME_NLP1(me, growth_key='mu')
muopt, hs, xopt, cache = me_nlp.bisectmu(precision=MU_PREC, mumin=MU_MIN, mumax=MU_MAX)

# write solution file 
pd_pg_str='_pd_050_pg_005_'
if ML_keff == 1:
    filename = 'T_'+str(temperature)+'_ps_'+str(ps)+pd_pg_str+'ETCfixi_MLkeff' 
elif ML_keff == 2:
    filename = 'T_'+str(temperature)+'_ps_'+str(ps)+pd_pg_str+'ETCfixi_AliKeff' 
else:
    filename = 'T_'+str(temperature)+'_ps_'+str(ps)+pd_pg_str+'ETCfix'

with open('TEST_Solutions/'+filename+'_sol.pickle','wb') as f:
    pickle.dump(me.solution,f)
#
sol=me.solution
print("Growth rate at temperature %d = %f" %(temperature,sol.x_dict['biomass_dilution']))
