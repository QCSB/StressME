###########################################
# StressME for heat_evolved e.coli
###########################################

import cobrame
import ecolime
import AcidifyME
import pickle
from cobrame.core.reaction import TranscriptionReaction
from cobrame.io.json import load_json_me_model
from cobra import Reaction
from cobra import Metabolite
from AcidifyME import *
from AcidifyME.periplasmic_proteome import *
from AcidifyME.membrane_lipid_constraint import *
from AcidifyME.membrane_protein_activity import *
from AcidifyME.proton_influx import *
from oxidizeme.model import *
from ecolime import chaperones
from oxidizeme.model import *

import argparse
import os
import os.path
from os.path import abspath, dirname
import numpy as np
from numpy.random import normal, uniform
import scipy
import math
import pandas as pd
import re
import re
import glob
import math

from qminospy.me1 import ME_NLP1

here = dirname(abspath(__file__))
############################
# input setup
############################
parser = argparse.ArgumentParser(description='Run temperature simulations')
parser.add_argument('temperature', type=int, nargs=1, help='temperature')
parser.add_argument('pH', type=float, nargs=1, help='pH')
parser.add_argument('ROS', type=float, nargs=1, help='ROS')
args = parser.parse_args()

temperature = args.temperature
temperature = temperature[0]

pH = args.pH
pH = pH[0]

ROS = args.ROS
ROS = ROS[0]

#################################
#load the FoldME using wild type Keffs
#################################
with open("FoldME_Ali_keff.pickle", "rb") as f:
  me = pickle.load(f)
#################################
#To change a key Keff for DXPRIi so that the strain can be heat-evolved
#################################
rxn = me.reactions.get_by_id('DXPRIi_FWD_DXPREDISOM-CPLX_mod_cobalt2')
rxn.keff=88.72
rxn.update()
#################################
# To set up thermal stress
#################################
ecolime.chaperones.change_temperature(me,temperature)
#################################
# to set up acid stress
#################################
# fold change of ATP synthesis rate at different external pH values
# Regression based on Figure 4A https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007525

fold_change=5.666666666669649*pH**4-1.396666666667337e+02*pH**3+1.286583333333885e+03*pH**2-5.254083333335287e+03*pH**1+8.037000000002481e+03
print(fold_change)

#add different constraints into the ME-model to simulate acid stress
temperatureK=273.15+temperature
add_membrane_constraint(me, 'MJR_habituated')
add_protein_stability_constraint(me, pH, temperatureK, constrain_periplasm_biomass=False)
add_proton_leak_rxn(me, pH)
modify_membrane_protein_activity(me, fold_change)
#################################
# unmodeled protein fraction
#################################
me.unmodeled_protein_fraction = 0.1
#################################
# to set up oxidative stress
#################################
#print("Making and saving stressME")
print("Making stressME")

solver = ME_NLP1(me)

stress = StressME(solver)

stress.make_stressme()

subs_dict = {
 'h2o2_c': 5e-08, 'o2s_c': 2e-10,
 'h2o2_e': 1.5e-06, 'o2s_e': 0,
 'pq2_c': 0.0}
 
subs_dict['h2o2_c'] = 5e-08*ROS
subs_dict['o2s_c'] = 2e-10*ROS
 
stress.substitute_ros(solver, subs_dict)
stress.substitute_metal(solver)

#### Unlock ROS detox enzymes, as they might be blocked (UB=0) for some reason by default
for r in me.process_data.SPODM.parent_reactions:
    print(r.id, '\t', r.upper_bound)
    r.upper_bound = 1000
    print(r.id, '\t', r.upper_bound)
    
for r in me.process_data.SPODMpp.parent_reactions:
    print(r.id, '\t', r.upper_bound)
    r.upper_bound = 1000
    print(r.id, '\t', r.upper_bound)
    
for r in me.process_data.CAT.parent_reactions:
    print(r.id, '\t', r.upper_bound)
    r.upper_bound = 1000
    print(r.id, '\t', r.upper_bound)
    
me.reactions.MOX_REV_CPLX_dummy.upper_bound = 0

stress.force_o2s_from_pq()

# save_json_me_model(me, StressME)

for r in me.reactions.query('PQ2RED'):
    print(r.id, '\t', r.lower_bound)
    
#################################
# test scenario
#################################
# DXPRIi keff is stable between 40 °C and 50 °C
# based on "The maximum activity was observed at 40–60 °C ... 
# The purified enzyme was not heat stable above 50 °C ..."
# J Biol Chem. 2000;275(26):19928-32
# rxn = me.reactions.get_by_id('DXPRIi_FWD_DXPREDISOM-CPLX_mod_cobalt2')
# rxn.keff=88.72
# rxn.update()  
    
#################################
# warm-start from saved basis
#################################
saved_basis = np.load('basis.npy')
MU_PREC = 1e-04
sol, hs, xopt, cache  = solver.bisectmu(precision=MU_PREC, basis = saved_basis)    
#################################
# write solution file 
#################################
filename = 'StressME_heatevolved'+'_T_'+str(temperature)+'_pH_'+str(pH)+'_ROS_'+str(ROS)+'X'
with open(filename+'_sol.pickle','wb') as f:
    pickle.dump(me.solution,f)
#
sol=me.solution

#################################
# Calculation for phenotypes; 
# proteome and fluxome
#################################
def find_nonzero_flux(met_flux,log=0):
    met_flux_nonzero={}
    for rxn,flux in met_flux.items():
        if abs(flux) > 1e-20 and 'DM_RNA' not in rxn:
            if log:
                met_flux_nonzero[rxn]=math.log(flux)
            else:
                met_flux_nonzero[rxn]=flux
    return met_flux_nonzero

def find_overlap(list1,list2):
    return [i for i in list1 if i in list2]

def ID_Algn(met_flux):
   met_flux['L_LACtex']=met_flux['L__LACtex']
   met_flux['L_LACt2rpp']=met_flux['L__LACt2rpp']
   met_flux['L_LACD3']=met_flux['L__LACD3']
   met_flux['L_LACD2']=met_flux['L__LACD2']

   met_flux['D_LACtex']=met_flux['D__LACtex']
   met_flux['D_LACt2pp']=met_flux['D__LACt2pp']

   met_flux['DSBAO2']=met_flux['DSBAO21']

MW_dict={}
protein_length={}
for gene in me.metabolites.query(re.compile('RNA_b[0-9]')):
    prot_id=gene.id.replace('RNA','protein')
    try:
        met = me.metabolites.get_by_id(prot_id)
        MW_dict[prot_id.replace('protein_','')] = met.formula_weight
        protein_length[prot_id.replace('protein_','')] = len(met.amino_acid_sequence)
    except:
        continue # this is a RNA gene
model_genes=sorted(MW_dict.keys())
MW_dict['dummy']=33210.96300000001 # average molecular weight for E. coli proteins
mass_glc=0.1800634

# Calculation procedures
def get_total_protein(sol):
    total_protein = 0
    mf_dict={}
    if type(sol) != type(None):
        for a,b in sol.x_dict.items():
            if a.startswith('translation_') and b > 1e-20:
                gene = a.replace('translation_','')
                flux = sol.x_dict['translation_'+gene]
                mw=MW_dict[gene]
                total_protein += flux*mw
                mf_dict[gene] = flux*mw
        return total_protein,mf_dict
    else:
        print("No feasible solution!")
        return -1

# Output procedures

def write_protein_abundance(solution_file_path,output_filename):
    import glob
    mf_matrix = []
    allExpressedGenes = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        condition = f[0:-2]
        output_dict = {'condition':f}

        if os.path.isfile(filename):

            # load solution
            sol = pickle.load(open(filename,'rb')) 
            
            # calculate mass fraction 
            total_protein,mf_dict = get_total_protein(sol)
            for gene in mf_dict.keys():
                output_dict[gene] = mf_dict[gene]/total_protein
                allExpressedGenes.append(gene)
            allExpressedGenes=list(set(allExpressedGenes))
            mf_matrix.append(output_dict)

    allExpressedGenes=sorted(allExpressedGenes)
    mf_matrix = pd.DataFrame(mf_matrix)
    mf_matrix = mf_matrix.fillna('0.0')
    #mf_matrix = mf_matrix.sort_values(by=['T','sample'],ascending=True)
    #mf_matrix = mf_matrix.reset_index(drop=True)
    cols = ['condition'] + allExpressedGenes
    mf_matrix = mf_matrix[cols]
    mf_matrix = mf_matrix.sort_values(['condition'])
    #mf_matrix = mf_matrix.set_index('simID')
    print ("")
    print ('Number of expressed genes = ',len(allExpressedGenes))
    mf_matrix.to_csv(output_filename)
      
solution_file_path='/home/meuser'
output_filename='StressME_heatevolved_proteome.csv'
write_protein_abundance(solution_file_path,output_filename)

def write_phenotype(solution_file_path,output_filename):
     
    """ write out phenotypic data from sampling simulations.

        Parameters
        ----------
        solution_file_path : path to where the solution files are located.
                             All solution file should be in named the format 'T_'+str(T)+'_sample_'+str(sample)+'_sol.pickle'
                             where 'T' is the simulated temperature, 'sample' is the sample ID
        output_filename : output file name
    """
    import glob
    phenotype_matrix = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
      f = os.path.basename(filename)
      f = f.replace('_sol.pickle','')
      sample = f[-1]
      condition = f[0:-2]
      output_dict = {'condition':f}
      if os.path.isfile(filename):
          # load solution
          sol = pickle.load(open(filename,'rb'))
          
          # growth and exchange rate
          mu = sol.x_dict['biomass_dilution']
          gur = -1*sol.x_dict['EX_glc__D_e']
          if abs(gur) < 1e-20:
              gur = 0.0
              biomass_yield = 0.0
          else:
              biomass_yield = mu/gur/mass_glc
          apr = sol.x_dict['EX_ac_e']
          if abs(apr) < 1e-20:
              apr = 0.0
          our = -1 * sol.x_dict['EX_o2_e']
          if abs(our) < 1e-20:
              our = 0.0
          lac = sol.x_dict['EX_lac__D_e']
          if abs(lac) < 1e-20:
              lac = 0.0
          output_dict['growth_rate'] = mu
          output_dict['GUR'] = gur
          output_dict['biomass_yield'] = biomass_yield
          output_dict['APR'] = apr
          output_dict['OUR'] = our
          output_dict['LPR'] = lac

          # proteome complexity
          total_protein,mf_dict = get_total_protein(sol)
          output_dict['Nexpressed'] = len(mf_dict)
         
 
          phenotype_matrix.append(output_dict)
    
    phenotype_matrix = pd.DataFrame(phenotype_matrix)
    cols = ['condition','growth_rate','GUR','OUR','APR','LPR','biomass_yield','Nexpressed']
    phenotype_matrix = phenotype_matrix[cols]
    #phenotype_matrix = phenotype_matrix.sort_values(by=['sample'],ascending=True)
    #phenotype_matrix = phenotype_matrix.reset_index(drop=True)
    #phenotype_matrix = phenotype_matrix.set_index('simID')
    phenotype_matrix = phenotype_matrix.sort_values(['condition'])
    phenotype_matrix.to_csv(output_filename)

solution_file_path='/home/meuser'
output_filename='StressME_heatevolved_phenotypes.csv'
write_phenotype(solution_file_path,output_filename)

def write_metabolic_flux(solution_file_path,output_filename):
    import glob
    metabolic_flux_matrix = []
    metabolic_rxns = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        condition = f[0:-2]

        if os.path.isfile(filename):
            # load solution
            sol = pickle.load(open(filename,'rb'))
            sol.x_dict['DM_RNA_b0439'] = 0
            # calculate flux
            try:
                met_flux = me.get_metabolic_flux(sol)
                ID_Algn(met_flux)
                output_dict = find_nonzero_flux(met_flux)
                output_dict['condition']=f
                output_dict['sample']=sample
                metabolic_rxns += output_dict.keys()
                metabolic_rxns = list(set(metabolic_rxns))
                metabolic_flux_matrix.append(output_dict)
            except:
                print(filename)

    metabolic_rxns = sorted(metabolic_rxns)
    print ('Number of non-zero metabolic fluxes = ',len(metabolic_rxns))
    metabolic_flux_matrix = pd.DataFrame(metabolic_flux_matrix)
    metabolic_flux_matrix = metabolic_flux_matrix.fillna('0.0')
    cols = ['condition']+metabolic_rxns
    metabolic_flux_matrix = metabolic_flux_matrix[cols]
    #metabolic_flux_matrix = metabolic_flux_matrix.sort_values(by=['T','sample'],ascending=True)
    #metabolic_flux_matrix = metabolic_flux_matrix.reset_index(drop=True)
    #metabolic_flux_matrix = metabolic_flux_matrix.set_index('simID')
    metabolic_flux_matrix.to_csv(output_filename)
    
solution_file_path='/home/meuser'
output_filename='StressME_heatevolved_fluxes.csv'
write_metabolic_flux(solution_file_path,output_filename)

