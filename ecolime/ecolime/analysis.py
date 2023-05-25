import pickle
import pandas as pd
import os.path
import re
import glob
import argparse
import math

# load general data type
me = pickle.load(open('Models/FoldME_py36_ps_1.pickle','rb'))

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
MW_dict['dummy']=33710.96300000001 # average molecular weight for E. coli proteins
mass_glc=0.1800634
parent_complex_dict=pickle.load(open('Models/parent_complex_dict.pickle','rb'))
complex_stoi_dict=pickle.load(open('Models/complex_stoi_dict.pickle','rb'))

# Calculation procedures
def get_complex_formation_flux(sol):
    complex_formation_flux_dict={}
    for rxn,flux in sol.x_dict.items():
        if rxn.startswith('formation_') and abs(flux)>1e-20:
            cplx = rxn.replace('formation_','')
            count=1
            if cplx in complex_stoi_dict.keys():
                if cplx in parent_complex_dict.keys():
                    parents=parent_complex_dict[cplx]
                    for pr in parents:
                        try:
                            if abs(sol.x_dict['formation_'+pr])>1e-20:
                                count=0
                        except:
                            continue
            if count:
                complex_formation_flux_dict[cplx]=flux

    return complex_formation_flux_dict

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

def get_met_stoichiometry(rxn,stoi,mu):
    try:
        float(stoi)
        return stoi
    except:
        try:
            a=stoi.args[0]
            b=float(str(stoi.args[1]).replace('*mu',''))
            stoi_new = a+mu*b
        except:
            a = stoi.args[0]
            b = float(str(stoi.args[1]).replace('*(mu + 0.3915)/mu',''))
            stoi_new = a + b*(mu+0.3915)/mu
        return stoi_new

def get_metabolite_balance(sol,met_id):
    mu=sol.x_dict['biomass_dilution']
    met=me.metabolites.get_by_id(met_id)
    met_production_rxn_dict={}
    met_consumption_rxn_dict={}
    production=0
    consumption=0
    for rxn in met.reactions:
        ID=rxn.id
        stoi=rxn.get_coefficient(met_id)
        stoi=get_met_stoichiometry(rxn,stoi,mu)
        
        flux=sol.x_dict[ID]
        if abs(flux) > 1e-20: # precision limit
            if stoi*flux>0:
                met_production_rxn_dict[ID]={'flux':flux,'stoichiometry':stoi}
                production+=stoi*flux
            elif stoi*flux<0:
                met_consumption_rxn_dict[ID]={'flux':flux,'stoichiometry':stoi}
                consumption+=stoi*flux

    return met_consumption_rxn_dict,met_production_rxn_dict


def check_met_production_reactions(solution_file_path,met_id):

    met_production_rxns = []
    met_consumption_rxns = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        T = int(f.split('_')[1])
        sample = int(f.split('_')[3])
        if os.path.isfile(filename):

            # load solution
            sol = pickle.load(open(filename,'rb'))
            
            mu = sol.x_dict['biomass_dilution']
            met_consumption_rxn_dict,met_production_rxn_dict = get_metabolite_balance(sol,met_id)

            met_production_rxns += met_production_rxn_dict.keys()
            met_consumption_rxns += met_consumption_rxn_dict.keys()

    met_production_rxns = list(set(met_production_rxns))
    met_consumption_rxns = list(set(met_consumption_rxns))
    if met_id == 'atp_c':
        # consider metabolic reactions, others may appear because of precision/numerical issue in solution  
        met_production_rxns = [r for r in met_production_rxns if 'folding_' not in r and 'transcription_' not in r and 'translation' not in r]
        met_production_rxns = list(set(met_production_rxns))
   
    print("")
    print(met_id,"producing reactions:")
    for r in met_production_rxns:
        print (r)

    if met_id != 'atp_c':
        print ("")
        print (met_id, "consuming reactions:")
        for r in met_consumption_rxns:
            print (r)

def get_net_flux(sol,rxnid):
    flux_a = sol.x_dict[rxnid]
    if '_REV_' in rxnid:
        rxnid_b = rxnid.replace('_REV_','_FWD_')
        flux_b = sol.x_dict[rxnid_b]
    elif '_FWD_' in rxnid and not rxnid.startswith('PYK_'):
        rxnid_b = rxnid.replace('_FWD_','_REV_')    
        flux_b = sol.x_dict[rxnid_b]
    else:
        flux_b = 0.0
    result=flux_a-flux_b
    if abs(result) < 1e-20:
        result=0.0
    return result

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

def if_production(flux):
    if flux > 0:
        return flux
    else:
        return 0.0

# Output procedures
def write_phenotype(solution_file_path,output_filename):
     
    """ write out phenotypic data from sampling simulations.

        Parameters
        ----------
        solution_file_path : path to where the solution files are located.
                             All solution file should be in named the format 'T_'+str(T)+'_sample_'+str(sample)+'_sol.pickle'
                             where 'T' is the simulated temperature, 'sample' is the sample ID
        output_filename : output file name
    """
    phenotype_matrix = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
      f = os.path.basename(filename)
      f = f.replace('_sol.pickle','')
      sample = f[-1]
      strain = f[0:-2]
      output_dict = {'strain':f}
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
          lac = -1 * sol.x_dict['EX_lac__D_e']
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
         
          peptide_per_complex = []
          cofactor_per_complex = []
          complex_formation_flux_dict = get_complex_formation_flux(sol) 
          for cplx,flux in complex_formation_flux_dict.items():
              genes = complex_stoi_dict[cplx]
              L = 0
              cofactor = 0
              for gene in genes:
                  if gene.startswith('b'):
                      try:
                          L += complex_stoi_dict[cplx][gene]
                          expressed_genes.append(gene)
                      except:
                          cofactor += float(complex_stoi_dict[cplx][gene])
                          continue
                  else:
                      cofactor += float(complex_stoi_dict[cplx][gene])
              peptide_per_complex.append(float(L))
              cofactor_per_complex.append(cofactor)
          output_dict['peptide_per_complex'] = sum(peptide_per_complex)/len(peptide_per_complex)
          output_dict['cofactor_per_complex'] = sum(cofactor_per_complex)/len(cofactor_per_complex)
          phenotype_matrix.append(output_dict)
    
    phenotype_matrix = pd.DataFrame(phenotype_matrix)
    cols = ['strain','growth_rate','GUR','OUR','APR','LPR','biomass_yield','Nexpressed','peptide_per_complex','cofactor_per_complex']
    phenotype_matrix = phenotype_matrix[cols]
    #phenotype_matrix = phenotype_matrix.sort_values(by=['sample'],ascending=True)
    #phenotype_matrix = phenotype_matrix.reset_index(drop=True)
    #phenotype_matrix = phenotype_matrix.set_index('simID')
    phenotype_matrix = phenotype_matrix.sort_values(['strain'])
    phenotype_matrix.to_csv(output_filename)


def write_ATP_production(solution_file_path,output_filename):

    energy_matrix = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        strain = f[0:-2]
        output_dict = {'strain':f}

        if os.path.isfile(filename):

            # load solution
            sol = pickle.load(open(filename,'rb'))
            sol.x_dict['DM_RNA_b0439'] = 0
 
            # get metabolic reaction fluxes
            met_flux = me.get_metabolic_flux(sol)
            ID_Algn(met_flux)
            
            # get ATP-production fluxes
            output_dict['ATPS4rpp'] = if_production(met_flux['ATPS4rpp'])
            output_dict['PGK'] = if_production(-1.* met_flux['PGK'])
            output_dict['ACKr'] = if_production(-1.* met_flux['ACKr'])
            output_dict['PYK'] = if_production(met_flux['PYK'])
            #output_dict['PPKr'] = if_production(-1.* met_flux['PPKr'])
            output_dict['PPK2r'] = if_production(-1.* met_flux['PPK2r'])
            #output_dict['SUCOAS'] = if_production(-1.* met_flux['SUCOAS'])
            #output_dict['PRPPS'] = if_production(-1.* met_flux['PRPPS'])
            output_dict['GLNS'] = if_production(met_flux['GLNS'])
            energy_matrix.append(output_dict)

    energy_matrix = pd.DataFrame(energy_matrix)
    #cols = ['strain','sample','ATPS4rpp','PGK','ACKr','PYK','PPKr','PPK2r','SUCOAS','PRPPS']
    cols = ['strain','ATPS4rpp','PGK','ACKr','PYK','PPK2r','GLNS']
    energy_matrix = energy_matrix[cols]
    #energy_matrix = energy_matrix.sort_values(by=['T','sample'],ascending=True)
    #energy_matrix = energy_matrix.reset_index(drop=True)
    #energy_matrix = energy_matrix.set_index('simID')
    energy_matrix.to_csv(output_filename)
            

def write_protein_abundance(solution_file_path,output_filename):

    mf_matrix = []
    allExpressedGenes = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        strain = f[0:-2]
        output_dict = {'strain':f}

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
    cols = ['strain'] + allExpressedGenes
    mf_matrix = mf_matrix[cols]
    mf_matrix = mf_matrix.sort_values(['strain'])
    #mf_matrix = mf_matrix.set_index('simID')
    print ("")
    print ('Number of expressed genes = ',len(allExpressedGenes))
    mf_matrix.to_csv(output_filename)
    
def write_metabolic_flux(solution_file_path,output_filename):

    metabolic_flux_matrix = []
    metabolic_rxns = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        strain = f[0:-2]

        if os.path.isfile(filename):
            # load solution
            sol = pickle.load(open(filename,'rb'))
            sol.x_dict['DM_RNA_b0439'] = 0
            # calculate flux
            met_flux = me.get_metabolic_flux(sol)
            ID_Algn(met_flux)
            output_dict = find_nonzero_flux(met_flux)
            output_dict['strain']=f
            output_dict['sample']=sample
            metabolic_rxns += output_dict.keys()
            metabolic_rxns = list(set(metabolic_rxns))
            metabolic_flux_matrix.append(output_dict)

    metabolic_rxns = sorted(metabolic_rxns)
    print ('Number of non-zero metabolic fluxes = ',len(metabolic_rxns))
    metabolic_flux_matrix = pd.DataFrame(metabolic_flux_matrix)
    metabolic_flux_matrix = metabolic_flux_matrix.fillna('0.0')
    cols = ['strain']+metabolic_rxns
    metabolic_flux_matrix = metabolic_flux_matrix[cols]
    #metabolic_flux_matrix = metabolic_flux_matrix.sort_values(by=['T','sample'],ascending=True)
    #metabolic_flux_matrix = metabolic_flux_matrix.reset_index(drop=True)
    #metabolic_flux_matrix = metabolic_flux_matrix.set_index('simID')
    metabolic_flux_matrix.to_csv(output_filename)
  

def write_q8_flux(solution_file_path,output_filename):

    '''

    First define q8 producing/consuming reactions:
    -- The list of q8 production and consumption reactions are generated by running check_met_production_reactions(solution_file_path,'q8_c')
    -- We are interested in the usage of different ETC enzymes to achieve a certain phenotype, 
       hence unlike metabolic reactions, here the same reaction catalyzed by different enzymes are separated into multiple entries. 
    -- We list out all non-zero q8 production/consumption reactions in our sampling simulations as a reference.

    '''

    #q8_production_reactions=['CYTBDpp_FWD_APP-UBIOX-CPLX_mod_pheme_mod_hemed',\
    #                         'CYTBDpp_FWD_CYT-D-UBIOX-CPLX_mod_pheme_mod_hemed',\
    #                         'CYTBO3_4pp_FWD_CYT-O-UBIOX-CPLX_mod_pheme_mod_hemeO_mod_cu2']
    #q8_consumption_reactions=['NADH5_FWD_NADH-DHII-MONOMER_mod_mg2_mod_cu_mod_fad',\
    #                          'NADH16pp_FWD_NADH-DHI-CPLX_mod_2fe2s_mod_4fe4s_mod_fmn',\
    #                          'FDH4pp_FWD_FORMATEDEHYDROGN-CPLX_mod_bmocogdp',\
    #                          'FDH4pp_FWD_FORMATEDEHYDROGO-CPLX_mod_bmocogdp',\
    #                          'SUCDi_FWD_SUCC-DEHASE_mod_3fe4s_mod_fad_mod_2fe2s_mod_4fe4s',\
    #                          'LDH_D2_FWD_DLACTDEHYDROGFAD-MONOMER_mod_fad',\
    #                          'NADPHQR2_FWD_CPLX0-253_mod_fad',\
    #                          'POX_FWD_PYRUVOXID-CPLX_mod_mg2_mod_thmpp_mod_fad',\
    #                          'GLCDpp_FWD_GLUCDEHYDROG-MONOMER_mod_pqq',\
    #                          'GLCDpp_FWD_G6437-MONOMER_mod_ca2_mod_pqq',\
    #                          'MDH2_FWD_EG12069-MONOMER_mod_fad',\
    #                          'GLYCTO2_FWD_CPLX0-7458',\
    #                          'ASPO3_FWD_L-ASPARTATE-OXID-MONOMER_mod_fad',\
    #                          'DHORD2_FWD_DIHYDROOROTOX-MONOMER_mod_fmn']
    q8_production_reactions=['NTRIR3pp_FWD_NRFMULTI-CPLX',\
                            'NO3R1pp_FWD_NITRATREDUCTZ-CPLX_mod_bmocogdp_mod_2:pheme_mod_3:4fe4s_mod_1:3fe4s',\
                            'NO3R1pp_FWD_NITRATREDUCTA-CPLX_mod_bmocogdp_mod_3fe4s_mod_4fe4s'] 
    q8_consumption_reactions=['ASPO3_FWD_L-ASPARTATE-OXID-MONOMER_mod_fad',\
                              'SUCDi_FWD_SUCC-DEHASE_mod_3fe4s_mod_fad_mod_2fe2s_mod_4fe4s',\
                              'FDH4pp_FWD_FORMATEDEHYDROGO-CPLX_mod_bmocogdp',\
                              'DHORD2_FWD_DIHYDROOROTOX-MONOMER_mod_fmn',\
                              'NADH16pp_FWD_NADH-DHI-CPLX_mod_2fe2s_mod_4fe4s_mod_fmn',\
                              'FDH4pp_FWD_FORMATEDEHYDROGN-CPLX_mod_bmocogdp',\
                              'MDH2_FWD_EG12069-MONOMER_mod_fad',\
                              'GLYCTO2_FWD_CPLX0-7458']


    q8_consumption_reaction_dict={}
    q8_consumption_stoi_dict={}
    for rxn in q8_consumption_reactions+q8_production_reactions:
        if rxn=='FDH4pp_FWD_FORMATEDEHYDROGN-CPLX_mod_bmocogdp':
            rxn_short = 'FDH4pp_fdn'
        elif rxn=='FDH4pp_FWD_FORMATEDEHYDROGO-CPLX_mod_bmocogdp':
            rxn_short = 'FDH4pp_fdo'
        elif rxn=='CYTBDpp_FWD_APP-UBIOX-CPLX_mod_pheme_mod_hemed':
            rxn_short = 'CYTBDpp_APP'
        elif rxn=='CYTBDpp_FWD_CYT-D-UBIOX-CPLX_mod_pheme_mod_hemed':
            rxn_short = 'CYTBDpp_CYD'
        elif rxn=='GLCDpp_FWD_G6437-MONOMER_mod_ca2_mod_pqq':
            rxn_short = 'GLCDpp_ylil'
        elif rxn=='GLCDpp_FWD_GLUCDEHYDROG-MONOMER_mod_pqq':
            rxn_short = 'GLCDpp_gcd'
        elif rxn=='NO3R1pp_FWD_NITRATREDUCTZ-CPLX_mod_bmocogdp_mod_2:pheme_mod_3:4fe4s_mod_1:3fe4s':
            rxn_short = 'NO3R1pp_pheme'
        else:
            rxn_short = rxn.split('_')[0]
        q8_consumption_reaction_dict[rxn]=rxn_short
        reaction=me.reactions.get_by_id(rxn)
        q8_consumption_stoi_dict[rxn_short]=reaction.get_coefficient('q8_c')
    
    # calculation
    q8_flux_matrix = []
    q8_rxns = []
    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        strain = f[0:-2]

        if os.path.isfile(filename):

            output_dict = {'strain':strain,'sample':sample}
            # load solution
            sol = pickle.load(open(filename,'rb'))

            mu = sol.x_dict['biomass_dilution']
            for rxn,rxn_short in q8_consumption_reaction_dict.items():
                flux = sol.x_dict[rxn] * q8_consumption_stoi_dict[rxn_short]
                output_dict[rxn_short]=flux
            q8_flux_matrix.append(output_dict)

    #cols=['simID','T','sample','NADH5','NADH16pp','FDH4pp_fdn','FDH4pp_fdo','SUCDi',\
    #      'LDH','NADPHQR2','POX','GLCDpp_gcd','GLCDpp_ylil','MDH2','GLYCTO2','ASPO3','DHORD2',\
    #      'CYTBDpp_APP','CYTBDpp_CYD','CYTBO3']
    cols=['strain','sample','NTRIR3pp','NO3R1pp_pheme','NO3R1pp','ASPO3','SUCDi',\
          'NADH16pp','FDH4pp_fdn','FDH4pp_fdo','DHORD2','MDH2','GLYCTO2']
    q8_flux_matrix = pd.DataFrame(q8_flux_matrix)
    q8_flux_matrix = q8_flux_matrix[cols]
    #q8_flux_matrix = q8_flux_matrix.sort_values(by=['T','sample'],ascending=True)
    #q8_flux_matrix = q8_flux_matrix.reset_index(drop=True)
    #q8_flux_matrix = q8_flux_matrix.set_index('simID')
    q8_flux_matrix.to_csv(output_filename)

#def write_mqn8_flux(solution_file_path,output_filename):
#
#    '''
#
#    First define mqn8 producing/consuming reactions:
#    -- The list of mqn8 production and consumption reactions are generated by running check_met_production_reactions(solution_file_path,'mqn8_c')
#    -- We are interested in the usage of different ETC enzymes to achieve a certain phenotype,
#       hence unlike metabolic reactions, here the same reaction catalyzed by different enzymes are separated into multiple entries.
#    -- We list out all non-zero mqn8 production/consumption reactions in our sampling simulations as a reference.
#
#    '''
#
#    mqn8_production_reactions=['NO3R2bpp_FWD_NAPAB-CPLX_NAPC-MONOMER_mod_bmocogdp',\
#                             'NTRIR4pp_FWD_NRFMULTI-CPLX',\
#                             'NO3R2pp_FWD_NITRATREDUCTA-CPLX_mod_bmocogdp_mod_3fe4s_mod_4fe4s',\
#                             'NO3R2pp_FWD_NITRATREDUCTZ-CPLX_mod_bmocogdp_mod_2:pheme_mod_3:4fe4s_mod_1:3fe4s']
#    mqn8_consumption_reactions=['MDH3_FWD_EG12069-MONOMER_mod_fad',\
#                                'ASPO4_FWD_L-ASPARTATE-OXID-MONOMER_mod_fad',\
#                                'GLYCTO3_FWD_CPLX0-7458',\
#                                'FDH5pp_FWD_FORMATEDEHYDROGO-CPLX_mod_bmocogdp',\
#                                'NADH17pp_FWD_NADH-DHI-CPLX_mod_2fe2s_mod_4fe4s_mod_fmn',\
#                                'FDH5pp_FWD_FORMATEDEHYDROGN-CPLX_mod_bmocogdp']

def write_mqn8_flux(solution_file_path,output_filename):

    mqn8_rxn_dict=pickle.load(open('mqn8_rxn_dict.pickle','rb'))

    # calculation
    mqn8_flux_matrix = []

    for filename in glob.glob(solution_file_path+"/*sol.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        strain = f[0:-2]
        output_dict = {'strain':strain,'sample':sample}

        if os.path.isfile(filename):

            # load solution
            sol = pickle.load(open(filename,'rb'))
            mu = sol.x_dict['biomass_dilution']
            for rxn_short in mqn8_rxn_dict.keys():
                rxn = mqn8_rxn_dict[rxn_short]['rxn_id']
                flux = sol.x_dict[rxn] * mqn8_rxn_dict[rxn_short]['coeff_mqn8_c']
                if abs(flux)<1e-20:
                    flux = 0.0
                output_dict[rxn_short]=flux
            mqn8_flux_matrix.append(output_dict)

    cols=['strain','sample']
    for col in mqn8_rxn_dict.keys():
        cols.append(col)
    mqn8_flux_matrix=pd.DataFrame(mqn8_flux_matrix)

    # check if column is zero
    cols=['strain','sample']
    for col in mqn8_rxn_dict.keys():
        if (mqn8_flux_matrix[col] == 0).all():
            print (col,)
        else:
           cols.append(col)

    mqn8_flux_matrix=mqn8_flux_matrix[cols]
    #mqn8_flux_matrix = mqn8_flux_matrix.sort_values(by=['T','sample'],ascending=True)
    #mqn8_flux_matrix = mqn8_flux_matrix.reset_index(drop=True)
    #mqn8_flux_matrix = mqn8_flux_matrix.set_index('simID')
    mqn8_flux_matrix.to_csv(output_filename)

def write_mutation(solution_file_path,output_prefix):
    keff_matrix = []
    ddG_matrix = []
    for filename in glob.glob(solution_file_path+"/*_perturbation_dict.pickle"):
        f = os.path.basename(filename)
        f = f.replace('_sol.pickle','')
        sample = f[-1]
        strain = f[0:-2]
        output_dict_keff = {'strain':strain,'sample':sample}
        output_dict_ddG = {'strain':strain, 'sample':sample}
        if os.path.isfile(filename):

            # load solution
            mutation_dict = pickle.load(open(filename,'rb'))
            # go through mutations
            for gene in model_genes:
                if gene in mutation_dict.keys():
                    mutation = mutation_dict[gene]
                    output_dict_keff[gene] = mutation['keff_factor']
                    output_dict_ddG[gene] = mutation['ddG']
                else:
                    output_dict_keff[gene] = 1.0
                    output_dict_ddG[gene] = 0.0
            keff_matrix.append(output_dict_keff)
            ddG_matrix.append(output_dict_ddG)

    cols=['strain','sample'] + model_genes
    keff_matrix=pd.DataFrame(keff_matrix)
    keff_matrix
    keff_matrix=keff_matrix[cols]
    #keff_matrix = keff_matrix.sort_values(by=['T','sample'],ascending=True)
    #keff_matrix = keff_matrix.reset_index(drop=True)
    #keff_matrix = keff_matrix.set_index('simID')
    keff_matrix.to_csv(output_prefix+'_keff.csv')

    ddG_matrix=pd.DataFrame(ddG_matrix)
    ddG_matrix=ddG_matrix[cols]
    #ddG_matrix = ddG_matrix.sort_values(by=['T','sample'],ascending=True)
    #ddG_matrix = ddG_matrix.reset_index(drop=True)
    #ddG_matrix = ddG_matrix.set_index('simID')
    ddG_matrix.to_csv(output_prefix+'_ddG.csv')


def main(solution_file_path,output_prefix,check_ATP):
     
    #write_mqn8_flux(solution_file_path,output_prefix+'mqn8_flux.csv') 
    write_phenotype(solution_file_path,output_prefix+'_phenotype.csv')
    write_protein_abundance(solution_file_path,output_prefix+'_mf.csv')
    #write_metabolic_flux(solution_file_path,output_prefix+'_metabolic_flux.csv')
    #write_mutation(solution_file_path,output_prefix)
 
    #if check_ATP:
    #    print ""
    #    print "Check if additional ATP-production metabolic reactions are activated,"
    #    print "except: ATPS4rpp,PGK,ACKr,PYK,PPKr,PPK2r,SUCOAS,PRPPS."
    #    print "After confirming, modify write_ATP_production() with additional reactions, and rerun."
    #    check_met_production_reactions(solution_file_path,'atp_c') 
    #else:
    #    write_ATP_production(solution_file_path,output_prefix+'_ATP.csv')
    #
    #if check_q8:
    #    print ""
    #    print "Check if additional q8-producing/consuming metabolic reactions are activated."
    #    print "After confirming, modify write_q8_production() with additional reactions, and rerun."
    #    check_met_production_reactions(solution_file_path,'q8_c')
    #else:
    #    write_q8_flux(solution_file_path,output_prefix+'_q8.csv')
    #
    #if check_mqn8:
    #    print ""
    #    print "Check if additional mqn8-producing/consuming metabolic reactions are activated."
    #    print "After confirming, modify write_mqn8_production() with additional reactions, and rerun."
    #    check_met_production_reactions(solution_file_path,'mqn8_c')
    
if __name__=="__main__":
   
    '''
    Example: python analysis.py 'TEST_Solutions/' 'Analysis_Output/propensity_scaling' 0
    '''
    
    parser = argparse.ArgumentParser(description='Process ME-model solution file, and write out data')
    parser.add_argument('path', type=str, nargs=1, help='solution_file_path')
    parser.add_argument('prefix', type=str, nargs=1, help='output_file_prefix')
    parser.add_argument('atp', type=int, nargs=1, help='if necessary to recheck ATP-production reactions')
    #parser.add_argument('q8', type=int, nargs=1, help='if necessary to recheck q8-producing/consuming reactions')
    #parser.add_argument('mqn8', type=int, nargs=1, help='if necessary to recheck mqn8-producing/consuming reactions')
    args = parser.parse_args()

    path = args.path
    solution_file_path = path[0]
    prefix = args.prefix
    output_prefix = prefix[0]
    atp = args.atp
    check_ATP = atp[0]
    #q8 = args.q8
    #check_q8 = q8[0]
    #mqn8 = args.mqn8
    #check_mqn8 = mqn8[0]

    main(solution_file_path,output_prefix,check_ATP)
