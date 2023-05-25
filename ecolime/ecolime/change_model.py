## Example ##
## python simulation.py 37##

import cobrame
import ecolime
from cobra import Reaction
from cobra import Metabolite
from qminospy.me1 import ME_NLP1

import scipy
import math
import pickle
import pandas as pd

# load FoldME model
modelfile = 'Models/FoldME_new_sasa_david.pickle'
with open(modelfile,'rb') as iofile:
    me = pickle.load(iofile)

# use David's ML keffs, should be used before change temperature 
#ecolime.chaperones.update_ML_keffs(me)

# #reactions for genes that only expressed in sasa model, not in david model
#rxn_keff_dict={'MPTSS1_FWD_EG10154-MONOMER':36.650302321444784, 'BMOCOS_FWD_EG10153-MONOMER':53.34083665106091, 'MOCOS_FWD_EG10153-MONOMER':53.34083665106091} # b0826
#rxn_keff_dict['CPMPS_FWD_EG11595-MONOMER_EG11666-MONOMER']=62.83982200999657 # b0781,b0783
#rxn_keff_dict['MOADSUx1_FWD_CPLX_dummy']=43.640419650213026 # b0784 moaD
#rxn_keff_dict['GLYCL_FWD_GCVMULTI-CPLX_mod_lipo']=260.0390765173197 # gcvHPT
#rxn_keff_dict['SUCDi_FWD_SUCC-DEHASE_mod_3fe4s_mod_fad_mod_2fe2s_mod_4fe4s']=113.06527038180563 # sdhABCD
#rxn_keff_dict['DHAPT_FWD_CPLX0-2081']=167.04671464514306 # b1198-b1200 dhaMLK
#rxn_keff_dict['GLYCTO2_FWD_CPLX0-7458']=125.08840691836666 # b4467-8 glcDE(F)
#rxn_keff_dict['GLYCTO3_FWD_CPLX0-7458']=125.08840691836666 # b4467-8 glcDE(F)
#rxn_keff_dict['GLYCTO4_FWD_CPLX0-7458']=125.08840691836666 # b4467-8 glcDE(F)

#for rxn_name,keff in rxn_keff_dict.items():
#    rxn = me.reactions.get_by_id(rxn_name)
#    rxn.keff = keff
#    rxn.update()

# #genes that are more abundant than BOP27/BOP1000 experimental measurements
#rxn_keff_dict['']=
#rxn_keff_dict['']=
for rxn in me.reactions.query('S-ADENMETSYN-CPLX'): # metK, b2942
    if hasattr(rxn, 'keff'):
        #rxn.keff = 145.5284555731781/100. # 1/100 of sasa keff
        keff_old = rxn.keff
        rxn.keff = keff_old * 40. # 40x of david keff 
        rxn.update()

for rxn in me.reactions.query('EG11817-MONOMER'): # umpG, b2744
    if hasattr(rxn, 'keff'):
        #rxn.keff = 368.98535798686396 # 10x of the sasa keff, can be more accurately set for different reactions
        keff_old = rxn.keff
        rxn.keff = keff_old * 500. # 80x of david keff 
        rxn.update()

for rxn in me.reactions.query('TRYPTOPHAN-CPLX'): # tnaA, b3708
    if hasattr(rxn, 'keff'):
        #rxn.keff =  1727.3550547267016 # 10x of the sasa keff, can be more accurately set for different reactions
        keff_old = rxn.keff
        rxn.keff = keff_old * 100.  # 60x of the david keff
        rxn.update()

cplx_list=['CPLX0-7643','FABB-CPLX','EG12712-MONOMER','LUMAZINESYN-CPLX','3-METHYL-2-OXOBUT-OHCH3XFER-CPLX','ACYLCOADEHYDROG-MONOMER']
# ilvC, b3774
# fabB, b2323
# luxS, b2687
# ribE, b0415
# panB, b0134
# fadE, b0221
for cplx in cplx_list:
    for rxn in me.reactions.query(cplx):
        if hasattr(rxn, 'keff'):
            keff_old = rxn.keff
            rxn.keff = keff_old*10. # 10x of the david keff
            rxn.update()

cplx_dict={'AAS-MONOMER':30.,'CPLX0-1666':50.,'CPLX0-1667':360.,'FADB-CPLX':10.,'G6700-MONOMER':10000.}
# aas,  b2836
# fadJ, b2341
# fadI, b2342
# fadB, b3846
# ompN, b1377
cplx_dict['CPLX0-761']=20. # xdhABC, b2866-8
cplx_dict['KETOBUTFORMLY-INACT-MONOMER_mod_glycyl']=18. # tdcE, b3114
#cplx_dict['PROPKIN-MONOMER']=300. # tdcD, b3115
cplx_dict['THREDEHYDCAT-CPLX_mod_pydx5p']=1300. # tdcB, b3117
cplx_dict['EG10723-MONOMER']= 100. # phnN, b4094
cplx_dict['ASPCARBCAT-TRIMER']=10. #pyrB, b4245
cplx_dict['ASPCARBTRANS-CPLX_mod_zn2']=10. #pyrB, b4245 (in complex with b4244 pyrI)
cplx_dict['GLUTAMINESYN-OLIGOMER_mod_mn2']=15. #glnA, b3870
cplx_dict['GALP-MONOMER']=30. #galP, b2943
cplx_dict['G7493-MONOMER']=25. # b2874, yqeA
cplx_dict['CPLX0-7912']=1000. #eutD, b2458
cplx_dict['GLUCOKIN-MONOMER']=20. #glk, b2388
cplx_dict['G6960-MONOMER']=25. #ydjI, b1773 # FBA
cplx_dict['CARBODEHYDRAT-CPLX_mod_zn2']=300. # cynT, b0339
cplx_dict['GPT-CPLX']=17. #gpt, b0238
#cplx_dict['KDPGALDOL-4OH2OXOGLUTARALDOL-CPLX']=20. #eda, b1850
#cplx_dict['PGLUCONDEHYDRAT-MONOMER']=55. #edd, b1851
cplx_dict['CPLX0-303_mod_mn2']=30. # glpX, b3925
cplx_dict['ORNCARBAMTRANSFERI-CPLX']=15. #argI, b4254
cplx_dict['CPLX0-341']=10. #dfp, b3639

for cplx,foldchange in cplx_dict.items():
    for rxn in me.reactions.query(cplx):
        if hasattr(rxn, 'keff'):            
            keff_old = rxn.keff
            rxn.keff = keff_old * foldchange
            rxn.update()

#rxn=me.reactions.get_by_id('EDA_FWD_KDPGALDOL-4OH2OXOGLUTARALDOL-CPLX')
#rxn.keff = 72.92370998596122 # sasa keff, david=22.834830097677397
#rxn.update()
#
#rxn=me.reactions.get_by_id('EDD_FWD_PGLUCONDEHYDRAT-MONOMER')
#rxn.keff = 71.10047771061795 # sasa keff, david=63.6584626602816
#rxn.update()

rxn=me.reactions.get_by_id('SUCDi_FWD_SUCC-DEHASE_mod_3fe4s_mod_fad_mod_2fe2s_mod_4fe4s')
rxn.keff = 113.06527038180563 # sasa keff, david=42.4644256884733
rxn.update()

# To solve problem of too high acetate secretion and low aero-type
# change ACKr, PTAr, PGI reaction rate back to SASA keff (meaning that they will have the same FWD and REV rate)
for rxn in me.reactions.query('PGI_'):
    if hasattr(rxn, 'keff'):            
        rxn.keff = 115.23498588199664 #sasa
        rxn.update()

for rxn in me.reactions.query('G6PDH2r_'):
    if hasattr(rxn, 'keff'):
        rxn.keff = 63.59636766381115  #sasa
        rxn.update()

for rxn in me.reactions.query('PGL_'):
    if hasattr(rxn, 'keff'):            
        rxn.keff = 46.11897913497253 #sasa
        rxn.update()

for rxn in me.reactions.query('ACKr'):
    if '_FWD_' in rxn.id:
        rxn.keff = 52.65787290142323 #/2.60866230357314 # david=39.103402172382395
    elif '_REV_' in rxn.id:
        rxn.keff = 52.65787290142323 # david=102.007571188554
    rxn.update()

# the following two complexes catalyze PTAr, which is upstream of ACKr
for rxn in me.reactions.query('CPLX0-7912'):
    if hasattr(rxn, 'keff'):
        rxn.keff = 77.20522054049931/10. #sasa
        rxn.update()

for rxn in me.reactions.query('PHOSACETYLTRANS-CPLX'):
    if hasattr(rxn, 'keff'):
        rxn.keff = 311.323023418915/10. #sasa
        rxn.update()

#me.unmodeled_protein_fraction = 0.1
me.gam = 24.2 
me.ngam = 16.35

with open('Models/FoldME_new_sasa_david_17.pickle','wb') as f:
    pickle.dump(me,f)
