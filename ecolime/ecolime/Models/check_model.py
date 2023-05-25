import pickle
import argparse
import pandas as pd

# input
parser = argparse.ArgumentParser(description='Check Folding Parameters')
parser.add_argument('model',type=str, nargs=1,help='model file name')
args = parser.parse_args()

model = args.model
model = model[0]

me=pickle.load(open(model,'rb'))

#rxn_keff_dict={'MPTSS1_FWD_EG10154-MONOMER':36.650302321444784, 'BMOCOS_FWD_EG10153-MONOMER':53.34083665106091, 'MOCOS_FWD_EG10153-MONOMER':53.34083665106091}
#rxn_keff_dict['CPMPS_FWD_EG11595-MONOMER_EG11666-MONOMER']=62.83982200999657 # b0781,b0783
#rxn_keff_dict['MOADSUx1_FWD_CPLX_dummy']=43.640419650213026 # b0784 moaD
#rxn_keff_dict['GLYCL_FWD_GCVMULTI-CPLX_mod_lipo']=260.0390765173197 # gcvHPT
#rxn_keff_dict['SUCDi_FWD_SUCC-DEHASE_mod_3fe4s_mod_fad_mod_2fe2s_mod_4fe4s']=113.06527038180563 # sdhABCD
#rxn_keff_dict['']=
#rxn_keff_dict['']=

#for rxn_name in rxn_keff_dict.keys():
#    rxn = me.reactions.get_by_id(rxn_name)
#    print (rxn,id,rxn.keff)
#
for rxn in me.reactions.query('PGL_'): 

    if hasattr(rxn, 'keff'):
        print (rxn.id,rxn.keff)
