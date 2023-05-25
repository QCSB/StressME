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
modelfile = 'Models/FoldME_new_sasa_Q40.pickle'
with open(modelfile,'rb') as iofile:
    me = pickle.load(iofile)

# from "kappmax_KO_ALE_davidi_per_pp_per_s_repl"
rxn_keff_dict={'PTAr_REV_CPLX0-7912':331.303239232142,
        'PTAr_FWD_CPLX0-7912':651.39322018908,
        'PTAr_REV_PHOSACETYLTRANS-CPLX':993.909717696426,
        'PTAr_FWD_PHOSACETYLTRANS-CPLX':1954.17966056724,
        'ACKr_REV_ACETATEKINA-MONOMER_mod_mg2':320.762753520128,
        'ACKr_FWD_ACETATEKINA-MONOMER_mod_mg2':154.296915465516,
        'ACKr_REV_PROPKIN-MONOMER':320.762753520128,
        'ACKr_FWD_PROPKIN-MONOMER':154.296915465516,
        'ACKr_REV_GARTRANSFORMYL2-MONOMER':320.762753520128,
        'ACKr_FWD_GARTRANSFORMYL2-MONOMER':154.296915465516}
# from "kappmax_davidi_per_pp_per_s_repl.pickle"
#rxn_keff_dict={'PTAr_REV_CPLX0-7912':168.62458785879141,
#        'PTAr_FWD_CPLX0-7912':273.703119195478,
#        'PTAr_REV_PHOSACETYLTRANS-CPLX':505.87376357637424,
#        'PTAr_FWD_PHOSACETYLTRANS-CPLX':821.109357586434,
#        'ACKr_REV_ACETATEKINA-MONOMER_mod_mg2':175.52702777854302,
#        'ACKr_FWD_ACETATEKINA-MONOMER_mod_mg2':91.3463434596863,
#        'ACKr_REV_PROPKIN-MONOMER':175.52702777854302,
#        'ACKr_FWD_PROPKIN-MONOMER':91.3463434596863,
#        'ACKr_REV_GARTRANSFORMYL2-MONOMER':175.52702777854302,
#        'ACKr_FWD_GARTRANSFORMYL2-MONOMER':91.3463434596863}
# from "kappmax_KO_ALE_per_pp_per_s_repl"
#rxn_keff_dict={'PTAr_REV_CPLX0-7912':55.112817615189,
#        'PTAr_FWD_CPLX0-7912':92.89185306870961,
#        'PTAr_REV_PHOSACETYLTRANS-CPLX':165.33845284556702,
#        'PTAr_FWD_PHOSACETYLTRANS-CPLX':278.6755592061288,
#        'ACKr_REV_ACETATEKINA-MONOMER_mod_mg2':139.780070170674,
#        'ACKr_FWD_ACETATEKINA-MONOMER_mod_mg2':57.6480788258165,
#        'ACKr_REV_PROPKIN-MONOMER':139.780070170674,
#        'ACKr_FWD_PROPKIN-MONOMER':57.6480788258165,
#        'ACKr_REV_GARTRANSFORMYL2-MONOMER':139.780070170674,
#        'ACKr_FWD_GARTRANSFORMYL2-MONOMER':57.6480788258165}

for rxn,keff in rxn_keff_dict.items():
    reaction = me.reactions.get_by_id(rxn)
    reaction.keff = keff
    reaction.update()

#me.unmodeled_protein_fraction = 0.55
#me.gam = 36.9
#me.ngam = 17.5

me.unmodeled_protein_fraction = 0.50
me.gam = 62.0
me.ngam = 0.

with open('Models/FoldME_new_sasa_PTAr_ACKr_Q55.pickle','wb') as f:
    pickle.dump(me,f)
