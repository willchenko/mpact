# combine paths
from mpact import write_ov_path_eqn
from mpact import write_rxn_string_from_mets_and_stoich
from mpact import parseRxnString
from mpact import metID2Name
import numpy as np
import re

# if paths = 1
#path = 'Ethyl acetate_Pyruvate_0'
#[ov_string,num_rxns] = getOvStringWithPathAsRxn(path,pathway_dict,reaction_dict,metabolite_dict)

def getOvStringWithPathsAsRxn(path,pathway_dict,reaction_dict,metabolite_dict):
    ov_string = ''
    num_rxns = 0 
    for rxn in pathway_dict[path]['rxns']:
        if rxn in reaction_dict['reactions']:
            #p.append(1)
            rxn_string = reaction_dict[rxn]['bigg_string']
            num_rxns = num_rxns + 1
        else:
            #p.append(0)
            rxn_string = write_ov_path_eqn(rxn,metID2Name(pathway_dict[rxn]['precursor_ids'],metabolite_dict),pathway_dict,reaction_dict,metabolite_dict)
            num_rxns = num_rxns + len(pathway_dict[rxn]['rxns'])
        if ov_string == '':
            ov_string = rxn_string
        else:
            ov_string = combineRxnStrings(ov_string,rxn_string)
    return ov_string, num_rxns

def combineRxnStrings(str1,str2):
    # example inputs
    #   str1 = 'pyr_c + nadh_c + 2 h_c => etoh_c + co2_c + nad_c'
    #   str2 = 'etoh_c + accoa_c =>  etylace_c + coa_c'
    [mets1,stoich1,rev1] = parseRxnString(str1)
    [mets2,stoich2,rev2] = parseRxnString(str2)
    ov_mets = list(set(mets1+mets2))
    ov_stoich = list(np.zeros((len(ov_mets),), dtype=int))
    for met in mets1:
        met_stoich = stoich1[mets1.index(met)]
        ov_met_index = ov_mets.index(met)
        ov_stoich[ov_met_index] = ov_stoich[ov_met_index] + met_stoich
    for met in mets2:
        met_stoich = stoich2[mets2.index(met)]
        ov_met_index = ov_mets.index(met)
        ov_stoich[ov_met_index] = ov_stoich[ov_met_index] + met_stoich
    ov_string = write_rxn_string_from_mets_and_stoich(ov_mets,ov_stoich)
    return ov_string
