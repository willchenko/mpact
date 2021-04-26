# fatty acid metadata calculation
from mpact import ModCellInput2Mpact
from mpact import write_rxn_string_from_mets_and_stoich
from mpact import parseRxnString
from mpact import metID2Name
import numpy as np
import re

metadata = get_FA_metadata(folder_path,metadata,metabolite_dict,reaction_dict,pathway_dict,mmdb_path)
writeMetaData(metadata,'metadata_table.csv',folder_path)


def get_FA_metadata(folder_path,metadata,metabolite_dict,reaction_dict,pathway_dict,mmdb_path):
    [metabolite_dict,reaction_dict,pathway_dict]= ModCellInput2Mpact(folder_path,metabolite_dict,reaction_dict,pathway_dict,mmdb_path)
    #[mets,stoich,rev] = parseRxnString('accoa_c + nad_c => bla_c')
    for path in pathway_dict['pathways']:
        fa = metID2Name(pathway_dict[path]['product_id'],metabolite_dict)
        fa_id = re.split(' ',fa)[1][1:-1]
        x = int(re.split(':',fa_id)[0][3:])
        y = int(re.split(':',fa_id)[1])
        n = (x-2)/2
        num_rxns = n*4 + y
        ov_mets = ['accoa_c','atp_c','nadph_c','h_c',metabolite_dict[fa]['bigg_id'],'coa_c','nadp_c','h2o_c']
        ov_stoich = list(np.zeros((len(ov_mets),), dtype=int))
        ov_stoich[ov_mets.index('accoa_c')] = -1*(n + 1)
        ov_stoich[ov_mets.index('atp_c')] = -1*n
        ov_stoich[ov_mets.index('nadph_c')] = -1*(2*n)
        ov_stoich[ov_mets.index('h_c')] = -2*n
        ov_stoich[ov_mets.index(metabolite_dict[fa]['bigg_id'])] = 1
        ov_stoich[ov_mets.index('coa_c')] = n + 1
        ov_stoich[ov_mets.index('nadp_c')] = 2*n
        ov_stoich[ov_mets.index('h2o_c')] = n - 1
        if y >= 1:
            ov_stoich[ov_mets.index('coa_c')] = ov_stoich[ov_mets.index('coa_c')] + 1
        if y >= 2:
            ov_stoich[ov_mets.index('nadph_c')] = ov_stoich[ov_mets.index('nadph_c')] - 1
            ov_stoich[ov_mets.index('nadp_c')] = ov_stoich[ov_mets.index('nadp_c')] + 1
        if y >= 3:
            ov_stoich[ov_mets.index('h2o_c')] = ov_stoich[ov_mets.index('h2o_c')] + 1
        if y >= 4:
            ov_stoich[ov_mets.index('nadph_c')] = ov_stoich[ov_mets.index('nadph_c')] - 1
            ov_stoich[ov_mets.index('nadp_c')] = ov_stoich[ov_mets.index('nadp_c')] + 1
        if 'P' in metabolite_dict[fa]['formula']:
            ov_mets.append('pi_c')
            ov_stoich.append(-1)
            ov_stoich[ov_mets.index('h_c')] = ov_stoich[ov_mets.index('h_c')] - 1
            num_rxns = num_rxns + 1
        ov_string = write_rxn_string_from_mets_and_stoich(ov_mets,ov_stoich)
        metadata[fa] = {}
        metadata[fa]['precursor'] = 'accoa'
        metadata[fa]['target'] = fa
        metadata[fa]['number_of_rxns'] = num_rxns
        metadata[fa]['ov_string'] = ov_string
        metadata['pathways'].append(fa)
    return metadata
    







    return
