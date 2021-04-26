# add metadata for pathways
from mpact import write_ov_path_eqn
from mpact import findMetNameFromID
from mpact import getOvStringWithPathsAsRxn

def add_to_metadata(path,metadata,pathway_dict,reaction_dict,metabolite_dict):
    precursor_id = pathway_dict[path]['precursor_ids']
    precursor = findMetNameFromID(precursor_id,metabolite_dict)
    #ov_string = write_ov_path_eqn(path,precursor,pathway_dict,reaction_dict,metabolite_dict)
    #num_rxns = len(pathway_dict[path]['rxns'])
    target_id = pathway_dict[path]['product_id']
    target = findMetNameFromID(target_id,metabolite_dict)
    metadata[path] = {}
    metadata[path]['precursor'] = precursor
    metadata[path]['target'] = target
    if 'paths' in list(pathway_dict[path].keys()):
        if pathway_dict[path]['paths'] == '0':
            ov_string = write_ov_path_eqn(path,precursor,pathway_dict,reaction_dict,metabolite_dict)
            num_rxns = len(pathway_dict[path]['rxns'])
        elif pathway_dict[path]['paths'] == '1':
            [ov_string,num_rxns] = getOvStringWithPathsAsRxn(path,pathway_dict,reaction_dict,metabolite_dict)
    else:
        ov_string = write_ov_path_eqn(path,precursor,pathway_dict,reaction_dict,metabolite_dict)
        num_rxns = len(pathway_dict[path]['rxns'])
    metadata[path]['number_of_rxns'] = num_rxns
    metadata[path]['ov_string'] = ov_string
    metadata['pathways'].append(path)
    return metadata

def add_all_paths_to_metadata(metadata,pathway_dict,reaction_dict,metabolite_dict):
    for path in pathway_dict['pathways']:
        metadata = add_to_metadata(path,metadata,pathway_dict,reaction_dict,metabolite_dict)
    return metadata
    
    
