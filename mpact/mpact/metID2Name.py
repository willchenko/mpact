# met ID to name

def metID2Name(ID,metabolite_dict):
    IDs = []
    for met in metabolite_dict['metabolites']:
        IDs.append(metabolite_dict[met]['bigg_id'])
    name = metabolite_dict['metabolites'][IDs.index(ID)]    
    return name
