# ModCell input folder --> mpact
import os
import csv
import re
from mpact import get_list_of_mets_in_db
from mpact import retrieve_metabolite_info_from_db
from mpact import metID2Name
from mpact import solveForMissingChargesinRecursiveLoop


def ModCellInput2Mpact(path,metabolite_dict,reaction_dict,pathway_dict,mmdb_path):
    os.chdir(path)
    metabolite_table = get_table(path,'metabolite_table.csv')
    metabolite_dict = get_metabolites_from_table(mmdb_path,metabolite_table,metabolite_dict)
    reaction_table = get_table(path,'reaction_table.csv')
    node_table = get_table(path,'node_table.csv')
    reaction_dict = get_reactions_from_table(mmdb_path,reaction_table,node_table,metabolite_dict,reaction_dict)
    pathway_table = get_table(path,'pathway_table.csv')
    pathway_dict = get_pathways_from_table(mmdb_path,pathway_table,pathway_dict,metabolite_dict)
    metabolite_dict = solveForMissingChargesinRecursiveLoop(metabolite_dict,reaction_dict)
    metabolite_dict = get_neutral_formulas(metabolite_dict)
    return metabolite_dict,reaction_dict,pathway_dict

def get_metabolites_from_table(mmdb_path,metabolite_table,metabolite_dict):
    db_mets = get_list_of_mets_in_db(mmdb_path)
    for row in metabolite_table:
        if row[0] != 'id':
            name = row[1]
            if name in db_mets:
                metabolite_dict = retrieve_metabolite_info_from_db(row[1],mmdb_path,metabolite_dict)
            else:
                metabolite_dict[name] = {}
                metabolite_dict[name]['kegg_id'] = row[4]
                metabolite_dict[name]['bigg_id'] = row[0]+'_c'
                if row[3] == '':
                    charge = 0
                else:
                    charge = int(row[3])
                metabolite_dict[name]['charge'] = charge
                metabolite_dict[name]['formula'] = row[2]
                metabolite_dict[name]['bigg_formula'] = row[2]
                metabolite_dict[name]['bigg_charges'] = [charge]
                metabolite_dict['metabolites'].append(name)
    return metabolite_dict

def get_reactions_from_table(mmdb_path,reaction_table,node_table,metabolite_dict,reaction_dict):
    node_rxns = [i[0] for i in node_table]
    for row in reaction_table:
        if row[0] != 'id':
            node_index = node_rxns.index(row[0])
            #node1 = metID2Name(node_table[node_index][1]+'_c',metabolite_dict)
            node1 = node_table[node_index][1] + '_c'
            #node2 = metID2Name(node_table[node_index][2]+'_c',metabolite_dict)
            node2 = node_table[node_index][2] + '_c'
            nodes = [node1,node2]
            reaction_dict[row[0]] = {}
            reaction_dict[row[0]]['bigg_rxn'] = row[4]
            reaction_dict[row[0]]['bigg_name'] = row[1]
            reaction_dict[row[0]]['kegg_rxn'] = row[3]
            reaction_dict[row[0]]['bigg_string'] = row[2]
            reaction_dict[row[0]]['BRENDA_enzymes'] = row[5]
            reaction_dict[row[0]]['nodes'] = nodes
            reaction_dict['reactions'].append(row[0])
    return reaction_dict

def get_pathways_from_table(mmdb_path,pathway_table,pathway_dict,metabolite_dict):
    for row in pathway_table:
        if row[0] != 'id':
            target = metID2Name(row[3]+'_c',metabolite_dict)
            precursor = row[5]
            precursor = metID2Name(precursor+'_c',metabolite_dict)
            i = get_pathway_naming_index(pathway_dict,target,precursor,metabolite_dict)
            name = target+'_'+precursor+'_'+str(i)
            pathway_dict[name] = {}
            pathway_dict['pathways'].append(name)
            pathway_dict[name]['id'] = name
            pathway_dict[name]['name'] = name
            rxns = re.split("', '",row[2][2:-2])
            pathway_dict[name]['rxns'] = rxns
            pathway_dict[name]['paths'] = row[6]
            pathway_dict[name]['product_id'] = metabolite_dict[target]['bigg_id']
            pathway_dict[name]['precursor_ids'] = metabolite_dict[precursor]['bigg_id']
    return pathway_dict

def get_table(path,filename):
    table = []
    os.chdir(path)
    with open(filename, 'r') as csvfile: 
        csvreader = csv.reader(csvfile) 
        for row in csvreader: 
            table.append(row)
    return table

def get_all_dict_paths_between_mets(target,precursor,pathway_dict,metabolite_dict):
    precursor_id = metabolite_dict[precursor]['bigg_id']
    target_id = metabolite_dict[target]['bigg_id']
    paths = []
    for path in pathway_dict['pathways']:
        if pathway_dict[path]['product_id'] == target_id and pathway_dict[path]['precursor_ids'] == precursor_id:
            paths.append(path)
    return paths

def get_pathway_naming_index(pathway_dict,target,precursor,metabolite_dict):
    paths = get_all_dict_paths_between_mets(target,precursor,pathway_dict,metabolite_dict)
    if paths == []:
        current_indices = [-1]
    else:    
        current_indices = []
        for path in paths:
            current_indices.append(int(re.split('_',path)[2]))
    return max(current_indices) + 1

def get_neutral_formula(charged_formula,charge):
    if charge == 'NA':
        neutral_formula = charged_formula
    else:
        elements = " ".join(re.split("[^a-zA-Z]*", charged_formula))
        elements = elements.replace(' ','')
        [h_charged_coeff,h_start_index,h_end_index] = get_element_coeff(charged_formula,'H')
        if h_charged_coeff == 0:
            neutral_formula = charged_formula
        #elif h_charged_coeff == 1:
        #    charged_formula = neutral_formula.replace('H','')
        else:
            h_neutral_coeff = int(h_charged_coeff - charge)
            neutral_formula = charged_formula[0:h_start_index] + str(h_neutral_coeff) + charged_formula[h_end_index:len(charged_formula)]
    return neutral_formula

def get_neutral_formulas(metabolite_dict):
    for met in metabolite_dict['metabolites']:
        charge = metabolite_dict[met]['charge']
        charged_formula = metabolite_dict[met]['formula']
        neutral_formula = get_neutral_formula(charged_formula,charge)
        metabolite_dict[met]['neutral_formula'] = neutral_formula
    return metabolite_dict

def get_element_coeff(formula,element):
    elements = " ".join(re.split("[^a-zA-Z]*", formula))
    elements = elements.replace(' ','')
    if elements == element or element not in elements:
        coeff = 0
        start_index = 0
        end_index = 0
    else:
        element_index = elements.index(element)
        element_index_in_formula = formula.index(element)
        start_index = element_index_in_formula + 1
        if element_index == len(elements)-1:
            end_index = len(formula)
            coeff = formula[start_index:end_index]
        else:
            element_after_h = elements[element_index + 1]
            end_index = formula.index(element_after_h)
            coeff = formula[start_index:end_index]
        if coeff == '':
            coeff = 1
    return int(coeff), start_index, end_index
        
