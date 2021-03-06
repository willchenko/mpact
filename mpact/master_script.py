import os
mpact_path = r"/mnt/c/Users/Owner/Documents/PhD_stuff/research/mpact"
os.chdir(mpact_path)
import mpact
from mpact import *

# add metabolites
compounds = ['Glutarate']
[metabolite_dict,bigg2kegg_db] = mpact.add_to_metabolite_dict(metabolite_dict,compounds,bigg2kegg_db,mpact_path)

# add reactions
compounds = ["Glutarate","Glutarate semialdehyde"]
[reaction_dict, metabolite_dict] = mpact.add_to_reaction_dict(reaction_dict,compounds,metabolite_dict,bigg2kegg_rxn_db,bigg2kegg_db,mpact_path)

# solve for unknown charges
metabolite_dict = mpact.solveForMissingChargesinRecursiveLoop(metabolite_dict,reaction_dict)

# check balance
reaction_dict = mpact.mass_check_reaction_dict(reaction_dict,metabolite_dict)
reaction_dict = mpact.charge_check_reaction_dict(reaction_dict,metabolite_dict)
mpact.find_unbalanced_rxns(reaction_dict)
mpact.rxnChargeBalance('HP',reaction_dict,metabolite_dict)
mpact.check_rxn_mass_sum('R5PDIP',reaction_dict,metabolite_dict)

# find and add all pathways between two metabolites
precursor = 'Glutarate'
target = 'Glutarate semialdehyde'
[pathway_dict,reaction_dict,metabolite_dict] = mpact.add_to_pathway_dict(pathway_dict,precursor,target,reaction_dict,metabolite_dict,mpact_path)

# ModCell input --> mpact
folder_path = r"C:\Users\Owner\Documents\PhD_stuff\research\ModCell2-master\ModCell2-master\problems\ecoli-gem-part1\input"
[metabolite_dict,reaction_dict,pathway_dict]= mpact.ModCellInput2Mpact(folder_path,metabolite_dict,reaction_dict,pathway_dict,mpact_path)


# add data to local db
mpact.add_all_data_to_local_db(mpact_path,metabolite_dict,reaction_dict,pathway_dict)

### retrieve all data from local db
# retrieve all data from a target's pathway
target = 'Glutarate'
[metabolite_dict,reaction_dict,pathway_dict] = mpact.add_all_info_of_paths_that_make_target(target,mpact_path,metabolite_dict,reaction_dict,pathway_dict)

#retrieve one path
path = 'cinnamic_acid_D-Erythrose-4-phosphate_0'
[metabolite_dict,reaction_dict,pathway_dict] = mpact.get_all_info_from_db_path(path,pathway_dict,reaction_dict,metabolite_dict,mpact_path)


# retrieve all reactions 
reaction_dict = mpact.get_all_reactions_from_db(mpact_path,reaction_dict)
metabolite_dict = mpact.get_all_metabolites_from_db(mpact_path,metabolite_dict)

# add metadata
metadata = mpact.add_all_paths_to_metadata(metadata,pathway_dict,reaction_dict,metabolite_dict)

mpact.writeMetaData(metadata,'metadata_table.csv',folder_path)


#manual metabolite addition
metabolite_dict = mpact.manually_add_to_metabolite_dict(metabolite_dict,'Isobutyl acetate')
reaction_dict = mpact.manually_add_to_reaction_dict('AATibutylace',reaction_dict)

# manually add pathways to pathway_dict
precursor = '2-oxoglutarate'
product = 'L-hydroxyproline'
pathway_rxns = ['GLUDxi','ACOTA','ORNTAC_1','ORNCD','PROAKGOX1']
pathway_name = 'L-hydroxyproline 1'
pathway_dict = mpact.add_to_pathway_dict(pathway_dict,pathway_name,product,pathway_rxns,precursor,metabolite_dict)


## exporting

#excel
# create a new problem
problem_name = 'adpac-glutarate'
mpact.create_new_problem(problem_name,mpact_path)
# add 
mpact.add_to_pathway_csv_table(problem_name,pathway_dict, mpact_path)
mpact.add_to_reaction_csv_table(problem_name,reaction_dict, mpact_path)
mpact.add_to_metabolite_csv_table(problem_name,metabolite_dict, mpact_path)

# json
mpact.export2json('e4p_fm.json',reaction_dict,metabolite_dict)

                                                                                                                                                                                  
