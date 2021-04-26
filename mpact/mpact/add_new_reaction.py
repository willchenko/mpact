# functions to add a new reaction in mpact_webapp

from mpact_webapp.add_new_metabolite import get_kegg_compound_information

[cpd_id, formula, r1, name, dblinks] = get_kegg_compound_information(kegg_id1)
[cpd_id, formula, r2, name, dblinks] = get_kegg_compound_information(kegg_id2)


common_rxns = [value for value in r1 if value in r2]
display_rxns = display_common_rxns(common_rxns)
display_rxns = np.vstack(display_rxns)

def display_common_rxns(common_rxns):
    base = "http://rest.kegg.jp/get/"
    i = 0
    display_rxns = []
    for rxn in common_rxns:
        url = base + rxn
        response = requests.get(url)
        tsv_data = response.text
        data = re.split('\n',tsv_data)
        defn = [s for s in data if "DEFINITION" in s]
        defn = re.split('DEFINITION ',defn[0])
        defn = [s for s in defn if s][0]
        display_rxns.append([i,defn])
        i = i + 1
    return display_rxns