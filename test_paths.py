import model_data
import paths


node1 = 'f6p_c'
node2 = 'h_c'

md = model_data.model_data()
md.load_model_cobra_sbml("e_coli_core.xml")
md.run_FBA()
"""
P = paths.paths(md, ignore_weights = 0, max_metabolite_occurance = 50)


print(P.get_shortest_path(node1, node2))

all_paths = P.get_all_shortest_paths(node1, node2)
all_paths = [" --> ".join([i for i in p if i in md.reactions]) for p in all_paths]
print(all_paths)
"""

P = paths.paths(md, ignore_weights = 0, max_metabolite_occurance = 50)
print(P.get_sub_graph([node1,node2]))

            


