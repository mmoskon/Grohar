#from cobra import Model, Reaction, Metabolite, io
import model_data as md
import make_pairs
import networkx as nx
import matplotlib.pyplot as plt






model = md.model_data()
#model.load_model_sbml("iCHOv1_final.xml")
#model.load_model_cobra_sbml("iCHOv1_final.xml")
#model.load_model_mat("gimmeDG44.mat")
model.load_model_mat("gimmeS.mat")
#model.load_model_mat("gimmeK1.mat")

model.run_FBA()
model.run_perturbation([['EX_gln_L_e_', 'l', 0], ['EX_asn_L_e_', 'l', 0]], objective_id = "biomass_cho_producing")
#model.run_perturbation(objective_id = "ASNN")

#model.set_boundaries(['EX_gln_L_e_', 'l', 0], ['EX_asn_L_e_', 'l', 0], objective_id = "biomass_cho_producing")

#met_name = "M_asn_L_c" # SBML
#met_name = "asn_L_c" # COBRA SBML
#met_names = {"asn_L[c]"}#, "nh4[c]"} # Matlab SBML

              
met_names = {"asn_L[c]", "gln_L[c]"}
# analiza laktata        
#met_names = {'lac_D[c]', 'lac_L[c]', 'lac_D[e]', 'lac_L[e]'}
             
             
#react_name = "R_ASNN" # SBML
#react_name = "ASNN" # COBRA
#print(model.reactions[react_name])





reaction_types = 0  # producing = 1, consuming = -1, all = 0

#(pairs, fluxes) = make_pairs.make_pairs(reaction_types, model, met_name, 1, compartment = 'c e', max_react_size = 25, max_metabolite_occurance = 150, min_flux = 0.0001)


comparison_type = 1 # new reactions = 1, reaction removed = -1, reactions in both model= 0
distance = 1

(pairs, fluxes) = make_pairs.make_pairs_perturbation(comparison_type, reaction_types, model, met_names, distance, compartment = 'c', max_react_size = 25, max_metabolite_occurance = 150, min_flux = 0.00000001)
#(pairs, fluxes) = make_pairs.make_pairs(reaction_types, model, met_names, distance, compartment = '', max_react_size = 25, max_metabolite_occurance = 150, min_flux = 0.00000001)

G = nx.DiGraph()
G.add_edges_from(pairs)



producing = set(make_pairs.direct_producers(model, met_names, fluxes))
consuming = set(make_pairs.direct_consumers(model, met_names, fluxes))

val_map = {'reaction': 'b',
           'metabolite': 'c',
           'producing': 'g',
           'both': 'y',
           'consuming': 'r',
           'observed': 'm'}

node_colors = [val_map['observed'] if node in met_names else val_map['both'] if node in producing & consuming else val_map['producing'] if node in producing else val_map['consuming'] if node in consuming else val_map['reaction'] if node in model.reactions else val_map['metabolite'] for node in G.nodes()]
                             
node_labels = {a:a + ":" + "{:.5f}".format(fluxes[a]) if a in model.reactions else a for a in G.nodes()}
node_sizes = [500 if a in model.reactions else 250 for a in G.nodes()]


#pos = nx.nx_pydot.graphviz_layout(G, prog="neato")              
#pos = nx.nx_pydot.graphviz_layout(G, prog="twopi")              
#pos = nx.nx_pydot.graphviz_layout(G, prog="circo")              
pos = nx.nx_pydot.graphviz_layout(G, prog="dot")              
#pos = nx.spring_layout(G)
#pos = nx.spectral_layout(G)

plt.clf()


nx.draw_networkx_labels(G, pos , node_labels, font_size=8, font_weight = "bold", font_color = "#1c2833")
#nx.draw_networkx_labels(G,labels=node_labels)

nx.draw(G, pos, node_color = node_colors, node_size=node_sizes, edge_color = "#abb2b9", dge_cmap=plt.cm.Reds, alpha = 0.4)
#nx.draw_networkx_nodes(G, pos, nodelist = [a for a in G if a in model.reactions], node_color = node_colors, node_size=250,  alpha=0.4)
#nx.draw_networkx_nodes(G, pos, nodelist = [a for a in G if a in model.metabolites], node_color = node_colors, node_size=100,  alpha=0.4)
#nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)

#pylab.show()


plt.axis('off')
plt.savefig("labels_and_colors.png") # save as png
plt.show() # display