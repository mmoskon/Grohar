#import numpy as np
import networkx as nx
#import matplotlib.pyplot as plt
import make_pairs
import SteinerTree as st


class paths:
    
    """
    # uporabiti matriko povezanosti reakcij
    # dilema: kje konstruirati matriko - konstruiramo jo že v alignmentu - ali jo bi bilo smiselno konstruirati kje drugje
    def get_reaction_connectivity(model, compartment = 'c', max_metabolite_occurence = 0, reversibility = 1):
        
        
        if compartment:
            reactions = list(model.compartment_reactions[compartment])            
        else:
            reactions = list(model.reactions)
        r = len(reactions)
        # matrika povezanosti reakcij
        M = np.zeros((r, r))
    
        metabolites_to_remove = set()
        if max_metabolite_occurence:
            for m in model.metabolites: #model.metabolite_reactions:
                # ven pomecemo prepogoste metabolite!!!
                if len(model.metabolite_reactions[m]) > max_metabolite_occurence:
                    metabolites_to_remove.add(m)
       
        reaction_names = []
        for (i1, r1) in enumerate(reactions):
            reaction_names.append(model.reactions[r1])
    
            for (i2, r2) in enumerate(reactions):
                if reversibility:
                    lb1, ub1 = model.reaction_bounds[r1]
                    lb2, ub2 = model.reaction_bounds[r2]
                        
                                        
                    if ub1 > 0:
                        if ub2 > 0:
                            if (model.reaction_products[r1] & model.reaction_reactants[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1
                               continue
                               #print("0: ",r1, r2)
                        if lb2 < 0 and r1 != r2:
                            if (model.reaction_products[r1] & model.reaction_products[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1                           
                               continue
                               #print("1: ",r1,r2)
                    
                    if lb1 < 0:
                        if ub2 > 0 and r1 != r2:
                            if (model.reaction_reactants[r1] & model.reaction_reactants[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1
                               continue
                               #print("3: ",r1, r2)
                        if lb2 < 0:
                            if (model.reaction_reactants[r1] & model.reaction_products[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1       
                               continue
                               #print("4: ",r1, r2)
                else:
                    if (model.reaction_products[r1] & model.reaction_reactants[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                       M[i1,i2] = 1
                        
        
    
        return (reactions, M)
    
    ###
    # NI OK - TREBA JE IZHAJATI IZ GRAFA ORIGINALNEGA MODELA, KER SICER NE MOREMO OPAZOVATI FLUKSOV!!!
    ###
    def shortest_path_reaction_graph(model, reaction1, reaction2, compartment = "c"):
        
        
        
        reactions, M = get_reaction_connectivity(model, compartment = compartment)
        #if reaction1 in reactions and reaction2 in reactions:
        if model.fluxes:
            weights = [model.fluxes[r] for r in reactions]
            print(weights)
        else:
            weights = [1 for r in reactions]
            print(weights)
            
        G = nx.from_numpy_matrix(M, create_using=nx.MultiDiGraph())
        
        pos = nx.circular_layout(G)
        nx.draw_circular(G)
        labels = {i : i + 1 for i in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=15)
        plt.show()
    
    """
    ########################################
    ########################################
    ########################################
    
    def __init__(self, model, ignore_weights, max_metabolite_occurance = 0, perturbed = 0):
        self.G = self.prepare_graph(model, ignore_weights, max_metabolite_occurance, perturbed)
    
    
    
    def prepare_graph(self, model, ignore_weights, max_metabolite_occurance, perturbed):
    
        if not perturbed:
            pairs, fluxes = make_pairs.make_pairs_reactions(model, model.reactions, max_metabolite_occurance = max_metabolite_occurance, change_reversed = 1)
        else:
            pairs, fluxes = make_pairs.make_pairs_reactions_perturbation(model, model.reactions, max_metabolite_occurance = max_metabolite_occurance, change_reversed = 1)
        G = nx.DiGraph() 
        
        
        if fluxes:
            M = max(fluxes.values()) + 1
        
        for a,b in pairs:
            if fluxes and not ignore_weights:
                if a in fluxes:         
                    if fluxes[a]:
                        weight = M - abs(fluxes[a])
                    else:
                        weight = 100000
                    #weight = abs(fluxes[a])
                elif b in fluxes:
                    if fluxes[b]:
                        weight = M - abs(fluxes[b])
                    else:
                        weight = 100000
                    #weight = abs(fluxes[b])
                else:
                    weight = M
            else:
                weight = 1
                
            a = a.replace(":","---")
            b = b.replace(":","---")
                
            G.add_edge(a, b, weight = weight)
        
        return G
    
    def get_all_shortest_paths(self, node1, node2):
               
        
        node1 = node1.replace(":","---")
        node2 = node2.replace(":","---")
        
        
      
        paths = nx.all_shortest_paths(self.G, source = node1, target = node2, weight = 'weight')
        #for p in paths:
        #    print(p)
        
        
        
        return paths
        
    
    
    
    def get_shortest_path(self, node1, node2):
               
        
        node1 = node1.replace(":","---")
        node2 = node2.replace(":","---")
        print(node1)
        print(node2)
        
        path = nx.shortest_path(self.G, source = node1, target = node2, weight = 'weight')        
        
        print(path)
        
        return path
        
    def get_sub_graph(self, nodes):
        G_S = st.make_steiner_tree(self.G, nodes)
        path = [n for n in G_S.nodes()]
    
        return(path)
