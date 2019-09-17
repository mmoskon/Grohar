#from model_data import *
#from FBA_data import *


def get_metabolite_occurrence(model, M):
    simple = 1
    if simple or not model.fluxes:
        return len(model.metabolite_reactions[M])
    else:
        c = 0
        for r in model.metabolite_reactions[M]:
            if model.fluxes[r] != 0:
                c += 1
        return c
        
        
        



def get_reactions_connected(model, metabolite, dist, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, perturbation = 0):
    tmp_metabolites = metabolite.copy()
    reactions = set()
    
    if model.fluxes:
        if not perturbation:
            fluxes = model.fluxes
        else:
            fluxes = model.perturbed_fluxes
    
    for d in range(dist):
        metabolites = tmp_metabolites.copy()
        tmp_metabolites.clear()
        for M in metabolites:
            for R in model.metabolite_reactions[M]:
                if not compartment or R in model.compartment_reactions[compartment]:    # compratment nastavljen
                    if not max_react_size or len(model.reaction_metabolites[R]) <= max_react_size:  # max reaction size
                        if not min_flux or abs(fluxes[R]) >= min_flux:    # min flux
                            reactions.add(R)                
                            if not max_metabolite_occurance:    # max metabolite occurance nastavljen ali ne
                                tmp_metabolites |= model.reaction_metabolites[R]
                            else:
                                for MM in model.reaction_metabolites[R]:
                                    if len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                     
    return reactions                                        

    
def get_reactions_producing(model, metabolite, dist, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, ignore_fluxes = 0, perturbation = 0):    
    tmp_metabolites = metabolite.copy()
    reactions = set()
     
    if model.fluxes:
        if not perturbation:
            fluxes = model.fluxes
        else:
            fluxes = model.perturbed_fluxes
       
    for d in range(dist):
        metabolites = tmp_metabolites.copy()      
        tmp_metabolites.clear()
        for M in metabolites:
            for R in model.metabolite_reactions[M]:
                if not compartment or R in model.compartment_reactions[compartment]:    # compratment nastavljen
                    if not max_react_size or len(model.reaction_metabolites[R]) <= max_react_size:  # max reaction size
                        if not min_flux or abs(fluxes[R]) >= min_flux:    # min flux                            
                            if (ignore_fluxes or fluxes[R] >= 0) and (R in model.metabolite_as_product[M]): 
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati reaktante
                                for MM in model.reaction_reactants[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                                                             
                            if (not ignore_fluxes and fluxes[R] < 0 and R in model.metabolite_as_reactant[M]):
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati produkte
                                for MM in model.reaction_products[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)           
    return reactions                                        

def get_reactions_consuming(model, metabolite, dist, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, ignore_fluxes = 0, perturbation = 0):    
    tmp_metabolites = metabolite.copy()
    reactions = set()
            
    if model.fluxes:
       if not perturbation:
           fluxes = model.fluxes
       else:
           fluxes = model.perturbed_fluxes
    
    for d in range(dist):
        metabolites = tmp_metabolites.copy()
        tmp_metabolites.clear()
        for M in metabolites:
            for R in model.metabolite_reactions[M]:
                if not compartment or R in model.compartment_reactions[compartment]:    # compratment nastavljen
                    if not max_react_size or len(model.reaction_metabolites[R]) <= max_react_size:  # max reaction size
                        if not min_flux or abs(fluxes[R]) >= min_flux:    # min flux
                            if (ignore_fluxes or fluxes[R] >= 0) and (R in model.metabolite_as_reactant[M]): 
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati produkte
                                for MM in model.reaction_products[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                                                             
                            if (not ignore_fluxes and fluxes[R] < 0 and R in model.metabolite_as_product[M]):
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati reaktante
                                for MM in model.reaction_reactants[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)      
    return reactions                                        
                                        
def make_pairs(reaction_types, model, metabolite, dist = 1, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, change_reversed = 0, return_pairs = True):
    pairs = []
    
    if reaction_types == 2: 
        reactions1 = get_reactions_producing(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)
        reactions2 = get_reactions_consuming(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)
        reactions = reactions1 & reactions2        
    
        
        
    elif reaction_types == 1:
        reactions = get_reactions_producing(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)              
    elif reaction_types == -1:
        reactions = get_reactions_consuming(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)              
    else:    
        reactions = get_reactions_connected(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, perturbation = 0)              
  
    
    if not return_pairs:
        return (reactions, model.fluxes)

    
    
    for R in reactions:
        for M in model.reaction_products[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                if change_reversed:
                    if model.fluxes[R] >= 0:
                        pairs.append((R, M))
                    else:
                        pairs.append((M, R))
                else:
                    pairs.append((R, M))
                
        for M in model.reaction_reactants[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                if change_reversed:
                    if model.fluxes[R] >= 0:
                        pairs.append((M, R))
                    else:
                        pairs.append((R, M))
                else:
                    pairs.append((M, R))
                
                              
    fluxes = model.fluxes
  
      
    

              
    return (pairs, fluxes)

def make_pairs_reactions(model, reactions, max_metabolite_occurance = 0, change_reversed = 0, perturbation = 0):
    pairs = []
    
    if not perturbation:
        fluxes = model.fluxes
    else:
        fluxes = model.perturbed_fluxes
    
    for R in reactions:
        for M in model.reaction_products[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                if fluxes and change_reversed:
                    if fluxes[R] >= 0:
                        pairs.append((R, M))
                    else:
                        pairs.append((M, R))
                else:
                    pairs.append((R, M))
                
        for M in model.reaction_reactants[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                if fluxes and change_reversed:
                    if fluxes[R] >= 0:
                        pairs.append((M, R))
                    else:
                        pairs.append((R, M))
                else:
                    pairs.append((M, R))
                
                              
    
    return (pairs, fluxes)
    
def make_pairs_reactions_perturbation(model, reactions, max_metabolite_occurance = 0, change_reversed = 0):
    return make_pairs_reactions(model, reactions, max_metabolite_occurance, change_reversed, perturbation = 1)


    
def make_pairs_perturbation(comparison_type, reaction_types, model, metabolite, dist = 1, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, return_pairs = True):    
   
    if comparison_type == 2:
        min_flux = 0
        
        
    if reaction_types == 2:
        reactions_p1 = get_reactions_producing(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)              
        reactions_p2 = get_reactions_producing(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 1)              
        
        reactions_c1 = get_reactions_consuming(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)              
        reactions_c2 = get_reactions_consuming(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 1)              
        
        reactions1 = reactions_p1 & reactions_c1
        reactions2 = reactions_p2 & reactions_c2
    elif reaction_types == 1: # producing
        reactions1 = get_reactions_producing(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)              
        reactions2 = get_reactions_producing(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 1)              
    elif reaction_types == -1: # consuming
        reactions1 = get_reactions_consuming(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 0)              
        reactions2 = get_reactions_consuming(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, ignore_fluxes = 0, perturbation = 1)              
    else:    # all
        reactions1 = get_reactions_connected(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, perturbation = 0)              
        reactions2 = get_reactions_connected(model, metabolite, dist, compartment, max_react_size, max_metabolite_occurance, min_flux, perturbation = 1)              
    
    return compare_perturbation_reactions(comparison_type, model, reactions1, reactions2, max_metabolite_occurance, return_pairs)      
    
   

def compare_perturbation_reactions(comparison_type, model, reactions1, reactions2, max_metabolite_occurance = 0, return_pairs = True):
    pairs = []
    
    if comparison_type == 1: # reactions that are new in the perturbed model
        reactions = reactions2 - reactions1
        fluxes = model.perturbed_fluxes
    elif comparison_type == -1: # reactions that are removed in the perturbed model
        reactions = reactions1 - reactions2
        fluxes = model.fluxes
    elif comparison_type == 2: # reactions that are present in both models but with modified flux values
        reactions = set()
        for r1 in reactions1:
            if r1 in reactions2 and model.fluxes[r1] != model.perturbed_fluxes[r1]:
                reactions.add(r1)                
        fluxes = {i:(model.fluxes[i], model.perturbed_fluxes[i]) for i in model.perturbed_fluxes}
    elif comparison_type == 3: # reactions that are different in both models
        reactions = set()
        for r1 in reactions1:
            if model.fluxes[r1] != model.perturbed_fluxes[r1]:
                reactions.add(r1)                
        for r2 in reactions2:
            if model.fluxes[r2] != model.perturbed_fluxes[r2]:
                reactions.add(r2)                                                            
        #fluxes = {i:(model.fluxes[i], model.perturbed_fluxes[i]) for i in reactions}
        #print(fluxes)
        fluxes = {i:(model.fluxes[i], model.perturbed_fluxes[i]) for i in model.perturbed_fluxes}
    elif comparison_type == 4: #reactions in both!        
        reactions = set(reactions1) | set(reactions2)
        fluxes = {i:(model.fluxes[i], model.perturbed_fluxes[i]) for i in model.perturbed_fluxes}    
    else: # reactions that are present in both models   
        reactions = set(reactions1) & set(reactions2)
        #fluxes = {i:model.perturbed_fluxes[i] - model.fluxes[i] for i in model.fluxes}
        fluxes = model.perturbed_fluxes #{i:model.perturbed_fluxes[i] for i in model.perturbed_fluxes}
        
    if not return_pairs:
        return(reactions, fluxes)
        
    for R in reactions:
        for M in model.reaction_products[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                pairs.append((R, M))
                
        for M in model.reaction_reactants[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                pairs.append((M, R))                 
                
    return (pairs, fluxes)
    
    
    
def make_pairs_producing(model, metabolite, dist, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, ignore_fluxes = 1):
    pairs = []
    tmp_metabolites = set([metabolite])
    reactions = set()
     
    if model.fluxes:
        fluxes = model.fluxes
       
    for d in range(dist):
        metabolites = tmp_metabolites.copy()
        tmp_metabolites.clear()
        for M in metabolites:
            for R in model.metabolite_reactions[M]:
                if not compartment or R in model.compartment_reactions[compartment]:    # compratment nastavljen
                    if not max_react_size or len(model.reaction_metabolites[R]) <= max_react_size:  # max reaction size
                        if not min_flux or abs(fluxes[R]) >= min_flux:    # min flux                            
                            if (ignore_fluxes or fluxes[R] >= 0) and (R in model.metabolite_as_product[M]): 
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati reaktante
                                for MM in model.reaction_reactants[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                                                             
                            if (not ignore_fluxes and fluxes[R] < 0 and R in model.metabolite_as_reactant[M]):
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati produkte
                                for MM in model.reaction_products[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                                       
    for R in reactions:
        for M in model.reaction_products[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                pairs.append((R, M))
                
        for M in model.reaction_reactants[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                pairs.append((M, R))
              
    return pairs            

def make_pairs_consuming(model, metabolite, dist, compartment = 0, max_react_size = 0, max_metabolite_occurance = 0, min_flux = 0, ignore_fluxes = 1):
    pairs = []
    tmp_metabolites = set([metabolite])
    reactions = set()
            
    if model.fluxes:
        fluxes = model.fluxes
    
    for d in range(dist):
        metabolites = tmp_metabolites.copy()
        tmp_metabolites.clear()
        for M in metabolites:
            for R in model.metabolite_reactions[M]:
                if not compartment or R in model.compartment_reactions[compartment]:    # compratment nastavljen
                    if not max_react_size or len(model.reaction_metabolites[R]) <= max_react_size:  # max reaction size
                        if not min_flux or abs(fluxes[R]) >= min_flux:    # min flux
                            if (ignore_fluxes or fluxes[R] >= 0) and (R in model.metabolite_as_reactant[M]): 
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati produkte
                                for MM in model.reaction_products[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                                                             
                            if (not ignore_fluxes and fluxes[R] < 0 and R in model.metabolite_as_product[M]):
                                reactions.add(R)                
                                # v nadaljevanju je potrebno gledati reaktante
                                for MM in model.reaction_reactants[R]:
                                    if not max_metabolite_occurance or len(model.metabolite_reactions[MM]) <= max_metabolite_occurance:
                                        tmp_metabolites.add(MM)                                       
    for R in reactions:
        for M in model.reaction_products[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                pairs.append((R, M))
                
        for M in model.reaction_reactants[R]:
            if not max_metabolite_occurance or get_metabolite_occurrence(model, M) <= max_metabolite_occurance:
                pairs.append((M, R))
              
    return pairs        

    
def direct_producers(model, metabolite, fluxes = 0):
    if not fluxes:
        return [R for M in metabolite for R in model.metabolite_as_product[M]]
    else:
        return [R for M in metabolite for R in model.metabolite_as_product[M] if fluxes[R] >= 0] +  [R for M in metabolite for R in model.metabolite_as_reactant[M] if fluxes[R] < 0]

def direct_consumers(model, metabolite, fluxes = 0):
    if not fluxes:
        return [R for M in metabolite for R in model.metabolite_as_reactant[M]]
    else:        
        return [R for M in metabolite for R in model.metabolite_as_reactant[M] if fluxes[R] >= 0] +  [R for M in metabolite for R in model.metabolite_as_product[M] if fluxes[R] < 0]
