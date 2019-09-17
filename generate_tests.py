import model_data as md
import pathway_data as pd
import bigg_data
import alignment
#from paths import paths
import make_pairs
import draw_plotly

import plotly

def generate_EC_file(model, prefix):
    f = open("data/" + prefix + "_EC.txt", 'w')
    EC1 = 1
    EC2 = 1
    EC3 = 1
        
    for r in model.reactions:
        EC4 = r
        f.write(r + ":::" + str(EC1) + "." + str(EC2) + "." + str(EC3) + "." + str(EC4) + "\n")
                
    f.close()

def generate_KEGG_file(model, prefix):
    f = open("data/" + prefix + "_KEGG.txt", 'w')
       
    for m in model.metabolites:
        KEGG_ID = model.metabolites[m]
        f.write(m.split("_")[0] + " " + KEGG_ID + "\n")        
    f.close()
    



def test_model1():
    model = md.model_data()
   
    
    model.create_model("test_model")
    model.add_metabolite("m1_c", "metabolite1", "c")
    model.add_metabolite("m2_c", "metabolite2", "c")
    model.add_metabolite("m3_c", "metabolite3", "c")
    model.add_metabolite("m4_c", "metabolite4", "c")
    model.add_metabolite("m5_c", "metabolite5", "c")
    model.add_metabolite("m6_c", "metabolite6", "c")
    
    
    model.add_reaction("r_in_1", "r in 1", [], ["m1_c"], [-1000, 1000])
    model.add_reaction("r_in_2", "r in 2", [], ["m2_c"], [-1000, 1000])
    model.add_reaction("r1", "reaction 1", ["m1_c", "m2_c"], ["m3_c"], [0, 1000])
    model.add_reaction("r2", "reaction 2", ["m3_c"], ["m4_c"], [-1000, 1000])
    model.add_reaction("r3", "reaction 3", ["m4_c"], ["m5_c", "m6_c"], [0, 1000])
    model.add_reaction("r_out_5", "r out 5", ["m5_c"], [], [0, 1000])
    model.add_reaction("r_out_6", "r out 6", ["m6_c"], [], [0, 1000])
    
    model.set_objective("r_out_5")
    
    #print(model.model.objective)
    #model.run_FBA()
    #print(model.sol)
    
    generate_EC_file(model,model.model.id)
    generate_KEGG_file(model,model.model.id)

    model.save_model_cobra_sbml(model.model.id + ".xml")
    
    return model

def test_model2():
    model = md.model_data()
   
    
    model.create_model("test_model2")
    model.add_metabolite("m1_c", "metabolite1", "c")
    model.add_metabolite("m2_c", "metabolite2", "c")
    model.add_metabolite("m3_c", "metabolite3", "c")
    model.add_metabolite("m4_c", "metabolite4", "c")
    model.add_metabolite("m5_c", "metabolite5", "c")
    model.add_metabolite("m6_c", "metabolite6", "c")
    
    
    model.add_reaction("r_in_1", "r in 1", [], ["m1_c"], [-1000, 1000])
    model.add_reaction("r_in_2", "r in 2", [], ["m2_c"], [-1000, 1000])
    model.add_reaction("r1", "reaction 1", ["m1_c", "m2_c"], ["m3_c"], [-1000, 1000])
    model.add_reaction("r2", "reaction 2", ["m3_c"], ["m4_c"], [-1000, 1000])
    model.add_reaction("r3", "reaction 3", ["m4_c"], ["m5_c", "m6_c"], [0, 1000])
    model.add_reaction("r_out_5", "r out 5", ["m5_c"], [], [0, 1000])
    model.add_reaction("r_out_6", "r out 6", ["m6_c"], [], [0, 1000])
    
    model.set_objective("r_out_5")
    
    #print(model.model.objective)
    #model.run_FBA()
    #print(model.sol)
    
    generate_EC_file(model,model.model.id)
    generate_KEGG_file(model,model.model.id)

    model.save_model_cobra_sbml(model.model.id + ".xml")
    
    return model

def test_model3():
    model = md.model_data()
   
    
    model.create_model("test_model2")
    model.add_metabolite("m1_c", "metabolite1", "c")
    model.add_metabolite("m2_c", "metabolite2", "c")
    model.add_metabolite("m3_c", "metabolite3", "c")
    model.add_metabolite("m4_c", "metabolite4", "c")
    model.add_metabolite("m5_c", "metabolite5", "c")
    model.add_metabolite("m6_c", "metabolite6", "c")
    
    
    model.add_reaction("r_in_1", "r in 1", [], ["m1_c"], [-1000, 1000])
    model.add_reaction("r_in_2", "r in 2", [], ["m2_c"], [-1000, 1000])
    model.add_reaction("r1", "reaction 1", ["m1_c", "m2_c"], ["m3_c"], [-1000, 1000])
    model.add_reaction("r2", "reaction 2", ["m3_c"], ["m4_c"], [-1000, 1000])
    model.add_reaction("r2_x", "reaction 2_x", ["m3_c"], ["m4_c"], [-1000, 1000])
    model.add_reaction("r3", "reaction 3", ["m4_c"], ["m5_c", "m6_c"], [0, 1000])
    model.add_reaction("r_out_5", "r out 5", ["m5_c"], [], [0, 1000])
    model.add_reaction("r_out_6", "r out 6", ["m6_c"], [], [0, 1000])
    
    model.set_objective("r_out_5")
    
    #print(model.model.objective)
    #model.run_FBA()
    #print(model.sol)
    
    generate_EC_file(model,model.model.id)
    generate_KEGG_file(model,model.model.id)

    model.save_model_cobra_sbml(model.model.id + ".xml")
    
    model.fluxes = {}
    for r in model.reactions:
        model.fluxes[r] = 1000
    model.fluxes['r2_x'] = 1000
    

    
    
    return model


def test_pathway1():
    pathway = pd.pathway_data("test_pathway")
    pathway.add_metabolite("m1", "metabolite1")
    pathway.add_metabolite("m2", "metabolite2")
    pathway.add_metabolite("m3", "metabolite3")
    pathway.add_metabolite("m4", "metabolite4")
    pathway.add_metabolite("m5", "metabolite5")
    pathway.add_metabolite("m6", "metabolite6")
    
    pathway.add_reaction("r1", "reaction 1", ["m1", "m2"], ["m3"], 0)
    pathway.add_reaction("r2", "reaction 2", ["m3"], ["m4"], 1)
    pathway.add_reaction("r3", "reaction 3", ["m4"], ["m5", "m6"], 0)
    
    generate_EC_file(pathway, pathway.pathway_id)
    
    
    return pathway    


def test_model_pathway():
    m = test_model1()
    bd = bigg_data.bigg_data(m.model.id)
    bd.load_metabolite_kegg_ids()
    bd.load_reaction_ECs()
    
    
    p = test_pathway1()
    p.load_ECs("test_pathway")
    
    align = alignment.alignment()
    align.load_model_pathway(p, m, bd)#, max_metabolite_occurance, compartment = "c")
    #align.load_pathway_pathway(p, p)
        
    distance_to_remove = 3
    #b, score = align.make_alignment(alignment_type = "Greedy", K_m = 3, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 0, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0, remove_too_distant = 1, remove_too_close = 1, distance_to_remove = distance_to_remove)
    b, score = align.make_alignment(alignment_type = "Greedy", K_m = 3, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 10, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0, remove_too_distant = 1, remove_too_close = 0, distance_to_remove = distance_to_remove)    
    align.print_alignment(b)
    #print()
    #print(align.KEGG_metabolites_dst)    
    #print()
    #print(align.KEGG_metabolites_src)
    #print()
    #print()
    align.print_T_sim()
    
    c = align.get_anchors()
    align.print_alignment(c)


    
    #align.T_sim[:,:] = 1
    #align.print_T_sim()
    #align.adjust_similarities(align.T_sim, 0, 5, 0, 0, remove_too_close=0, remove_too_distant=1, distance_to_remove=distance_to_remove)
    #align.print_T_sim()

        
    
    print()

    #paths.get_reaction_connectivity(m)
    #paths.shortest_path(m, "r1", "r_out_5", compartment = "c")
    
    #align.print_matrix("dst", "in")
    #align.print_matrix("dst", "out")
    
    #align.print_matrix("src", "in")
    #align.print_matrix("src", "out")
    
    #print(align.reactions_dst)
    #print()
    #for i in align.M_dst:
    #    for j in i:
    #        print(j, end= " ")
    #    print()
    return (align, m)
"""    
def test_model_model():
    m_dst = test_model1()
    bd_dst = bigg_data.bigg_data(m_dst.model.id)
    bd_dst.load_metabolite_kegg_ids()
    bd_dst.load_reaction_ECs()
    
    m_src = test_model2()
    bd_src = bigg_data.bigg_data(m_src.model.id)
    bd_src.load_metabolite_kegg_ids()
    bd_src.load_reaction_ECs()
    
    align = alignment.alignment()   
    align.load_model_model(m_dst, m_src, bd_dst, bd_src, max_metabolite_occurence = 0, compartment_dst = 'c', compartment_src ='c', reversibility = 1)
    b, score = align.make_alignment(alignment_type = "Greedy", K_m = 10, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 0, homological_weight_EC = 10, homological_weight_KEGG = 0, homological_additive = 0, remove_too_distant = 1, remove_too_close = 1, distance_to_remove = 3)
    
    align.print_alignment(b)
    align.print_T_sim()
    return align
"""

a, model = test_model_pathway()

model.run_FBA()
reaction_types = 0
met_names = {'m1_c', 'm2_c', 'm6_c'}
distance = 2
compartments = 'c'
max_react_size = 100
max_metabolite_occurance = 100
min_flux = 0



(pairs, fluxes) = make_pairs.make_pairs(reaction_types, model, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux = min_flux)

draw_plotly.draw_plotly(model, met_names, pairs, fluxes)
#
# 
#




