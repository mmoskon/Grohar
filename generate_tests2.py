import model_data as md
import pathway_data as pd
import alignment
from model_meta_data import model_meta_data


model = md.model_data()
model.load_model_cobra_sbml("cge_model.xml")
#model.load_model_cobra_sbml("cge00020_model.xml")
model_id = model.model.id                
model_meta_data = model_meta_data(model_id, model = model, from_bigg = 0)
model_meta_data.load_reaction_ECs()


pathway = pd.pathway_data("cge00020")
pathway.load_kgml()


align = alignment.alignment()
align.load_model_pathway(pathway, model, model_meta_data = model_meta_data)


      
distance_to_remove = 3
b, score = align.make_alignment(alignment_type = "Greedy", K_m = 3, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 10, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0, remove_too_distant = 1, remove_too_close = 1, distance_to_remove = distance_to_remove)
b#, score = align.make_alignment(alignment_type = "Greedy", K_m = 3, to_many = 0, from_many = 0, remove_borders = 0, topological_weight = 1, homological_weight_EC = 1, homological_weight_KEGG = 0, homological_additive = 1, remove_too_distant = 1, remove_too_close = 0, distance_to_remove = distance_to_remove)    
align.print_alignment(b)
#align.print_T_sim("A.txt", tsv = 1)
print("*****************")
align.print_T_sim_max("A_max.txt")



#align.T_sim[:,:] = 1
#align.print_T_sim()
#align.adjust_similarities(align.T_sim, 0, 5, 0, 0, remove_too_close=0, remove_too_distant=1, distance_to_remove=distance_to_remove)
#align.print_T_sim()

    


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


"""   
    

def test_model2():
    model = md.model_data()
   
    
    model.create_model("test_model2")
    model.add_metabolite("m1_c", "metabolite1", "c")
    model.add_metabolite("m2_c", "metabolite2", "c")
    model.add_metabolite("m3_c", "metabolite3", "c")
    model.add_metabolite("m4_c", "metabolite4", "c")
    model.add_metabolite("m5_c", "metabolite5", "c")
    model.add_metabolite("m6_c", "metabolite6", "c")
    
    
    model.add_reaction("r_in_1", "r in 1", ["m1_c"], [], [-1000, 1000])
    model.add_reaction("r_in_2", "r in 2", ["m2_c"], [], [-1000, 1000])
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
    b, score = align.make_alignment(alignment_type = "Greedy", K_m = 3, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 10, homological_weight_EC = 0, homological_weight_KEGG = 0, homological_additive = 0, remove_too_distant = 1, remove_too_close = 0, distance_to_remove = distance_to_remove)    
    align.print_alignment(b)
    #print()
    #print(align.KEGG_metabolites_dst)    
    #print()
    #print(align.KEGG_metabolites_src)
    #print()
    #print()
    align.print_T_sim()


    
    align.T_sim[:,:] = 1
    align.print_T_sim()
    align.adjust_similarities(align.T_sim, 0, 5, 0, 0, remove_too_close=0, remove_too_distant=1, distance_to_remove=distance_to_remove)
    align.print_T_sim()

        
    
    print()
    
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
    return align

align = test_model_pathway()
"""

"""
f = open("matrike.txt", "w")
for i in align.src_out:
    for j in i:
        for k in j:
            f.write(str(k) + " ")
        f.write("\n")
    f.write("********\n")
f.close()


f = open("matrike_remove.txt", "w")
for i in align.src_out_remove:
    for j in i:
        for k in j:
            f.write(str(k) + " ")
        f.write("\n")
    f.write("********\n")
f.close()
"""



