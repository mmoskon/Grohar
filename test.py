import model_data as md
import pathway_data as pd
from alignment import alignment
from bigg_data import bigg_data

"""
model = md.model_data()
#model.load_model_sbml("iCHOv1_final.xml")
#model.load_model_cobra_sbml("iCHOv1_final.xml")
#model.load_model_mat("gimmeDG44.mat")
model.load_model_mat("gimmeS.mat")
#model.load_model_mat("gimmeK1.mat")

"""
align = alignment()
dummy = -1
reload = 0

if dummy == 1:
    align.load_dummy()
elif dummy == -1:
    model = md.model_data()
    #model.load_model_mat("gimmeS.mat") # omogociti uporabniku, da izbere tip modela iz seznama modelov, ki ga pridobimo iz bigg-a
    #model.load_model_cobra_sbml("iCHOv1_bigg.xml")
    model.load_model_cobra_sbml("e_coli_core.xml")
    
    
    #if reload:
    #    model.get_ECs()
    #    model.save_ECs("ec.txt")
    #else:    
    #    model.get_ECs_from_file("ec.txt")
    
    #pathway.load_kegg("cge00250")
    
    #pathway = pd.pathway_data("cge00250")
    pathway = pd.pathway_data("eco00020")
    pathway.load_kgml()
    
    
    max_metabolite_occurance = 10
    
    #bigg_id = "iCHOv1"
    bigg_id = model.model.id
    print("*******************")
    print(bigg_id)
    bd = bigg_data(bigg_id)
    bd.load_metabolite_kegg_ids()
    bd.load_reaction_ECs()
    
    align.load_model_pathway(pathway, model, bd, max_metabolite_occurance, compartment = "c")
        

    
else:
    pathway1 = pd.pathway_data("cge00250")
    pathway1.load_kgml()
    
    pathway2 = pd.pathway_data("hsa00250")
    pathway2.load_kgml()
    
    align.load_pathway_pathway(pathway1, pathway2)
    #<pathway.save_ECs("cge00250_ec.txt")

#align.load_model_pathway(pathway, model)

#print(align.align_probs(1))
a, score = align.make_alignment(3, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0)
align.print_alignment(a)
print(align.T_sim.max())
print(score)
b, score = align.make_alignment(3, to_many = 0, from_many = 0, topological = 0, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0)
align.print_alignment(b)
print(align.T_sim.max())
print(score)


align.print_alignment_kegg(b)




#model.save_pajek("test", 1)