import bigg_data
import model_data as md

model = md.model_data()
model.load_model_cobra_sbml("iCHOv1_bigg.xml")

bd = bigg_data(model)

#bd.get_metabolite_kegg_ids()
#bd.save_metabolite_kegg_ids("iCHOv1_meta_keg.txt")
bd.load_metabolite_kegg_ids("iCHOv1_meta_keg.txt")

#bd.get_reaction_ECs()
bd.save_reaction_ECs("iCHOv1_EC.txt")
bd.load_reaction_ECs("iCHOv1_EC.txt")