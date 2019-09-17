import model_data as md
#import pathway_data as pd
#import bigg_data
#import alignment
#from paths import paths
#import make_pairs
#import draw_plotly

#import plotly

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

 