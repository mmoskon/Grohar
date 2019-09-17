import pathway_data as pd
from bioservices.kegg import KEGG 


organism = "cge"

kegg = KEGG()

pathway_multi = pd.pathway_data(organism + "_CH")


#kegg.organism = organism

"""
pathways = kegg.pathwayIds


for p in pathways:
    try:
        entry = p.replace("path:","")
        print("loading pathway", entry)
        pathway_multi.load_kgml(entry)
    except Exception as inst:
        print(inst)
"""
# Carbohydrate metabolism

pathways = { "00010": "Glycolysis / Gluconeogenesis", 
      "00020" : "Citrate cycle (TCA cycle)",
      "00030" : "Pentose phosphate pathway",
      "00040" : "Pentose and glucuronate interconversions",
      "00051" : "Fructose and mannose metabolism",
      "00052" : "Galactose metabolism",
      "00053" : "Ascorbate and aldarate metabolism",
      "00500" : "Starch and sucrose metabolism",
      "00520" : "Amino sugar and nucleotide sugar metabolism",
      "00620": "Pyruvate metabolism",
      "00630" : "Glyoxylate and dicarboxylate metabolism",
      "00640" : "Propanoate metabolism",
      "00650" : "Butanoate metabolism",
      "00660" : "C5-Branched dibasic acid metabolism",
      "00562" : "Inositol phosphate metabolism"}

for p in pathways:
    try:
        entry = organism + p
        print("loading pathway", entry)
        pathway_multi.load_kgml(entry)
    except Exception as inst:
        print(inst)



#pathway_multi = pd.pathway_data("cge00020")
#pathway_multi.load_kgml("cge00020")
    
pathway_multi.save_to_sbml()
