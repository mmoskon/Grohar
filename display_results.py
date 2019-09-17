from bioservices.kegg import KEGG
from model_data import model_data

def load_alignment():
    pass

def print_alignment_kegg(model):
    f = open("cor.txt")
    f_o = open("cor_readable.txt", "w")
    kegg = KEGG()
    
    for i in f:
        if ":***:" in i:
            k, b = i.split(":***:")
            b = b.strip()
            
            if not k == "MULTIR":
                k = kegg.get(k)
            
                i1 = k.find("NAME") + 4
                i2 = k[i1:].find("\n")
                    
                k = k[i1:i1+i2].strip()             
            
            if not b == "MULTIR":
                b = model.reactions[b]
            
            print(k, ":***:", b)
            f_o.write(k + ":***:" + b + "\n")
    f.close()
    f_o.close()

    

model = model_data()
model.load_model_cobra_sbml("e_coli_core.xml")
print_alignment_kegg(model)