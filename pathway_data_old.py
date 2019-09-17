import xml.etree.cElementTree as ET
from collections import defaultdict    

# https://pythonhosted.org/bioservices/kegg_tutorial.html
# https://pythonhosted.org/bioservices/references.html#bioservices.kegg.KEGG.parse_kgml_pathway
from bioservices.kegg import KEGG
#from bioservices.kegg import KEGGParser


class pathway_data:
    def __init__(self, pathway_id):
        self.pathway_id = pathway_id
        
        
        # kljuc = ID_reakcije, vsebina = mnozica reakcij za id-jem (lahko jih je vec)
        self.reactions = defaultdict(set)
        
        # kljuc = ime reackije, vsebina = ID reackije
        self.reaction_ids = {}
        

        # kljuc = ID_metabolita, vsebina = ime metabolita
        self.metabolites = {}
        
         # kljuc = ime metabolita, vsebina = ID metabolita
        self.metabolite_ids = {}
        
        
        # kljuc = ID_gena, vsebina = ime gena
        #self.genes = {}
        
        
        # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki sodelujejo v reakciji
        self.reaction_metabolites = defaultdict(set)
        
        # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki jih reakcija producira
        self.reaction_products = defaultdict(set)
        
        # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki jih reakcija porablja
        self.reaction_reactants = defaultdict(set)
       
        # kljuc = ID_reakcije, vsebina = mnozica genov, ki reakcijo katalizirajo
        #self.reaction_genes = defaultdict(set)
        
        # kljuc = ID_reakcije, vsebina = mnozica EC-jev, ki reakcijo katalizirajo
        self.reaction_ECs = defaultdict(set)
        
        # kljuc je ime gena, vrednost je EC number
        #self.gene_EC = {}     
        
        # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih metabolit sodeluje
        #self.metabolite_reactions = defaultdict(set)
        
        # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih je metabolit reaktant
        #self.metabolite_as_reactant = defaultdict(set)
        
        # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih je metabolit produkt
        #self.metabolite_as_product = defaultdict(set)

        # kljuc = ID_reakcije, vsebina = reversibilnost reakcije
        self.reversibility_reactions = {}

        
        # ID reakcij urejen po indeksu 
        self.listed_reactions = [] 
             
        # ID metabolitov urejen po indeks
        self.listed_metabolites = [] 


    def add_metabolite(self, m_id, m_name):
        if m_id not in self.metabolites:
            self.metabolites[m_id] = m_name
            self.metabolite_ids[m_name] = m_id
            self.listed_metabolites.append(m_id)
            
    def add_reaction(self, r_id, r_name, reactants, products, reversibility = 0):
        if r_id not in self.reactions:
            self.reactions[r_id] = r_name
            self.reaction_ids[r_name] = r_id
            self.listed_reactions.append(r_id)
                                    
            for m_r in reactants:
                if m_r in self.metabolites:
                    self.reaction_metabolites[r_id].add(m_r)
                    self.reaction_reactants[r_id].add(m_r)
                 
            for m_p in products:
                if m_p in self.metabolites:
                    self.reaction_metabolites[r_id].add(m_p)
                    self.reaction_products[r_id].add(m_p)

            self.reversibility_reactions[r_id] = reversibility
                    
            
    def load_kgml(self, file_name ="", ec_file=""):
        if not file_name:
            file_name = self.pathway_id + ".xml"
        
        try:          
            f = open(file_name, 'r')        
            self.kgml = f.read()
            if not ec_file:
                ec_file = self.pathway_id + "_EC.txt"
            self.parse_kgml(ec_file)
        except:
            self.get_kegg()
        
        
        
    def get_kegg(self):
        #try:
            self.kegg = KEGG()
            kegg = self.kegg
                   
            #self.kgml = kegg.parse(kegg.get(pathway_id))
            #self.pathway = kegg.parse_kgml_pathway(pathway_id) 
            self.kgml = kegg.get(self.pathway_id, "kgml")
            self.parse_kgml()
            
            self.save_kegg()
            self.save_ECs()
            
            
       # except:
       #     print("not found")
       
    def save_kegg(self, file_name = ""):
        if not file_name:
            file_name = self.pathway_id + ".xml"
        
        try:
            f = open(file_name, "w")
            f.write(self.kgml)
            f.close()
        except:
            print("error")
                
    def parse_kgml(self, ec_file = ""):
        # http://biopython.org/DIST/docs/api/Bio.KEGG.KGML.KGML_parser-pysrc.html
        # https://github.com/deep-introspection/kegg-kgml-parser-python/blob/master/keggparser/parse_KGML.py
            
            
            
        tree = ET.fromstring(self.kgml)
            
                    
        for reaction in tree.getiterator('reaction'):
            r_id = reaction.get('id')
            r_name = reaction.get('name') # lahko je vec imen locenih s presledki
            r_names= set(reaction.get('name').split()) # mnozica imen
            
            self.reactions[r_id] = r_names
            self.reaction_ids[r_name] = r_id                   
            self.listed_reactions.append(r_id)                              
            
                          
            for sub in reaction.getiterator('substrate'):
                self.reaction_metabolites[r_id].add(sub.get('id'))
                self.reaction_reactants[r_id].add(sub.get('id'))
                #substrates.append(sub.get('name'))
            
            for prod in reaction.getiterator('product'):
                self.reaction_metabolites[r_id].add(prod.get('id'))
                self.reaction_products[r_id].add(prod.get('id'))
                #products.append(prod.get('name'))
                
            self.reversibility_reactions[r_id] = 1 if reaction.get('type') == 'reversible' else 0
            #reactions[i] = {'reaction': reaction, 'substrates': substrates, 'products': products, 'gene':[], 'reversible': reversible}

        EC_file_loaded = False
        if ec_file:
            try:
                self.load_ECs(ec_file)
                EC_file_loaded = True
            except:
                self.kegg = KEGG()
                    
        
        for entry in tree.getiterator('entry'):
            if not EC_file_loaded:                
                if entry.get('type') == 'gene': # or entry.get('type') == 'ortholog':
                    genes = entry.get('name').split()
                    gene_reaction_name = entry.get('reaction')
                    #print(gene_reaction_name)
                    gene_reaction_id = self.reaction_ids[gene_reaction_name]
                    for g in genes:
                        #self.reaction_genes[gene_reaction_id].add(g) 
                        EC = self.get_EC(g)
                        #self.gene_EC[g] = EC
                        for e in EC:
                            self.reaction_ECs[gene_reaction_id].add(e)
                
                
            if entry.get('type') == 'compound':
                metabolite = entry.get('name')
                metabolite_id = entry.get('id')
                self.metabolites[metabolite_id] = metabolite
                self.metabolite_ids[metabolite] = metabolite_id
                self.listed_metabolites.append(metabolite_id)
                  
        
        
    def get_EC(self, gene_id):
        # gene = kegg.parse(kegg.get("cge:100754434"))
        kegg = self.kegg
        try:
            gene = kegg.parse(kegg.get(gene_id))
            for i in gene['ORTHOLOGY'].values():
                EC = i.split('[')[1].split(']')[0].replace("EC:","")
            return set(EC.split())
        except:
            return set()
   

    def save_ECs(self, file_name = ""):
        
        file_name = self.pathway_id + "_EC.txt"
        try:
            f = open(file_name, "w")
            for re in self.reaction_ECs:
                f.write(re + ":")
                for ec in self.reaction_ECs[re]:
                    f.write(ec+" ")
                f.write("\n")
            f.close()
        except:
            pass
        

    def load_ECs(self, file_name = ""):
        file_name = self.pathway_id + "_EC.txt"
        self.reaction_ECs = defaultdict(set)
        f = open(file_name, "r")
        for l in f:
            l = l.split(":")
            reaction_id = l[0]
            ECs = set(l[1].strip().split(" "))
            ECs.discard('')
            self.reaction_ECs[reaction_id] = ECs                            
        f.close()
        