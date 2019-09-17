#http://bigg.ucsd.edu/data_access

import requests
from time import sleep
from collections import defaultdict    
from bioservices.kegg import KEGG

class model_meta_data:

    def __init__(self, model_id, model = None, from_bigg = 1):    
        self.meta_kegg_ids =  {}
        self.reaction_ECs = defaultdict(set)
        self.model_id = model_id
        #self.model_id = self.model.bigg_id
        self.from_bigg = from_bigg
        self.model = model
    

    def get_metabolite_kegg_ids_from_bigg(self):
        self.meta_kegg_ids =  {}
        queried = []
                
        try:
            metabolites = requests.get('http://bigg.ucsd.edu/api/v2/models/' + self.model_id + '/metabolites').json()
        except Exception as inst:
            print(inst)
        else:
            for m in metabolites['results']:
                m_id = m['bigg_id'] + '_' + m['compartment_bigg_id']
                m_id_short = m['bigg_id']
                if m_id_short not in queried:
                    queried.append(m_id_short)
                    
                    # we are respectful and do not send more than 10 requests per second
                    sleep(0.1)
                    try:
                        metabolite = requests.get('http://bigg.ucsd.edu/api/v2/models/' + self.model_id + '/metabolites/' + m_id).json()
                        
                        print("metabolite:",metabolite['bigg_id'])
                        if 'database_links' in metabolite and 'KEGG Compound' in metabolite['database_links']:
                            kegg_id = metabolite['database_links']['KEGG Compound'][0]['id']
                            m_id_short = m_id_short.replace("__", "_")
                            self.meta_kegg_ids[m_id_short] = kegg_id
                    except Exception as inst:
                        print(inst)
                    
                    
                                
    def save_metabolite_kegg_ids(self, file_name = ""):
        if not file_name:
            file_name = "data/" + self.model_id + "_KEGG.txt"
        
        f = open(file_name, 'w')
        for i in self.meta_kegg_ids:
            f.write(i + " " + self.meta_kegg_ids[i] + "\n")
        f.close()
        
    def load_metabolite_kegg_ids(self, file_name = ""):
        if not file_name:
            file_name =  "data/" + self.model_id + "_KEGG.txt"
                        
        try:
            f = open(file_name, 'r')
            self.meta_kegg_ids = {}
            
            for i in f:
                m_id, kegg_id = i.strip().split(" ")
                self.meta_kegg_ids[m_id] = kegg_id            
            f.close()
            
            print("Metabolite ID data loaded from a file")
        except:
            self.get_metabolite_kegg_ids_from_bigg()
            self.save_metabolite_kegg_ids("")
            print("Metabolite ID data loaded from BiGG")

    def get_reactions_ECs_from_bigg(self):
        self.reaction_ECs = defaultdict(set)
   
        try:
            reactions = requests.get('http://bigg.ucsd.edu/api/v2/models/' + self.model_id + '/reactions').json()
        except Exception as inst:
            print(inst)
        else:
            for r in reactions['results']:
                r_id = r['bigg_id']
                
                # we are respectful and do not send more than 10 requests per second
                sleep(0.1)
                try:
                    reaction = requests.get('http://bigg.ucsd.edu/api/v2/models/' + self.model_id + '/reactions/' + r_id).json()
                    print("reaction",reaction['bigg_id'])
                    if 'database_links' in reaction and 'EC Number' in reaction['database_links']:
                            ECs = reaction['database_links']['EC Number']
                            for e in ECs:
                                self.reaction_ECs[r_id].add(e['id'])              
                except Exception as inst:
                    print(inst)
            print("EC data loaded from BiGG")
        
    def get_reaction_ECs_from_kegg(self):
        self.reaction_ECs = defaultdict(set)
        
        kegg = KEGG()
        for r in self.model.reactions:
            ECs = []
            try:
                reacts = r.split(" ")
                for i in reacts:
                    if i not in self.reaction_ECs:
                        print("KEGG reaction", i)
                        ECs += kegg.parse(kegg.get(i))['ENZYME']                 
                        for e in ECs:
                            self.reaction_ECs[i].add(e)      
                
            except Exception as inst:
                print(inst)                
            #for e in ECs:
            #    self.reaction_ECs[r].add(e)      
       
        print("EC data loaded from KEGG")
            
    def get_reaction_EC_from_bigg(self, r_id):
        try:
            reaction = requests.get('http://bigg.ucsd.edu/api/v2/models/' + self.model_id + '/reactions/' + r_id).json()
            print(reaction['bigg_id'])
            if 'database_links' in reaction and 'EC Number' in reaction['database_links']:
                if r_id in self.reaction_ECs:
                    del self.reaction_ECs[r_id]
                ECs = reaction['database_links']['EC Number']
                for e in ECs:
                    self.reaction_ECs[r_id].add(e['id'])              
        except Exception as inst:
            print(inst)
            
    def save_reaction_ECs(self, file_name = ""):
        if not file_name:
            file_name =  "data/" + self.model_id + "_EC.txt"
        
        f = open(file_name, 'w')
        for i in self.reaction_ECs:
            f.write(i + ":::")
            for j in self.reaction_ECs[i]:
                f.write(j + " ")
            f.write("\n")
        f.close()

    def load_reaction_ECs(self, file_name = ""):
       
        if not file_name:
            file_name =  "data/" + self.model_id + "_EC.txt"
       
        try:
            f = open(file_name, 'r')
            self.reaction_ECs = defaultdict(set)
            for i in f:
                r_id, ECs = i.strip().split(":::")  
                ECs = set(ECs.split(" "))
                self.reaction_ECs[r_id] = ECs            
            f.close()
            print("EC data loaded from a file")
        except:
            
            if self.from_bigg:
                self.get_reactions_ECs_from_bigg()
            else:                
                self.get_reaction_ECs_from_kegg()
            self.save_reaction_ECs("")