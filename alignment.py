import numpy as np
from collections import defaultdict    
import operator

from bioservices.kegg import KEGG

class alignment:
   
    def __init__(self):
        self.K_m = 0
        self.distance_to_remove = 0
        self.last_alignment = None
          
    def load_model(self, model, model_meta_data = None, max_metabolite_occurence = 0, compartment = "", reversibility = 1, load_kegg = 1):
        
        if compartment:
            reactions = list(model.compartment_reactions[compartment])            
        else:
            reactions = list(model.reactions)
                    
     
        r = len(reactions)
        M = np.zeros((r, r))
        
        
        #print(model.bigg_id)
        #model.bd.load_reaction_ECs()
        #model.bd.load_metabolite_kegg_ids()
        
        
        EC = []
        
        KEGG_metabolites = []
      
        reaction_names = []
        metabolites_to_remove = set()
        if max_metabolite_occurence:
            for m in model.metabolites: #model.metabolite_reactions:
                # ven pomecemo prepogoste metabolite!!!
                if len(model.metabolite_reactions[m]) > max_metabolite_occurence:
                    metabolites_to_remove.add(m)
        """
        f = open("debug_KEGG.txt","w")
        for m in model.metabolites:
            if "[" in m:
                 m_id = m.split("[")[0]
            else:
                 m_id = m[:-(1 + len(m.split("_")[-1]))]
            
            m_id = m_id.replace("__","_")
            if m_id in model.bd.meta_kegg_ids:
                f.write(m_id + ": " + model.bd.meta_kegg_ids[m_id]+"\n")
            else:
                f.write(m_id + ": !!!"+"\n")
        f.close()
        
           
        f = open("debug_EC.txt","w")
        for r in reactions:        
            f.write(r + ": " + " ".join(s for s in model.bd.reaction_ECs[r])+"\n")
        f.close()
        """
        
        
        for (i1, r1) in enumerate(reactions):
            reaction_names.append(model.reactions[r1])
            """
            if r1 in model.reaction_genes:
                EC.append(model.reaction_genes[r1])
            else:            
                EC.append(set())
            """    
            
            if load_kegg:
                EC_set = set()
                if model_meta_data:                     
                    for rr in r1.split(" "):
                        if rr in model_meta_data.reaction_ECs:
                            EC_set |= set(model_meta_data.reaction_ECs[rr])
                EC.append(EC_set)
            else:
                EC.append({r1})
                    
            for (i2, r2) in enumerate(reactions):
                if reversibility:
                    lb1, ub1 = model.reaction_bounds[r1]
                    lb2, ub2 = model.reaction_bounds[r2]
                    
                                    
                    if ub1 > 0:
                        if ub2 > 0:
                            if (model.reaction_products[r1] & model.reaction_reactants[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1
                               continue
                               #print("0: ",r1, r2)
                        if lb2 < 0 and r1 != r2:
                            if (model.reaction_products[r1] & model.reaction_products[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1                           
                               continue
                               #print("1: ",r1,r2)
                    
                    if lb1 < 0:
                        if ub2 > 0 and r1 != r2:
                            if (model.reaction_reactants[r1] & model.reaction_reactants[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1
                               continue
                               #print("3: ",r1, r2)
                        if lb2 < 0:
                            if (model.reaction_reactants[r1] & model.reaction_products[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                               M[i1,i2] = 1       
                               continue
                               #print("4: ",r1, r2)
                else:
                    if (model.reaction_products[r1] & model.reaction_reactants[r2]) - metabolites_to_remove:  # vsaj en produkt prve je reaktant druge
                       M[i1,i2] = 1
                      
            KEGG_reactants = set()
            KEGG_products = set()
            
            for m in model.reaction_reactants[r1]:
                # remove the compartment information!
                             
                if "[" in m:
                     m_id = m.split("[")[0]
                elif "_" in m:
                     m_id = m[:-(1 + len(m.split("_")[-1]))]
                else:
                    m_id = m
                
                m_id = m_id.replace("__","_")
                m_id = m_id.replace("cpd:","")
              
                
                
                if load_kegg and model_meta_data and m_id in model_meta_data.meta_kegg_ids:
                    KEGG_reactants.add(model_meta_data.meta_kegg_ids[m_id])
                else:
                    KEGG_reactants.add(m_id)
                                         
                
            for m in model.reaction_products[r1]:
                # remove the compartment information!
                if "[" in m:
                     m_id = m.split("[")[0]
                elif "_" in m:
                     m_id = m[:-(1 + len(m.split("_")[-1]))]
                else:
                    m_id = m
                     
                m_id = m_id.replace("__","_")
                m_id = m_id.replace("cpd:","")
                
            
                if load_kegg and model_meta_data and m_id in model_meta_data.meta_kegg_ids:
                    KEGG_products.add(model_meta_data.meta_kegg_ids[m_id])
                else:
                    KEGG_products.add(m_id)

              

            KEGG_metabolites.append([KEGG_reactants, KEGG_products])            
            
        
        return (reactions, reaction_names, M, EC, KEGG_metabolites)
    
    def load_pathway(self, pathway, reversibility = 1):
        
        reactions = pathway.listed_reactions
        r = len(reactions)
        M = np.zeros((r,r))
        EC = []
        KEGG_metabolites = []
        reaction_names = []
       
        for (i1, r1) in enumerate(reactions):
            
            #print(r1, pathway.reversibility_reactions[r1])
            
            reaction_names.append(pathway.reactions[r1])
            EC_set = set()
            for rr in  r1.split(" "): # multiple reactions in one
                if rr in pathway.reaction_ECs:
                    EC_set |= pathway.reaction_ECs[rr]
            
            #print(EC_set)
            EC.append(EC_set)
          
            for (i2, r2) in enumerate(reactions):
                if pathway.reaction_products[r1] & pathway.reaction_reactants[r2]:  # vsaj en produkt prve je reaktant druge
                    M[i1,i2] = 1
                    continue
                if reversibility:
                    if r1 != r2:
                        if pathway.reversibility_reactions[r2] and (pathway.reaction_products[r1] & pathway.reaction_products[r2]):  # vsaj en produkt prve je reaktant druge
                            M[i1,i2] = 1
                            continue
                        if pathway.reversibility_reactions[r1] and (pathway.reaction_reactants[r1] & pathway.reaction_reactants[r2]):  # vsaj en produkt prve je reaktant druge
                            M[i1,i2] = 1
                            continue
                    if pathway.reversibility_reactions[r1] and pathway.reversibility_reactions[r2] and (pathway.reaction_reactants[r1] & pathway.reaction_products[r2]):
                        M[i1,i2] = 1
                        continue
                
                    
                     #if pathway.reversibility_reactions[r1]:
                     #   M[i2,i1] = 1
            
            KEGG_reactants = set()
            KEGG_products = set()
          
         
            for m in pathway.reaction_reactants[r1]:
                KEGG_reactants.add(pathway.metabolites[m].replace("cpd:",""))
          
            
            for m in pathway.reaction_products[r1]:
                KEGG_products.add(pathway.metabolites[m].replace("cpd:",""))
            
            KEGG_metabolites.append([KEGG_reactants, KEGG_products])            
            
        return (reactions, reaction_names, M, EC, KEGG_metabolites)                 
    
    def load_model_model(self, model_dst, model_src, model_meta_data_dst = None, model_meta_data_src = None, max_metabolite_occurence = 0, compartment_dst = 'c', compartment_src ='c', reversibility = 1):
        print("loading model dst")
        (self.reactions_dst, self.reaction_names_dst, self.M_dst, self.EC_dst, self.KEGG_metabolites_dst) = self.load_model(model_dst, model_meta_data_dst, max_metabolite_occurence, compartment_dst, reversibility, load_kegg = 0)                
        
        print("loading model src")
        (self.reactions_src, self.reaction_names_src, self.M_src, self.EC_src, self.KEGG_metabolites_src) = self.load_model(model_src, model_meta_data_src, max_metabolite_occurence, compartment_src, reversibility, load_kegg = 0)                
        
        #self.pathway = pathway
        #self.model = model
        
                 
    def load_model_pathway(self, pathway, model, model_meta_data = None, max_metabolite_occurence = 0, compartment = 'c', reversibility = 1):
        print("loading pathway")
        (self.reactions_src, self.reaction_names_src, self.M_src, self.EC_src, self.KEGG_metabolites_src) = self.load_pathway(pathway, reversibility)
        #print(self.reaction_names_src)
        print("loading model")
        (self.reactions_dst, self.reaction_names_dst, self.M_dst, self.EC_dst, self.KEGG_metabolites_dst) = self.load_model(model, model_meta_data, max_metabolite_occurence, compartment, reversibility)        
        
        #self.pathway = pathway
        #self.model = model
        
    def load_pathway_pathway(self, pathway1, pathway2, reversibility = 1):
        (self.reactions_src, self.reaction_names_src, self.M_src, self.EC_src, self.KEGG_metabolites_src) = self.load_pathway(pathway1, reversibility)
        (self.reactions_dst, self.reaction_names_dst, self.M_dst, self.EC_dst, self.KEGG_metabolites_dst) = self.load_pathway(pathway2, reversibility)
        
        #self.pathway1 = pathway1
        #self.pathway2 = pathway2
    
     
    def print_alignment(self, a, file_name = ""):
        tmp = ""
        
        for i in a:
            tmp += self.reaction_names_src[i] + " : " +  self.reaction_names_dst[a[i]] +  " - score = " + str(self.T_sim[i,a[i]]) + "\n"
        
        tmp += str(len(a)) + " reactions out of " + str(len(self.reactions_src)) + " aligned"
    
        if file_name:
            f = open(file_name, "w")
            f.write(tmp)
            f.close()
        else:
            print(tmp)
    
    def print_alignment_full(self, a, file_name = ""):
        kegg = KEGG()

        tmp = ""


        for i in a:
            name = ""
            for j in self.reaction_names_src[i].split(" "):
                if name:
                    name += " ;"
                kegg_id = j
                
                #print(kegg_id)
                
                r = kegg.get(kegg_id)
                i1 =  r.find("NAME") + 4
                i2 = r[i1:].find("\n")
                
                name += r[i1:i1+i2].strip()             
        
            tmp += str(kegg_id) + " :---: " + str(name) + " :***: " + str(self.reactions_dst[a[i]]) + " :---: " + str(self.reaction_names_dst[a[i]]) + "\n"
            
            """
                TODO: 
                    dodaj izpis reaktantov in produktov iz KEGG-a in iz BiGG-a
                    primerjaj!
           """
        
        if file_name:
            f = open(file_name, "w")
            f.write(tmp)
            f.close()
        else:
            print(tmp)

        
         
    def load_dummy(self):
        self.M_src = np.array(([0,0,1,0],[0,0,1,0],[0,0,0,1],[1,1,0,0]))
        self.M_dst = np.array(([0,1,0],[0,0,1],[1,0,0]))
          


    def generate_matrices(self, X, K_m):
        n = len(X)
      
        M = np.zeros([K_m, n, n])
        M[0] = X
        
        
        
        
        # matrični produkt: np.dot(M_p, M_p)
        for k in range(1, K_m):
            M[k] = np.array(np.array(np.dot(M[k-1], M[0]), dtype = bool), dtype = int)
                    
        M_degrees = np.zeros([K_m, n])
             
        for i in range(K_m):
            M_degrees[i,:] = M[i].sum(axis=1)
         
        neighbours_degrees = np.zeros([K_m, n])
    
            
        for k in range(K_m):
            for i1 in range(n):
                for j1 in range(n):
                    if M[k,i1,j1] == 1 and i1 != j1:
                        neighbours_degrees[k,i1] += M_degrees[0,j1]                                   
        
        num_k_neighbours = np.zeros([n])
        num_k_edges = np.zeros([n])
        for i in range(n):
            #num_k_neighbours[i] = np.sum(np.logical_or(M[k,i,:] for k in range(K_m)))
            num_k_neighbours[i] = np.sum(np.array(M[:,i,:].sum(axis=0), dtype = bool)) # stevilo razlicnih sosedov v sosescini razdalje k
            num_k_edges[i] = M[:,i,:].sum() # stevilo povezav v sosescini razdalje k
    
        #print(neighbours_degrees)
               
        return (M, neighbours_degrees, num_k_neighbours, num_k_edges)
    

    def generate_matrices_borders(self, X, K_m):
        n = len(X)
      
        M = np.zeros([K_m, n, n])
        M[0] = X
        
        # matrični produkt: np.dot(M_p, M_p)
        for k in range(1, K_m):
            M[k] = np.array(np.array(np.dot(M[k-1], M[0]), dtype = bool), dtype = int)
                    
            
        # stopnje vozlišč glede na širino okolice
        #               k (širina okolice) \ reaction_idx:   0    1    3 ....
        #                   0
        #                   1
        # M_degrees =       2
        #                   3
        #                   ...
        #
        # v nadaljevanju uporabljamo le M_degrees[0,:]
        #
        M_degrees = np.zeros([K_m, n])
             
        for i in range(K_m):
            M_degrees[i,:] = M[i].sum(axis=1)
        #print(M)
        #print(M_degrees)
         
        
        # stopnje vozlišl sosedov, glede na širino kolice
        
        
        neighbours_degrees = np.zeros([K_m, n])
    
            
        for k in range(K_m):
            for i1 in range(n):                
                for j1 in range(n):
                    if M[k,i1,j1] == 1 and i1 != j1: # if neighbours?
                        if M_degrees[0,j1]: # if neighbour is not on the border
                            neighbours_degrees[k,i1] += M_degrees[0,j1]                                   
                        else: # M_degrees[0,j1] == 0?
                            neighbours_degrees[k,i1] = -1   # if on the border of the pathway it should not be taken into consideration!!!
                            break
                if neighbours_degrees[k,i1] == 0:
                    neighbours_degrees[k,i1] = -1
        
        num_k_neighbours = np.zeros([n]) # stevilo razlicnih sosedov v sosescini razdalje k
        num_k_edges = np.zeros([n]) # stevilo povezav v sosescini razdalje k
        for i in range(n):            
            num_k_neighbours[i] = np.sum(np.array(M[:,i,:].sum(axis=0), dtype = bool)) # stevilo razlicnih sosedov v sosescini razdalje k
            num_k_edges[i] = M[:,i,:].sum() # stevilo povezav v sosescini razdalje k
    
    
    
        # if on the border???
        for i in range(n):
            for k in range(K_m):
                if M[k,i,:].sum() == 0: # I am not connected to anyone!!!
                    num_k_neighbours[i] = -1    # disregard the number of different neighbours
                    num_k_edges[i] = -1         # disregard the number of edges
                    break
                                   

               
        return (M, neighbours_degrees, num_k_neighbours, num_k_edges)



      
    def topological_prepare(self, K_m):
        
        (self.src, self.neighbours_degrees_src, self.num_k_neighbours_src, self.num_k_edges_src) = self.generate_matrices(self.M_src, K_m)
        (self.dst, self.neighbours_degrees_dst, self.num_k_neighbours_dst, self.num_k_edges_dst) = self.generate_matrices(self.M_dst, K_m)
        
        self.num_vertices_src = len(self.src[0])
        self.num_vertices_dst = len(self.dst[0])
        

    def topological_prepare_both(self, K_m, remove_borders):
        
         
        # if remove borders is set ignore the properties of reactions that are on the borders of the network (input connectivity or output connectivity is zero!)
        # remove borders only for the pathway!!!
        
        if remove_borders:
            (self.src_out, self.neighbours_degrees_src_out, self.num_k_neighbours_src_out, self.num_k_edges_src_out) = self.generate_matrices_borders(self.M_src, K_m)
            #(self.dst_out, self.neighbours_degrees_dst_out, self.num_k_neighbours_dst_out, self.num_k_edges_dst_out) = self.generate_matrices_borders(self.M_dst, K_m) 
            (self.dst_out, self.neighbours_degrees_dst_out, self.num_k_neighbours_dst_out, self.num_k_edges_dst_out) = self.generate_matrices(self.M_dst, K_m)
        
            (self.src_in, self.neighbours_degrees_src_in, self.num_k_neighbours_src_in, self.num_k_edges_src_in) = self.generate_matrices_borders(self.M_src.T, K_m)
            #(self.dst_in, self.neighbours_degrees_dst_in, self.num_k_neighbours_dst_in, self.num_k_edges_dst_in) = self.generate_matrices_borders(self.M_dst.T, K_m)
            (self.dst_in, self.neighbours_degrees_dst_in, self.num_k_neighbours_dst_in, self.num_k_edges_dst_in) = self.generate_matrices(self.M_dst.T, K_m)
            
        else:
            (self.src_out, self.neighbours_degrees_src_out, self.num_k_neighbours_src_out, self.num_k_edges_src_out) = self.generate_matrices(self.M_src, K_m)
            (self.dst_out, self.neighbours_degrees_dst_out, self.num_k_neighbours_dst_out, self.num_k_edges_dst_out) = self.generate_matrices(self.M_dst, K_m)
        
            (self.src_in, self.neighbours_degrees_src_in, self.num_k_neighbours_src_in, self.num_k_edges_src_in) = self.generate_matrices(self.M_src.T, K_m)
            (self.dst_in, self.neighbours_degrees_dst_in, self.num_k_neighbours_dst_in, self.num_k_edges_dst_in) = self.generate_matrices(self.M_dst.T, K_m)
                
        self.num_vertices_src = len(self.src_out[0])
        self.num_vertices_dst = len(self.dst_out[0])
        
    def topological(self, K_m):
        self.K_m = K_m
        self.topological_prepare(K_m)
        
        T_sim = np.zeros([self.num_vertices_src, self.num_vertices_dst])
               
        for i1 in range(self.num_vertices_src):
            for i2 in range(self.num_vertices_dst):
                T_sim[i1, i2] =  0.5 * sum(min([self.neighbours_degrees_src[k,i1], self.neighbours_degrees_dst[k,i2]]) for k in range(K_m))
                T_sim[i1, i2] += min(self.num_k_neighbours_src[i1], self.num_k_neighbours_dst[i2])
                
                denom = max(self.num_k_neighbours_src[i1], self.num_k_neighbours_dst[i2]) + max(self.num_k_edges_src[i1], self.num_k_edges_dst[i2])
                
                if denom:
                    T_sim[i1, i2] /= denom
                else:
                    T_sim[i1, i2] = 1
                
                
        
        self.T_sim = T_sim
    
    
    def topological_both(self, K_m, weight, remove_borders):
        self.K_m = K_m
        self.topological_prepare_both(K_m, remove_borders)
        
        T_sim = np.zeros([self.num_vertices_src, self.num_vertices_dst])
        
        #weight = 10
        
        
        for i1 in range(self.num_vertices_src):
            for i2 in range(self.num_vertices_dst):
                #x =  sum([((self.neighbours_degrees_src[k,i1] - self.neighbours_degrees_dst[k,i2])/max(self.neighbours_degrees_src[k,i1],self.neighbours_degrees_dst[k,i2]))**2 for k in range(K_m)])/K_m
                #x += ((self.num_k_neighbours_src[i1] - self.num_k_neighbours_dst[i2])/(max(self.num_k_neighbours_src[i1],self.num_k_neighbours_dst[i2])))**2
                #x += ((self.num_k_edges_src[i1] - self.num_k_edges_dst[i2])/max(self.num_k_edges_src[i1], self.num_k_edges_dst[i2]))**2
                   
                x = 0
                
                                
                # output connectivity
                
                for k in range(K_m):
                    a = min(self.neighbours_degrees_src_out[k,i1],self.neighbours_degrees_dst_out[k,i2])
                    b = max(self.neighbours_degrees_src_out[k,i1],self.neighbours_degrees_dst_out[k,i2])
                    if b != 0 and a >= 0:
                        x += (K_m - k) * weight * a/b
                    elif a == -1:
                        x += 1
                    else: #a = b = 0
                        x += (K_m - k) * weight
                        
                x /= K_m
                
                a = min(self.num_k_neighbours_src_out[i1],self.num_k_neighbours_dst_out[i2])                 
                b = max(self.num_k_neighbours_src_out[i1],self.num_k_neighbours_dst_out[i2])
                if b != 0 and a >= 0:
                    x += weight * a/b
                elif a == -1:
                    x += 1
                else: #a = b = 0
                    x += weight
                
                a = min(self.num_k_edges_src_out[i1],self.num_k_edges_dst_out[i2])
                b = max(self.num_k_edges_src_out[i1],self.num_k_edges_dst_out[i2])
                if b != 0 and a >= 0:
                    x += weight * a/b
                elif a == -1:
                    x += 1
                else: #a = b = 0
                    x += weight
                
                # input connectivity
                y = 0
                
                for k in range(K_m):
                    a = min(self.neighbours_degrees_src_in[k,i1],self.neighbours_degrees_dst_in[k,i2])
                    b = max(self.neighbours_degrees_src_in[k,i1],self.neighbours_degrees_dst_in[k,i2])
                    if b != 0 and a >= 0:
                        y += (K_m - k) * weight * a/b
                    elif a == -1:
                        y += 1
                    else: #a = b = 0
                        y += (K_m - k) * weight
                y /= K_m
                
                a = min(self.num_k_neighbours_src_in[i1],self.num_k_neighbours_dst_in[i2])                 
                b = max(self.num_k_neighbours_src_in[i1],self.num_k_neighbours_dst_in[i2])
                if b != 0 and a >= 0:
                    y += weight * a/b
                elif a == -1:
                    y += 1
                else: #a = b = 0
                    y += weight
                
                a = min(self.num_k_edges_src_in[i1],self.num_k_edges_dst_in[i2])
                b = max(self.num_k_edges_src_in[i1],self.num_k_edges_dst_in[i2])
                if b != 0 and a >= 0:
                    y += weight * a/b
                elif a == -1:
                    y += 1
                else: #a = b = 0
                    y += weight
                
                
                T_sim[i1,i2] = (x + y)/6
                
        
        self.T_sim = T_sim
        
        
    
    
    def compareEC(self, ec1, ec2):
       if (ec1.count(".") >= 2) or (ec2.count(".") >= 2):
            ec1 = ec1.split(".")
            ec2 = ec2.split(".")
            if ec1 == ec2:
                return 1
            elif ec1[:-1] == ec2[:-1]:           
                return 0.75
            elif ec1[:-2] == ec2[:-2]:
                return 0.5
            elif ec1[:-3] == ec2[:-3]:
                return 0.25
            else:
                return 0
       elif ec1 == ec2:
            return 1
       else:
            return 0
        
        
    def homological(self, weight_EC, weight_KEGG, additive):
        if weight_EC:
            self.homological_EC(weight_EC, additive)
        if weight_KEGG:
            self.homological_meta_KEGG(weight_KEGG, additive)
        #pass
    
    def homological_EC(self, weight, additive):
       
        
        #print("*****************")
        #print("EC src")
        #print(self.EC_src)
        #print("*****************")
        #print("EC dst")
        #print(self.EC_dst)
        #print("*****************") 
        
        T_sim = self.T_sim         
        """
        for i, x in enumerate(T_sim):
            EC1 = self.EC_src[i]
            for j, y in enumerate(x):
                if y or additive:
                    EC2 = self.EC_dst[j]                    
                    score = []
                    for e1 in EC1:
                        for e2 in EC2:
                            score.append(self.compareEC(e1, e2))
                    if score:
                        if additive:
                            T_sim[i,j] = T_sim[i,j] + max(score) * weight
                        else:
                            T_sim[i,j] = T_sim[i,j] * max(score) * weight
        """                     
                    #elif additive == 0 and (EC1 or EC2): # one is defined the other one is not
                    #    T_sim[i,j] = T_sim[i,j] * 0.5
       
        
        #f = open("data/ec_scores.txt", 'w')
        """
        for i, x in enumerate(T_sim):                                
            EC1 = self.EC_src[i]
            if EC1:
                for j, y in enumerate(x):                                    
                    if y or additive:
                        EC2 = self.EC_dst[j]                    
                        if EC2:    
                            
                            if len(EC1) > len(EC2):
                                EC_longer = EC1
                                EC_shorter = EC2
                            else:
                                EC_longer = EC2
                                EC_shorter = EC1                       
                                                        
                            score = 0
                            for e1 in EC_longer:
                                partial_score = []    
                                for e2 in EC_shorter:
                                    partial_score.append(self.compareEC(e1, e2))
                                score += max(partial_score)
                                                        
                            
                            score /= len(EC_longer)
                            #f.write("EC1 --- " + str(EC1) + ";" + "EC2 --- " + str(EC2) + "::: score ---" + str(score) + "\n")
                            
                            if additive:
                                T_sim[i,j] = T_sim[i,j] + score * weight
                            else:
                                T_sim[i,j] = T_sim[i,j] * score * weight
        
        """
        #f.close()
        
        
        #f = open("data/ec_scores.txt", 'w')
        for i, x in enumerate(T_sim):                                
            EC1 = self.EC_src[i]
            if EC1:
                for j, y in enumerate(x):                                    
                    if y or additive:
                        EC2 = self.EC_dst[j]                    
                        if EC2:    
                            EC_sim = np.zeros([len(EC1), len(EC2)])
                            
                            for i1, e1 in enumerate(EC1):
                                for i2, e2 in enumerate(EC2):
                                    EC_sim[i1, i2] = self.compareEC(e1, e2)
                           
                            
                            # sum the maximal elements of each row
                            # sum the maximal elements of each column
                            # divide by the sum of number of rows and number of columns
                            score = (EC_sim.max(axis = 1).sum() + EC_sim.max(axis = 0).sum())/sum(EC_sim.shape)
                                                               
                            
                            #f.write("EC1 --- " + str(EC1) + ";" + "EC2 --- " + str(EC2) + "::: score ---" + str(score) + "\n")
                            
                            if additive:
                                T_sim[i,j] = T_sim[i,j] + score * weight
                            else:
                                T_sim[i,j] = T_sim[i,j] * score * weight
        
        
        #f.close()
                        
                    
    def homological_meta_KEGG(self, weight, additive):
        T_sim = self.T_sim         
        
        #print("*****************")
        #print("KEGG metabolite src")
        #print(self.KEGG_metabolites_src)
        #print("*****************")
        #print("KEGG metabolite dst")
        #print(self.KEGG_metabolites_dst)
        #print("*****************")
        
        simple = 0
        
        if simple:        
            for i, x in enumerate(T_sim):
                reactants1, products1 = self.KEGG_metabolites_src[i]
                for j, y in enumerate(x):
                    if y or additive:
                        reactants2, products2 = self.KEGG_metabolites_dst[j]
                        score = []
                        #print("***")
                        #print(reactants1,reactants2)
                        if reactants1 == reactants2:
                            score.append(1)
                        elif reactants1 & reactants2 == reactants1:
                            score.append(0.75)
                        elif reactants1 & reactants2:
                            score.append(0.5)
                        else:
                            score.append(1/weight)
                        
                        #print("***")
                        #print(products1,products2)
                        if products1 == products2:
                            score.append(1)
                        elif products1 & products2 == products1:
                            score.append(0.75)
                        elif products1 & products2:
                            score.append(0.5)
                        else:
                            score.append(1/weight)
                            
                        #print(score)
                                                                    
                        if score:
                            if additive:
                                T_sim[i,j] = T_sim[i,j] + min(score) * weight
                            else:
                                T_sim[i,j] = T_sim[i,j] * min(score) * weight
        else:
        # upostevanje okolice
            #f1 = open("data/kegg_scores.txt", 'w')
            #f2 = open("data/kegg_scores_reactions.txt", 'w')
            for i, x in enumerate(T_sim):
                reactants1, products1 = self.KEGG_metabolites_src[i]
                for j, y in enumerate(x):
                        if y or additive:
                            score = 0
                            #score = []
                            
                            #print("******************")
                            #print("******************")
                            #print("******************")
                            #print(i,j)
                            #print("******************")
                            #print("******************")
                            #print("******************")
                            
                           
                            
                            # primerjava cele sosescine
                            reactants2 = self.KEGG_metabolites_dst[j][0].copy()
                            products2 = self.KEGG_metabolites_dst[j][1].copy()
                            for k in range(-1, self.K_m):
                                if k > -1:
                                    for ii,a in enumerate(self.dst_in[k][j]):  
                                        if a == 1:
                                            r2, _ = self.KEGG_metabolites_dst[ii]
                                            reactants2 |= r2
                                    for ii,a in enumerate(self.dst_out[k][j]):  
                                        if a == 1:
                                            _, p2 = self.KEGG_metabolites_dst[ii]
                                            products2 |= p2
                                
                                #print("***************")
                                #print("reactants1")
                                #print(reactants1)
                                #print("reactants2")
                                #print(reactants2)
                                #print("products1")
                                #print(products1)
                                #print("products2")
                                #print(products2)
                                
                                
                                
                                if reactants1 == reactants2 and products1 == products2:
                                    #score.append(1/(k+2))
                                    score = 10/(k+2)
                                    break
                                    
                                if reactants1 & reactants2 == reactants1 and products1 & products2 == products1:
                                    #score.append(0.75 * 1/(k+2))
                                    score = 7.5/(k + 2)
                                    break
                                
                                if reactants1 & reactants2 and products1 & products2:
                                    #score.append(0.25 * 1/(k+2))
                                    score = 5/(k + 2)
                                    break
                                    
                                    
                                #print(score)
                            if score:
                                if additive:
                                    #T_sim[i,j] = T_sim[i,j] + max(score) * weight
                                    T_sim[i,j] = T_sim[i,j] + score * weight
                                else:
                                    #T_sim[i,j] = T_sim[i,j] * max(score) * weight
                                    T_sim[i,j] = T_sim[i,j] * score * weight
                                
                            #f1.write("reactants1 --- " + str(reactants1) + ";" + "reactants2 --- " + str(reactants2) + " *** ")
                            #f1.write("products1 --- " + str(products1) + ";" + "products2 --- " + str(products2) + "::: score ---" + str(score) + "\n")
                            #f2.write(self.reactions_src[i] + " : " + self.reactions_dst[j] + " --- " + str(score) + "\n")
                 
                            
            #f1.close()
            #f2.close()
            
                     
        
    
    # ********************* #
    # ********************* #
    # ********************* #
    # ******** ! ********** #
    # ********************* #
    # ********************* #
    # ********************* #
    def get_anchors(self):
        alignment, _ = self.make_alignment("Greedy", K_m = 1, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 0, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0, remove_too_distant = 0, remove_too_close = 0, distance_to_remove = 1, max_from_neighbours = 0)
        print(alignment)
        return alignment
            
    
    
    # ********************* #
    # ********************* #
    # ********************* #
    # ******** ! ********** #
    # ********************* #
    # ********************* #
    # ********************* #
    def make_alignment(self, alignment_type, K_m, to_many = 0, from_many = 0, remove_borders = 1, topological_weight = 0, homological_weight_EC = 0, homological_weight_KEGG = 0, homological_additive = 0, remove_too_distant = 1, remove_too_close = 1, distance_to_remove = 3, max_from_neighbours = 1):
       
    # max from neighbours: next alignment is obtained from the K_m neighbourhood of the already aligned nodes!!!
        
        if topological_weight:
            if self.K_m != K_m:
                self.topological_both(K_m, topological_weight, remove_borders)
        else:
            if self.K_m != K_m:
                self.K_m = K_m
                self.topological_prepare_both(K_m, remove_borders)
            self.T_sim = np.ones([self.num_vertices_src, self.num_vertices_dst]) 
         
        
        
        self.homological(homological_weight_EC, homological_weight_KEGG, homological_additive)
            
            
        sims = self.T_sim.copy()
                
        alignment = {}
        
        if remove_too_distant or remove_too_close:
            if distance_to_remove != self.distance_to_remove:
                self.distance_to_remove = distance_to_remove
            
                if distance_to_remove == K_m:
                    self.src_in_remove = self.src_in
                    self.src_out_remove = self.src_out
                    self.dst_in_remove = self.dst_in
                    self.dst_out_remove = self.dst_out
                else:
                    
                    n_src = len(self.src_in[0])
                    print(n_src)
                    n_dst = len(self.dst_in[0])
                    print(n_dst)
                    
                    self.src_in_remove = self.src_in
                    self.src_out_remove = self.src_out
                    self.dst_in_remove = self.dst_in
                    self.dst_out_remove = self.dst_out
                    
                    self.src_in_remove = np.zeros([distance_to_remove, n_src, n_src])
                    self.src_out_remove = np.zeros([distance_to_remove, n_src, n_src])
                    self.dst_in_remove = np.zeros([distance_to_remove, n_dst, n_dst])
                    self.dst_out_remove = np.zeros([distance_to_remove, n_dst, n_dst])
                    
                    self.src_in_remove[0] = self.src_in[0]
                    self.src_out_remove[0] = self.src_out[0]
                    self.dst_in_remove[0] = self.dst_in[0]
                    self.dst_out_remove[0] = self.dst_out[0]
            
                    # matrični produkt: np.dot(M_p, M_p)
                    for k in range(1, distance_to_remove):
                        self.src_in_remove[k] = np.array(np.array(np.dot(self.src_in_remove[k-1], self.src_in_remove[0]), dtype = bool), dtype = int)
                        self.src_out_remove[k] = np.array(np.array(np.dot(self.src_out_remove[k-1], self.src_out_remove[0]), dtype = bool), dtype = int)
                        self.dst_in_remove[k] = np.array(np.array(np.dot(self.dst_in_remove[k-1], self.dst_in_remove[0]), dtype = bool), dtype = int)
                        self.dst_out_remove[k] = np.array(np.array(np.dot(self.dst_out_remove[k-1], self.dst_out_remove[0]), dtype = bool), dtype = int)
            
        
        while sims.any():
            
            # ["Simple", "Greedy", "Probabilistic"]
            if alignment_type == "Greedy":                
                (i,j) = self.align_pair_max_opt(sims, alignment, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove, max_from_neighbours)
            elif alignment_type =="Greedy old":
                (i,j) = self.align_pair_max(sims, alignment, from_many)
                self.adjust_similarities(sims, i, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove)
            elif alignment_type == "Simple":
                (i,j) = self.align_pair_simple(sims)
                self.adjust_similarities(sims, i, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove)
        
            else:
                (i,j) = self.align_random(sims)
                self.adjust_similarities(sims, i, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove)
        
            """
            print("*********")
            print("source", i, '=', self.reaction_names_src[i])
            print("destination", j, '=', self.reaction_names_dst[j])
            print("score = ", sims[i,j])
            print("*********")
            """
            alignment[i] = j 
            
        
        
        score = sum(self.T_sim[i,alignment[i]] for i in alignment)
        self.last_alignment = alignment
        
        return alignment, score
    
    
    # find the projection with a maximal score in the similarity matrix sims
    # if multiple projections have the same maximal score select the projections that have the least repetitions of the source node
    # among these select the one that reduces the similarity scores the least
    
    def align_pair_max_opt(self, sims, alignment, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove, max_from_neighbours):
        
        
        # candidates should be limited to the K_m neighbourhoud of the aligned nodes!!!
        if max_from_neighbours:
            if alignment:
                neighbours_aligned = []
                for a in alignment:
                    for s in range(self.num_vertices_src):
                        if s not in alignment or from_many: # ce je s ze vkljucen, pogledamo ce je from_many 1 - v tem primeru ga lahko preslikamo veckrat
                            if s not in neighbours_aligned and (any(self.src_out[:,a,s]) or any(self.src_in[:,a,s]) or s == a):
                                neighbours_aligned.append(s)
                
            else:
                neighbours_aligned = range(self.num_vertices_src)
            
            if neighbours_aligned and np.amax(sims[neighbours_aligned,:]) > 0:
                candidates = np.argwhere(sims[neighbours_aligned,:] == np.amax(sims[neighbours_aligned,:])) # find the indices with maximal scores
                candidates[:,0] = np.array(neighbours_aligned)[candidates[:,0]] # recalculate the indices             
            else:
                candidates = np.argwhere(sims == np.amax(sims)) # find the indices with maximal scores
        else:
            candidates = np.argwhere(sims == np.amax(sims)) # find the indices with maximal scores
        
        # candidates now includes the possible candidates for the next projection
        # canidates is a 2D ndarray
        # candidates = [[source, destination]]
        
        
        # dictionary: how many times does the same reaction map with the maximal value?
        l = defaultdict(int)
        for c in candidates[:,0]:
            l[c] += 1
        
        m = min(l.values())
        i = []
        for k in l:
            if l[k] == m:
                i.append(k)
        #i now includes the elementes with the least repetitions
        
        #i = min(l.items(), key=operator.itemgetter(1))[0]   # tisti, ki ima najmanj ponovitev
        
        # only one solution
        if len(i) == 1:     
            i = i[0]
            j = (candidates[np.argwhere(i == candidates[:,0])]).flatten()[1] 
            self.adjust_similarities(sims, i, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove)
        # multiple solutions
        # select the solution that reduces the similarities the least (call adjust_similarities with force_zeros = 1)
        else:
            best_avg = -1
            best_max = -1
            best_i = -1
            best_j = -1
            
            for i1 in i:
                j = (candidates[np.argwhere(i1 == candidates[:,0])]).flatten()[1] 
                s_tmp = sims.copy()
                self.adjust_similarities(s_tmp, i1, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove, force_zeros = 1)
                s_avg = np.average(s_tmp)
                s_max = np.max(s_tmp)
                if (s_max > best_max) or (s_max == best_max and s_avg > best_avg):
                    best_max = s_max
                    best_avg = s_avg
                    best_i = i1
                    best_j = j
                    
            
            i = best_i
            j = best_j
            """
            print("##########")
            print(i)
            print("##########")
            """
            #i = i[0]
            #j = (candidates[np.argwhere(i == candidates[:,0])]).flatten()[1] 
            self.adjust_similarities(sims, i, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove)
        """
        print("**************")
        print(candidates)   
        print(i,j)
        print("**************")
        """
        return i,j            
    

    def align_pair_max(self, sims, alignment, from_many):
        # candidates should be limited to the K_m neighbourhoud of the aligned nodes!!!
        if alignment:
            neighbours_aligned = []
            for a in alignment:
                for s in range(self.num_vertices_src):
                    if s not in alignment or from_many: # ce je s ze vkljucen, pogledamo ce je from_many 1 - v tem primeru ga lahko preslikamo veckrat
                        if s not in neighbours_aligned and (any(self.src_out[:,a,s]) or any(self.src_in[:,a,s]) or s == a):
                            neighbours_aligned.append(s)
            
        else:
            neighbours_aligned = range(self.num_vertices_src)
        
        if neighbours_aligned and np.amax(sims[neighbours_aligned,:]) > 0:
            candidates = np.argwhere(sims[neighbours_aligned,:] == np.amax(sims[neighbours_aligned,:]))
            candidates[:,0] = np.array(neighbours_aligned)[candidates[:,0]]                
        else:
            candidates = np.argwhere(sims == np.amax(sims))
            
        l = defaultdict(int)
        for c in candidates[:,0]:
            l[c] += 1
        
        i = min(l.items(), key=operator.itemgetter(1))[0]   # tisti, ki ima najmanj ponovitev
        j = (candidates[np.argwhere(i == candidates[:,0])]).flatten()[1] 
            
        return i,j            
    
    def align_pair_simple(self, sims):
        i,j = np.unravel_index(sims.argmax(), sims.shape)
        return i,j
    
    def align_random(self, sims):
        r = np.random.random() * sims.sum()
        cum = 0
        for i,x in enumerate(sims):
            for j,y in enumerate(x):
                if y > 0:
                    cum += y
                    if cum >= r:
                        return i,j
    
    
    def adjust_similarities(self, sims, i, j, from_many, to_many, remove_too_distant, remove_too_close, distance_to_remove, force_zeros = 0):        
        
        
        weight = 0.1
        
        if remove_too_distant:
            
            
            if force_zeros:
                #############################
                # ODSTRANI PREVEC ODDALJENE #
                #############################
                
                # 
                # postavi na 0 vse preslikave, ki so neskladne z izbrano: oddaljena vozlisca se preslikajo v bliznja in obratno
                # s je sosed od i
                
                #       ce je kandidat od j oddaljen vec kot k, mu daj ujemanje na 0
                
                
                # ZVEZNO MANJŠANJE UJEMANJA Z RAZDALJO OD ŽE PRESLIKANIH
               
                
                # pojdi cez vse sosede od izvora dodane preslikave (i)
                for s in range(self.num_vertices_src):
                    #print(i)
                    if self.src_out_remove[0,i,s] > 0 or i == s:   # ali je sosed ali on?
                        #   pojdi cez kanidate za preslikavo sosedov
                        for c,y in enumerate(sims[s,:]):
                            # ali je preslikava vecja od 0 in ali je od ponora dodane preslikave (j) oddaljen vec kot K_m
                            if y > 0 and self.dst_out_remove[:,j,c].sum() == 0:
                            #c not in neighbours(j):
                                #sims[s,c] = sims[c,s] = 0
                                sims[s,c] = 0
                              
                for s in range(self.num_vertices_src):
                    if self.src_in_remove[0,i,s] > 0 or i == s:   # ali je sosed ali on?
                        #   pojdi cez kanidate za preslikavo sosedov
                        for c,y in enumerate(sims[s,:]):
                            # ali je preslikava vecja od 0 in ali je od ponora dodane preslikave (j) oddaljen vec kot K_m
                            if y > 0 and self.dst_in_remove[:,j,c].sum() == 0:
                            #c not in neighbours(j):
                                #sims[s,c] = sims[c,s] =  0
                                sims[s,c] = 0
            else:
                # pojdi cez vse sosede od izvora dodane preslikave (i)
                for s in range(self.num_vertices_src):
                  if self.src_out_remove[0,i,s] > 0 or self.src_in_remove[0,i,s] > 0 or i == s:   # ali je sosed ali on?
                        #   pojdi cez kanidate za preslikavo sosedov
                        for c,y in enumerate(sims[s,:]):
                            # ali je preslikava vecja od 0 in ali je od ponora dodane preslikave (j) oddaljen manj kot K_m
                            if y > 0:
                                for k in range(distance_to_remove):
                                    if self.dst_out_remove[k,j,c] == 1 or self.dst_in_remove[k,j,c] == 1:
                                        sims[s,c] *= weight/(k+1) + 1
                                        break
                              
                
        if remove_too_close:
            if force_zeros:
                ##############################
                # ODSTRANI PREMALO ODDALJENE #
                ##############################            
                    
                # pojdi cez vse sosede, ki so od izvora dodane preslikave i oddaljeni vec kot K_m
                for s in range(self.num_vertices_src):
                    if s != i:
                        if self.src_out_remove[:,i,s].sum() == 0: # not in the neighbourhood
                            for c,y in enumerate(sims[s,:]):
                                if y > 0 and (self.dst_out_remove[0,j,c] == 1 or j == c): # ne morejo se preslikati v vozlisca, ki so sosedi od preslikanega vozlisca j
                                    #sims[s,c] = sims[c,s] = 0
                                    sims[s,c] = 0
                            
                for s in range(self.num_vertices_src):
                    if s != i:
                        if self.src_in_remove[:,i,s].sum() == 0: # not in the neighbourhood
                            for c,y in enumerate(sims[s,:]):
                                if y > 0 and (self.dst_in_remove[0,j,c] == 1 or j == c): # ne morejo se preslikati v vozlisca, ki so sosedi od preslikanega vozlisca j
                                    #sims[s,c] = sims[c,s] = 0
                                    sims[s,c] = 0
            else:
                # pojdi cez vse sosede, ki so od izvora dodane preslikave i oddaljeni vec kot K_m
                for s in range(self.num_vertices_src):
                    if s != i:
                        if self.src_out_remove[:,i,s].sum() == 0 or self.src_in_remove[:,i,s].sum() == 0: # not in the neighbourhood
                            for c,y in enumerate(sims[s,:]):
                                if y > 0 and (self.dst_out_remove[0,j,c] == 0) and (self.dst_in_remove[0,j,c] == 0): # not a neighbour
                                    #sims[s,c] = sims[c,s] = 0
                                    sims[s,c] *= (weight + 1)
            
            
        # preferiraj tista vozlisca, v katera preslikava se ni bila narejena - ce je vec enakih vrednosti, vzemi se nepreslikane
    
        sims[i,j] = 0
        
        if from_many == 0:
            sims[i,:] = 0
                
        if to_many == 0:
            sims[:,j] = 0
       
    def align_probs(self, K_m):
        if self.K_m != K_m:
            #self.topological(K_m)
            self.topological2(K_m)
            self.homological()
        
        sims = self.T_sim.copy()
        score = 0
        i_src = list(range(len(sims[:,0])))
        i_dst = list(range(len(sims[0,:])))
        
      
        alignment = {}
        
     
        while sims.any():
                        
            a,b = self.align_random(sims)            
            i,j = i_src[a], i_dst[b]
            
            # 
            # postavi na 0 vse preslikave, ki so neskladne z izbrano: oddaljena vozlisca se preslikajo v bliznja in obratno
            # s je sosed od i
            
            #       ce je kandidat od j oddaljen vec kot k, mu daj ujemanje na 0
            
            # pojdi cez vse sosede od izvora dodane preslikave (i)
            for s, s_full in enumerate(i_src):
                if self.src[0, i, s_full] > 0:  # ali je sosed?
                    #   pojdi cez kanidate za preslikavo sosedov
                    for c,y in enumerate(sims[s,:]):
                        c_full = i_dst[c]
                        # ali je preslikava vecja od 0 in ali je od ponora dodane preslikave (j) oddaljen vec kot K_m
                        if y > 0 and self.dst[:,j,c_full].sum() == 0:
                        #c not in neighbours(j):
                            sims[s,c] = 0
                                        
            # TODO
            # preferiraj tista vozlisca, v katera preslikava se ni bila narejena - ce je vec enakih vrednosti, vzemi se nepreslikane
            
           
            alignment[i] = j    
            score += self.T_sim[i,j]
            # T_sim[i,:] = 0
            sims = np.delete(sims, a, axis=0)
            del i_src[a]
        
        return alignment, score
    
    def print_T_sim(self, file_name = "", tsv = 0):
        tmp = ""
        if tsv:
            tmp += "\t"
            for i in self.reactions_dst:
                tmp += i + "\t"
            tmp +="\n"
            for i, x in enumerate(self.T_sim):
                
                tmp += self.reactions_src[i] + "\t"
                for j in x:
                    tmp += str(j) + "\t"
                tmp += "\n"
        
        else:
            tmp +=" " * 25
            for i in self.reactions_dst:
                tmp += "{:>25}".format(i)
            tmp +="\n"
            
            for i, x in enumerate(self.T_sim):
                
                tmp += "{:25}".format(self.reactions_src[i])
                for j in x:
                    tmp += "{:>25.4f}".format(j)
                tmp += "\n"
        if file_name:
            f = open(file_name, "w")
            f.write(tmp)
            f.close()
        else:
            print(tmp)
      
    def print_T_sim_max(self, file_name = ""):
        tmp = ""
        for i,x in enumerate(self.T_sim):
            M = -1
            i_M = []
            for j,y in enumerate(x):
                if y > M:
                    M = y
                    i_M = [j]
                elif y == M:
                    i_M.append(j)
            for j in i_M:
                tmp += self.reactions_src[i] + " : " + self.reactions_dst[j] + " - score = " + str(M) + "\n"
        if file_name:
            f = open(file_name, "w")
            f.write(tmp)
            f.close()
        else:
            print(tmp)
                    
        
        
        
    
    def print_matrix(self, destination, in_out):
        if destination == "dst" and in_out == "in":
            M = self.dst_in
            names = self.reactions_dst
        if destination == "dst" and in_out == "out":
            M = self.dst_out
            names = self.reactions_dst
        if destination == "src" and in_out == "in":
            M = self.src_in
            names = self.reactions_src
        if destination == "src" and in_out == "out":
            M = self.src_out
            names = self.reactions_src
        
        print("*"*50)   
        print(destination.upper() + " " + in_out.upper())
        print(" " * 8, end="")
        for i in names:
            print("{:>8}".format(i), end="")
        print()
        
        for i, x in enumerate(M[0]):
            print("{:8}".format(names[i]), end="")
            for j in x:
                print("{:>8}".format(j), end = "")
            print()
        print("*"*50)  
    
    def print_M(self, M):
        for i in M:
            for j in i:
                print(j, end=" ")
            print()
            
            