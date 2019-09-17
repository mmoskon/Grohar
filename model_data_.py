import cobra.test
import cobra
import cobra.flux_analysis.moma as moma
import cobra.flux_analysis.variability as fva
from libsbml import SBMLReader
from collections import defaultdict    
import numpy as np
#import matplotlib.pyplot as plt
#from bioservices.kegg import KEGG

from cobra import Model, Reaction, Metabolite

class model_data:

    def __init__(self):

        # kljuc = ID_reakcije, vsebina = ime reakcije
        self.reactions = {}

        # kljuc = ID_metabolita, vsebina = ime metabolita
        self.metabolites = {}
        
        # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki sodelujejo v reakciji
        self.reaction_metabolites = defaultdict(set)
        
        # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki jih reakcija producira
        self.reaction_products = defaultdict(set)
        
        # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki jih reakcija porablja
        self.reaction_reactants = defaultdict(set)
        
        # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih metabolit sodeluje
        self.metabolite_reactions = defaultdict(set)
        
        # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih je metabolit reaktant
        self.metabolite_as_reactant = defaultdict(set)
        
        # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih je metabolit produkt
        self.metabolite_as_product = defaultdict(set)

        
        # kljuc = ID_compartmenta, vsebina = mnozica ID_metabolitov v tem compartmentu
        self.compartment_metabolites = defaultdict(set)

        # kljuc = ID_compartmenta, vsebina = mnozica ID_reakcij v tem compartmentu
        self.compartment_reactions = defaultdict(set)

        # kljuc = ID reakcije, vsebina = zaporedna stevilka reakcije v cobra modelu
        self.reaction_position = {}
        
        # ID reakcij urejen po indeksu v cobra modelu
        self.listed_reactions = [] 

        # kljuc = ID metabolita, vsebina = zaporedna stevilka metabolita v cobra modelu
        self.metabolite_position = {}
        
        # ID metabolitov urejen po indeksu v cobra modelu
        self.listed_metabolites = [] 

        
  
        # kljuc = ID reackije, vsebina = [lower bound, upper bound]
        self.reaction_bounds = {}

        # kljuc = ID reakcije, vsebina = ime reakcije
        self.exchange_reactions = {}
        self.uptake_reactions = {}
        
        self.pFBA = False
        self.MOMA = False
        self.FVA = False
        
        self.fluxes = None
        self.perturbed_model = None
       
    def load_model_sbml(self, path):
        reader = SBMLReader()
        self.document = reader.readSBML(path)
        self.model = self.document.getModel() 
        self.fill_data_sbml()
        self.type = "sbml"
        
        
    
    def load_model_cobra_sbml(self, path):
        self.model = cobra.io.read_sbml_model(path)
            
        self.fill_data_cobra()
        self.type = "cobra"
        
        
            
    def load_model_mat(self, path):
        self.model =  cobra.io.load_matlab_model(path)                       
        self.fill_data_cobra()
        self.type = "cobra"
             
    def import_model(self, model):
        self.model = model
        self.fill_data_cobra()
        self.type = "cobra"
    
    
    def create_model(self, name):
        self.model = Model(name)
        self.type = "cobra"
                
    def save_model_cobra_sbml(self, file_name):
        try:
            cobra.io.write_sbml_model(self.model, file_name)
        except Exception as inst:
            print(inst)

               
    
    def save_model_cobra_mat(self, file_name):
        try:
            cobra.io.save_matlab_model(self.model, file_name)
        except Exception as inst:
            print(inst)
    
    def add_reaction_basic(self, r_id, r_name):
        if r_id not in self.reactions:
            ##############
            # model_data #
            ##############
            self.reactions[r_id] = r_name
            self.listed_reactions.append(r_id)
            self.reaction_position[r_id] = len(self.listed_reactions) - 1
            
            self.reaction_bounds[r_id] = (-1000, 1000)
            
            #########
            # cobra #
            #########
            reaction = Reaction(r_id)
            self.model.add_reactions([reaction])
            
            reaction.name = r_name
            
            r_metabolites = {}
            reaction.add_metabolites(r_metabolites)
                         
            reaction.lower_bound = -1000
            reaction.upper_bound = 1000
            return 1
        else:
            return 0
            
             
    def add_reaction(self, r_id, r_name, reactants, products, bounds = []):
        if r_id not in self.reactions:
            ##############
            # model_data #
            ##############
            self.reactions[r_id] = r_name
            self.listed_reactions.append(r_id)
            self.reaction_position[r_id] = len(self.listed_reactions) - 1
            R_compartment = set()                                  
                                  
            for m_r in reactants:
                if m_r in self.metabolites:
                    self.reaction_metabolites[r_id].add(m_r)
                    self.reaction_reactants[r_id].add(m_r)
                    self.metabolite_reactions[m_r].add(r_id)    
                    self.metabolite_as_reactant[m_r].add(r_id)    
                    R_compartment.add(self.model.metabolites[self.metabolite_position[m_r]].compartment)
                                      
                 
            for m_p in products:
                if m_p in self.metabolites:
                    self.reaction_metabolites[r_id].add(m_p)
                    self.reaction_products[r_id].add(m_p)
                    self.metabolite_reactions[m_p].add(r_id)
                    self.metabolite_as_product[m_p].add(r_id)
                    R_compartment.add(self.model.metabolites[self.metabolite_position[m_p]].compartment)
                    
            if bounds:
                self.reaction_bounds[r_id] = bounds
            else:
                self.reaction_bounds[r_id] = (-1000, 1000)
            
            if len(reactants) == 1 and len(products) == 0:
                self.exchange_reactions[r_id] = r_name
                if bounds and bounds[0] < 0:
                    self.uptake_reactions[r_id] = r_name                 
            
            R_compartment = list(R_compartment)
            R_compartment.sort()            
            self.compartment_reactions[" ".join(i for i in R_compartment)].add(r_id)
            
            #########
            # cobra #
            #########
            reaction = Reaction(r_id)
            self.model.add_reactions([reaction])
            
            reaction.name = r_name
            
            r_metabolites = {}
            for m_r in reactants:
                if m_r in self.metabolites:
                    r_metabolites[m_r] = -1
            for m_p in products:
                if m_p in self.metabolites:
                    r_metabolites[m_p] = 1
            
            reaction.add_metabolites(r_metabolites)
            
            if bounds:
                reaction.lower_bound = bounds[0]
                reaction.upper_bound = bounds[1]  
            else:                   
                reaction.lower_bound = -1000
                reaction.upper_bound = 1000
            
            return 1
        
        else:
            return 0
   
    def delete_reaction(self, r_id):
        if r_id in self.reactions:
            # TODO
            if self.model.reactions[self.reaction_position[r_id]] not in self.model.objective:
                self.delete_reaction_model_data(r_id)
                self.delete_reaction_cobra(r_id)
                return 1
            else:
                print("cannot delete the reaction composing objective function")
        return 0

     
    def delete_reaction_model_data(self, r_id):
        if r_id in self.reactions:
            # kljuc = ID_reakcije, vsebina = ime reakcije
            del self.reactions[r_id]
                        
            # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki sodelujejo v reakciji
            if r_id in self.reaction_metabolites:
                del self.reaction_metabolites[r_id]
            
            # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki jih reakcija producira
            if r_id in self.reaction_products:
                del self.reaction_products[r_id]
            
            # kljuc = ID_reakcije, vsebina = mnozica ID_metabolit-ov, ki jih reakcija porablja
            if r_id in self.reaction_reactants:
                del self.reaction_reactants[r_id]
            
            # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih metabolit sodeluje
            for m in self.metabolite_reactions:
                if r_id in self.metabolite_reactions[m]:
                    self.metabolite_reactions[m].remove(r_id)                
            
            # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih je metabolit reaktant
            for m in self.metabolite_as_reactant:
                if r_id in self.metabolite_as_reactant[m]:
                    self.metabolite_as_reactant[m].remove(r_id)
            
            # kljuc = ID_metabolita, vsebina = mnozica ID_reakcij, v katerih je metabolit produkt
            for m in self.metabolite_as_product:
                if r_id in self.metabolite_as_product[m]:
                    self.metabolite_as_product[m].remove(r_id)
               
            # kljuc = ID_compartmenta, vsebina = mnozica ID_reakcij v tem compartmentu
            for c in self.compartment_reactions:
                if r_id in self.compartment_reactions[c]:
                    self.compartment_reactions[c].remove(r_id)
                
            # kljuc = ID reackije, vsebina = [lower bound, upper bound]
            if r_id in self.reaction_bounds:
                del self.reaction_bounds[r_id]
    
            # kljuc = ID reakcije, vsebina = ime reakcije
            if r_id in self.exchange_reactions:
                del self.exchange_reactions[r_id]
            if r_id in self.uptake_reactions:
                del self.uptake_reactions[r_id]
          
        
        
    def delete_reaction_cobra(self, r_id):
        
        idx = self.reaction_position[r_id]
        #idx = self.listed_reactions.index(r_id)
        
        # ID reakcij urejen po indeksu v cobra modelu
        del self.listed_reactions[idx]
                
        # kljuc = ID reakcije, vsebina = zaporedna stevilka reakcije v cobra modelu
        self.reaction_position = {}
        for i,r in enumerate(self.listed_reactions):
            self.reaction_position[r] = i
                                  
        #del self.model.reactions[idx]                                  
        self.model.remove_reactions([r_id])
            
        
        
    
    def edit_reaction(self, r_id, r_name, reactants, products, bounds = []):
        if r_id in self.reactions:
            self.delete_reaction_model_data(r_id)
            
            ##############
            # model_data #
            ##############
            self.reactions[r_id] = r_name
            #self.listed_reactions.append(r_id)
            #self.reaction_position[r_id] = len(self.listed_reactions) - 1
            R_compartment = set()                                  
                                  
            for m_r in reactants:
                if m_r in self.metabolites:
                    self.reaction_metabolites[r_id].add(m_r)
                    self.reaction_reactants[r_id].add(m_r)
                    self.metabolite_reactions[m_r].add(r_id)    
                    self.metabolite_as_reactant[m_r].add(r_id)    
                    R_compartment.add(self.model.metabolites[self.metabolite_position[m_r]].compartment)
                                      
                 
            for m_p in products:
                if m_p in self.metabolites:
                    self.reaction_metabolites[r_id].add(m_p)
                    self.reaction_products[r_id].add(m_p)
                    self.metabolite_reactions[m_p].add(r_id)
                    self.metabolite_as_product[m_p].add(r_id)
                    R_compartment.add(self.model.metabolites[self.metabolite_position[m_p]].compartment)
                    
            if bounds:
                self.reaction_bounds[r_id] = bounds
            
            if len(reactants) == 1 and len(products) == 0:
                self.exchange_reactions[r_id] = r_name
                if bounds and bounds[0] < 0:
                    self.uptake_reactions[r_id] = r_name                 
            
            R_compartment = list(R_compartment)
            R_compartment.sort()            
            self.compartment_reactions[" ".join(i for i in R_compartment)].add(r_id)
            
            #########
            # cobra #
            #########
            reaction =  self.model.reactions[self.reaction_position[r_id]]
            reaction.clear_metabolites()
            
            reaction.name = r_name
            
            r_metabolites = {}
            for m_r in reactants:
                if m_r in self.metabolites:
                    r_metabolites[m_r] = -1
            for m_p in products:
                if m_p in self.metabolites:
                    r_metabolites[m_p] = 1
            
            reaction.add_metabolites(r_metabolites)
            #reaction.metabolites = r_metabolites
            
          
            
            if bounds:
                reaction.lower_bound = bounds[0]
                reaction.upper_bound = bounds[1]    

            return 1
        else:
            return 0                        
        
    def delete_metabolite(self, m_id):
        if m_id in self.metabolites:
            if m_id not in self.metabolite_reactions:
                # kljuc = ID_metabolita, vsebina = ime metabolita
                del self.metabolites[m_id]
                              
                # kljuc = ID_compartmenta, vsebina = mnozica ID_metabolitov v tem compartmentu
                for c in self.compartment_metabolites:
                    if m_id in self.compartment_metabolites[c]:
                        self.compartment_metabolites[c].remove(m_id)
        
                
                idx = self.metabolite_position[m_id]
                               
                # ID metabolitov urejen po indeksu v cobra modelu
                del self.listed_metabolites[idx] 
                  
                 # kljuc = ID metabolita, vsebina = zaporedna stevilka metabolita v cobra modelu
                self.metabolite_position = {}
                for i,m in enumerate(self.listed_metabolites):
                    self.metabolite_position[m_id] = i
                                          
                #del self.model.metabolites[idx]                                                  
                self.model.remove_metabolites([m_id])
                
                
                return 1
            else:
                print("cannot delete; metabolite in reactions", self.metabolite_reactions[m_id])
               
            
            return 0

        
    
    
    def add_metabolite(self, m_id, m_name, compartment):
        if m_id not in self.metabolites:
            
            ##############
            # model_data #
            ##############
            self.metabolites[m_id] = m_name
            self.listed_metabolites.append(m_id)
            self.metabolite_position[m_id] = len(self.listed_metabolites) - 1
            self.compartment_metabolites[compartment].add(m_id)
                                    
            #########
            # cobra #
            #########
            if compartment not in self.model.compartments:
                self.model.compartments[compartment] = compartment
            
            metabolite = Metabolite(m_id, formula="", name = m_name, compartment = compartment)                                    
            #self.model.metabolites.append(metabolite)
            self.model.add_metabolites([metabolite])
            
            return 1
        
        return 0
    
    def edit_metabolite(self, m_id, m_name):  
         if m_id in self.metabolites:
            
            ##############
            # model_data #
            ##############
            self.metabolites[m_id] = m_name
            #self.listed_metabolites.append(m_id)
            #self.metabolite_position[m_id] = len(self.listed_metabolites) - 1
            #self.compartment_metabolites[compartment].add(m_id)
                                    
            #########
            # cobra #
            #########
            #if compartment not in self.model.compartments:
            #    self.model.compartments[compartment] = compartment
            
            #metabolite = Metabolite(m_id, formula="", name = m_name, compartment = compartment)                                    
            
            self.model.metabolites[self.metabolite_position[m_id]].name = m_name
                                  
            return 1
        
         return 0                     
     
    def delete_compartemnt(self, c_id):
        if c_id in self.compartments:
            if self.compartment_metabolites[c_id] or self.compartment_reactions[c_id]:
                print("cannot delete; compartment is not empty")
            else:
                del self.compartments[c_id]
                del self.model.compartments[c_id]
                return 1
        return 0
        
    def add_compartment(self, c_id, c_name):
        if c_id not in self.compartments:
            self.compartments[c_id] = c_name
            self.model.compartments[c_id] = c_name
            return 1
        return 0
           
           
        
    def fill_data_sbml(self):
        model = self.model
        
        self.compartments = {i.getId():i.getName() for i in model.getListOfCompartments()}
               
        
        
        for R in model.getListOfReactions():
            self.reactions[R.getId()] = R.getName()
        
       
            
        for S in model.getListOfSpecies():
            self.metabolites[S.getId()] = S.getName()
            self.compartment_metabolites[S.getCompartment()].add(S.getId())

               
        
        for R in model.getListOfReactions():
            R_compartment = set()
            for M in R.getListOfReactants():
                self.reaction_metabolites[R.getId()].add(M.getSpecies())
                self.reaction_reactants[R.getId()].add(M.getSpecies())
                                
                self.metabolite_reactions[M.getSpecies()].add(R.getId())
                self.metabolite_as_reactant[M.getSpecies()].add(R.getId())               
                
                R_compartment.add(model.getListOfSpecies().getElementBySId(M.getSpecies()).getCompartment())
                
            
            for M in R.getListOfProducts():
                self.reaction_metabolites[R.getId()].add(M.getSpecies())
                self.reaction_products[R.getId()].add(M.getSpecies())
                
                self.metabolite_reactions[M.getSpecies()].add(R.getId())
                self.metabolite_as_product[M.getSpecies()].add(R.getId())               
                
                R_compartment.add(model.getListOfSpecies().getElementBySId(M.getSpecies()).getCompartment())
            
            R_compartment = list(R_compartment)
            R_compartment.sort()
            self.compartment_reactions[" ".join(i for i in R_compartment)].add(R.getId())
            
                    
    def fill_data_cobra(self):
        model = self.model
        
                 
        if (not self.model.compartments) or (self.model.compartments == None) or (None in self.model.compartments) or ('[' in self.model.compartments):              
           self.set_compartments() 
        else:
           self.compartments = self.model.compartments   
        
        
        
        for (i, R) in enumerate(model.reactions):
            self.reactions[R.id] = R.name
            self.reaction_position[R.id] = i
            self.listed_reactions.append(R.id)
        
        for (i, M) in enumerate(model.metabolites):
            self.metabolites[M.id] = M.name
            self.compartment_metabolites[M.compartment].add(M.id)
            self.metabolite_position[M.id] = i
            self.listed_metabolites.append(M.id)
        
        for R in model.reactions:
            R_compartment = set()
            for M in R.reactants:
                self.reaction_metabolites[R.id].add(M.id)
                self.reaction_reactants[R.id].add(M.id)
                                
                self.metabolite_reactions[M.id].add(R.id)
                self.metabolite_as_reactant[M.id].add(R.id)               
                
                R_compartment.add(M.compartment)
                
            
            for M in R.products:
                self.reaction_metabolites[R.id].add(M.id)
                self.reaction_products[R.id].add(M.id)
                
                self.metabolite_reactions[M.id].add(R.id)
                self.metabolite_as_product[M.id].add(R.id)               
                
                R_compartment.add(M.compartment)
            
            R_compartment = list(R_compartment)
            R_compartment.sort()
            self.compartment_reactions[" ".join(i for i in R_compartment)].add(R.id)

        # reaction boundaries
        for R in model.reactions:
            self.reaction_bounds[R.id] = [R.lower_bound, R.upper_bound]
        
        # uptake and exchange reactions
        for i in self.reaction_metabolites:
            if len(self.reaction_metabolites[i]) == 1:
                self.exchange_reactions[i] = self.reactions[i]
                if self.reaction_bounds[i][0] < 0:
                    self.uptake_reactions[i] = self.reactions[i]

        
        
    def set_lower_bound(self, R_id, l):
        self.reaction_bounds[R_id][0] = l
        i = self.reaction_position[R_id]
        self.model.reactions[i].lower_bound = l        
       
    def set_upper_bound(self, R_id, u):
        self.reaction_bounds[R_id][1] = u
        i = self.reaction_position[R_id]
        self.model.reactions[i].upper_bound = u 
    
    def set_bounds(self, R_id, l, u):
        self.set_lower_bound(R_id, l)
        self.set_upper_bound(R_id, u)
                           
    def set_compartments(self):
        model = self.model
        model.compartments = {}
        
        for M in model.metabolites:
            M.compartment = M.id.split('[')[1].split(']')[0] # asn_L[c]--> c
            model.compartments[M.compartment] = M.compartment                        
        
        
        self.compartments = model.compartments

                
    def run_FBA(self):
        
        self.perturbed_model = None
        
        if self.model.reactions:            
            if not self.model.objective:
                self.model.objective = self.model.reactions[0]
            
                       
            if self.pFBA:
                print("running pFBA")
                self.sol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(self.model)
                print(self.sol)                
            else:
                print("running FBA")
                self.sol = self.model.optimize()
            
            
            self.fluxes = self.sol.x_dict
        
    
    def run_FVA(self, r):
        print("Running FVA for reaction", r)
        if self.perturbed_model:
            b = fva.flux_variability_analysis(self.perturbed_model, [r])[r]
        else:
            b = fva.flux_variability_analysis(self.model, [r])[r]
        
        return [b['minimum'], b['maximum']]
        
                
    
    
    # boundaries = [reactionID, boundary_type, boundary_value]        
    def set_boundaries(self, *boundaries, **objective):
        model = self.model
               
        if self.type == "cobra":
            reaction_position = self.reaction_position
            
            # make perturbations
            for (R, boundary_type, boundary_value) in boundaries:
                i = reaction_position[R]
                if boundary_type == 'l':
                    model.reactions[i].lower_bound = boundary_value
                elif boundary_type == 'u':                    
                    model.reactions[i].upper_bound = boundary_value

            # TODO
            if objective:
                model.objective = model.reactions[reaction_position[objective["objective_id"]]]

            self.run_FBA()
            
            
    def set_objective(self, r_id):
        # TODO
        self.model.objective = self.model.reactions[self.reaction_position[r_id]]
        
        
    # boundaries = [reactionID, boundary_type, boundary_value]        
    def run_perturbation(self, boundaries, objectives = {}):
        
               
        if self.type == "cobra":
            if self.MOMA:
                self.run_perturbation_MOMA(boundaries, objectives)
            else:
                self.run_perturbation_FBA(boundaries, objectives)
            
    
    def run_perturbation_FBA(self, boundaries, objectives = {}):
      
        
        self.run_FBA()
        
        reaction_position = self.reaction_position
        
        perturbed_model = self.model.copy()
        
        """
        if self.pFBA:
            print("running pFBA")
            old_sol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(moma_model)                              
        else:
            print("running FBA")
            old_sol = moma_model.optimize()
        """        
        reaction_position = self.reaction_position
        # make perturbations        
        for (R, boundary_type, boundary_value) in boundaries:
            i = reaction_position[R]
            if boundary_type == 'l':
                perturbed_model.reactions[i].lower_bound = boundary_value
            elif boundary_type == 'u':                    
                perturbed_model.reactions[i].upper_bound = boundary_value

        if objectives:
            #model.objective = model.reactions[reaction_position[objective["objective_id"]]]
            # TODO
            perturbed_model.objective = {perturbed_model.reactions[self.reaction_position[r_id]]:int(objectives[r_id]) for r_id in objectives}
            #print(model.objective)     
                
        if self.pFBA:
            print("running pFBA")
            self.perturbed_sol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(perturbed_model)                              
        else:
            print("running FBA")
            self.perturbed_sol = perturbed_model.optimize()

        self.perturbed_fluxes = self.perturbed_sol.x_dict
        
        self.perturbed_model = perturbed_model
        
        """
        
        model = self.model
        reaction_position = self.reaction_position
        
        self.run_FBA()
        
        # make perturbations
        old_bounds = []
        for (R, boundary_type, boundary_value) in boundaries:
            i = reaction_position[R]
            if boundary_type == 'l':
                old_bounds.append(model.reactions[i].lower_bound)
                model.reactions[i].lower_bound = boundary_value
            elif boundary_type == 'u':
                old_bounds.append(model.reactions[i].upper_bound)
                model.reactions[i].upper_bound = boundary_value

        if objectives:
            old_objective = model.objective
            #model.objective = model.reactions[reaction_position[objective["objective_id"]]]
            model.objective = {model.reactions[self.reaction_position[r_id]]:int(objectives[r_id]) for r_id in objectives}
            #print(model.objective)               
            

       
        if self.pFBA:
            print("running pFBA")
            self.perturbed_sol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(model)                              
        else:
            print("running FBA")
            self.perturbed_sol = self.model.optimize()

        self.perturbed_fluxes = self.perturbed_sol.x_dict
        
        

    
        # recover previous state
        for ((R, boundary_type, boundary_value), bound) in zip(boundaries, old_bounds):
            i = reaction_position[R]
            if boundary_type == 'l':
                model.reactions[i].lower_bound = bound
            elif boundary_type == 'u':
                model.reactions[i].upper_bound = bound
        if objectives:
            model.objective = old_objective
    """
    
    #https://github.com/opencobra/cobratoolbox/blob/70b700944f07ea59227be1980e9f9ca65ecc7d7c/src/analysis/MOMA/MOMA.m
    def run_perturbation_MOMA(self, boundaries, objectives = {}):
        print("running MOMA...")
        
        moma_model = self.model.copy()
        
        """
        if self.pFBA:
            print("running pFBA")
            old_sol = cobra.flux_analysis.parsimonious.optimize_minimal_flux(moma_model)                              
        else:
            print("running FBA")
            old_sol = moma_model.optimize()
        """        
        reaction_position = self.reaction_position
        # make perturbations        
        for (R, boundary_type, boundary_value) in boundaries:
            i = reaction_position[R]
            if boundary_type == 'l':
                moma_model.reactions[i].lower_bound = boundary_value
            elif boundary_type == 'u':                    
                moma_model.reactions[i].upper_bound = boundary_value

        if objectives:
            #model.objective = model.reactions[reaction_position[objective["objective_id"]]]
            # TODO
            
            moma_model.objective = {moma_model.reactions[self.reaction_position[r_id]]:int(objectives[r_id]) for r_id in objectives}
            #print(model.objective)     
                   
    
        self.perturbed_sol = moma.moma(self.model, moma_model)
        self.perturbed_fluxes = self.perturbed_sol.x_dict
        self.perturbed_fluxes = {x:self.perturbed_fluxes[x] for x in self.perturbed_fluxes if x in self.fluxes}
    
        
        
        
    

    ######################################################################
    # function used in the review paper 10.1016/j.compbiomed.2017.07.005 #
    ######################################################################
    def objective_dependency(self, r1, r2, red_glc):    
        self.run_FBA()
        ref_flux1 = self.fluxes[r1]
        ref_flux2 = self.fluxes[r2]
       
        
        i1 = self.reaction_position[r1]
        i2 = self.reaction_position[r2]
        old_f1 = self.model.reactions[i1].lower_bound 
        old_f2 = self.model.reactions[i2].lower_bound
        
        #dh = 0.1                                    
        #reduction = [1 - i*dh for i in range(int(1/dh))]
        reduction = np.logspace(0,-5,num=6)
       
        
        sol = []
        
        red1 = red_glc
        f1 = ref_flux1 * red1 #red
        self.model.reactions[i1].lower_bound = f1
        # ali je smiselno vsakic prilagajati zgornjo mejo uptake-a glede na glukozo ali vedno uporabiti enako gornjo mejo?               
        #ref_flux2 = self.model.optimize().x_dict[r2]
        
        for red in reduction:
            
            f2 = ref_flux2 * red
        
            self.model.reactions[i2].lower_bound = f2
            print(red, f2)
            sol.append(self.model.optimize().f)
                           
        self.model.reactions[i1].lower_bound = old_f1
        self.model.reactions[i2].lower_bound = old_f2
                            
        sol_rel = []                            
        unperturbed = self.model.optimize().f
        for i in sol:
            sol_rel.append(i/unperturbed)
                                    
                            
        #plt.figure()
        #plt.plot(sol_rel)
        #plt.ylim(0,1)
        #plt.show()
        
        f = open(str(red1) + "_" + r1 + "_" + r2 + ".txt", 'w')
        f.write(r1 + ":" + str(red1) + "\n")
        f.write(r2 + ":" + str(reduction.tolist()) + "\n")
        f.write("relative solutions:" + str(sol_rel) + "\n")
        f.write("absolute solution:" + str(sol) + "\n")
        f.close()
        
    def save_pajek(self, filename, save_fluxes = 0):
        vertex_id = {}               
        f = open(filename + ".net", 'w')
        f_clu = open(filename + ".clu", 'w')
        f_vec = open(filename + ".vec", 'w')
        
        if save_fluxes:
            self.run_FBA()
        
        i = 1
        f.write("*Vertices " + str(len(self.reactions) + len(self.metabolites)) + "\n")
        f_clu.write("*Vertices " + str(len(self.reactions) + len(self.metabolites)) + "\n")
        if save_fluxes:
            f_vec.write("*Vertices " + str(len(self.reactions) + len(self.metabolites)) + "\n")
            
        for r in self.reactions:
            vertex_id[r] = i
            f.write(str(i) + " " + r + " " + "ellipse" + " " + "ic" + " " + "Green" + " " + "bc" + " " + "Green" + "\n")
            f_clu.write("1\n")
            i += 1
            if save_fluxes:
                f_vec.write(str(abs(self.fluxes[r])) + "\n")
        
        for m in self.metabolites:
            vertex_id[m] = i
            f.write(str(i) + " " + m + " " + "diamond" + " "+ "ic" + " " + "Red" + " " + "bc" + " " + "Red" + "\n")
            f_clu.write("2\n")
            #if save_fluxes:
            #    f_vec.write("1\n")
            i += 1
            
        f_clu.close()               
        f_vec.close()
        
        f.write("*Arcs\n")
     
            
        
        
        for r in self.reactions:
            id_r = vertex_id[r]
            
            reactants = self.reaction_reactants[r]
            for m in reactants:
                id_m = vertex_id[m]
                if not save_fluxes:
                    f.write(str(id_m) + " " + str(id_r)) 
                else:
                    fl = (self.fluxes[r])
                    if fl >= 0:
                        f.write(str(id_m) + " " + str(id_r) + " " + str(abs(fl)) + " c Green")                         
                    else:
                        f.write(str(id_r) + " " + str(id_m) + " " + str(abs(fl)) + " c Red") 
                    
                    #f.write(" " + str(abs(fl)))
                        
                f.write("\n")
            
            products = self.reaction_products[r]
            for m in products:
                id_m = vertex_id[m]
                if not save_fluxes:
                    f.write(str(id_r) + " " + str(id_m)) 
                else:
                    fl = (self.fluxes[r])
                    if fl >= 0:
                        f.write(str(id_r) + " " + str(id_m) + " " + str(abs(fl)) + " c Green")                         
                    else:
                        f.write(str(id_m) + " " + str(id_r) + " " + str(abs(fl)) + " c Red") 
                    
                    #f.write(" " + str(abs(fl)))
                        
                f.write("\n")
            
        
        f.close()
        
    """
    def get_ECs(self):
        self.kegg = KEGG()
        self.gene_EC = {}
        self.reaction_genes = defaultdict(set)
        
        for r in self.model.reactions:
            for g in r.genes:
                if g.id not in self.gene_EC:
                    self.gene_EC[g.id] = self.get_EC(g.id)
              
                #self.reaction_genes[r.id].add(self.gene_EC[g.id])
                self.reaction_genes[r.id] = self.reaction_genes[r.id] | (self.gene_EC[g.id])
        
    def get_ECs_from_file(self, file_name):
        self.load_ECs(file_name)          
        self.reaction_genes = defaultdict(set)
        
        for r in self.model.reactions:
            for g in r.genes:
                if g.id in self.gene_EC: 
                    self.reaction_genes[r.id] = self.reaction_genes[r.id] | (self.gene_EC[g.id])
               
    def get_EC(self, gene_id):
        # gene = kegg.parse(kegg.get("cge:100754434"))
        kegg = self.kegg
        gene_id = "cge:" + gene_id
        try:
            gene = kegg.parse(kegg.get(gene_id))
            for i in gene['ORTHOLOGY'].values():
                
                EC = i.split('[')[1].split(']')[0].replace("EC:","")
            return set(EC.split())
        except Exception as ex:
            template = "An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)
            print(gene_id)
            return set()
    
    
        
    def save_ECs(self, file_name):
        f = open(file_name, "w")
        for ge in self.gene_EC:
            f.write(ge+":")
            for ec in self.gene_EC[ge]:
                f.write(ec+" ")
            f.write("\n")
        f.close()
        
        
    def load_ECs(self, file_name):
        self.gene_EC = {}
        f = open(file_name, "r")
        for l in f:
            l = l.split(":")
            gene_id = l[0]
            ECs = set(l[1].strip().split(" "))
            ECs.discard('')
            self.gene_EC[gene_id] = ECs                            
        f.close()
    """