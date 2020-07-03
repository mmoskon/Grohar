from PyQt5 import uic

import model_data as md
#from bigg_data import bigg_data
from model_meta_data import model_meta_data
import pathway_data as pd
import make_pairs
from alignment import alignment
from paths import paths

from os.path import splitext
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtWidgets import QMessageBox
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.figure import Figure

from PyQt5.QtWidgets import QApplication

from bioservices.kegg import KEGG

from collections import defaultdict   

import cobra


#import requests

#import draw_plotly

import cobrababel

class top:
    def __init__(self):
       
        #self.app = QApplication([])
        
        self.dlg = uic.loadUi("top.ui")
        
        
        
        
        self.dlg.setWindowTitle("Grohar")
        
        
        
        
        self.dlg.cbModels.clear()
        self.dlg.cbModels.addItems(["gimmeS", "gimmeDG44", "gimmeK1"])
        self.dlg.btLoadModel.clicked.connect(self.load_model_fast)
        #
        #self.load_model()
        
        # BiGG Modesl
        self.BiGG_models = None
        self.dlg.cbModelsBiGG.clear()
        self.dlg.btLoadModelBiGG.clicked.connect(self.load_model_BiGG)
        
        
        
        self.dlg.btDraw.clicked.connect(self.draw_network)
        self.dlg.btExportNetwork.clicked.connect(self.export_network)
        self.dlg.btAddCompartment.clicked.connect(self.add_to_teCompartments)
        self.dlg.btAddMetabolite.clicked.connect(self.add_to_teMetabolites)
        #self.dlg.teMetabolites.setText("asn_L[c]\ngln_L[c]")
        self.dlg.teMetabolites.setText("lac_D[c]\nlac_D[e]\nlac_L[c]\nlac_L[e]\nlac_L[m]\nlac_L[x]")
       # self.dlg.leMinFlux
        
        
        self.dlg.cbCompartments.activated.connect(self.select_compartment)
        self.dlg.cbMetabolites.activated.connect(self.select_metabolite)
        self.dlg.cbReactions.activated.connect(self.select_reaction)
        self.dlg.cbFilterReactions.activated.connect(self.filter_reactions)
        
        
        
        self.dlg.btAddConstraint.clicked.connect(self.add_constraint)
        self.dlg.btAddObjective.clicked.connect(self.add_objective)
        self.dlg.btObjectiveCommit.clicked.connect(self.objective_commit)
        
    
        
        self.dlg.btDrawPerturbations.clicked.connect(self.draw_perturbations)
        self.dlg.btExportPerturbations.clicked.connect(self.export_perturbations)
        
        self.dlg.menuFileOpen.triggered.connect(self.openFileNameDialog)
        self.dlg.menuFileSaveAs.triggered.connect(self.openFileSaveAsDialog)
        
        
        #############
        #############
        # ALIGNMENT #
        #############
        #############
        
        
        self.dlg.btAlign.clicked.connect(self.align)
        self.dlg.btLoadBiGG.clicked.connect(self.load_bigg_data)
        self.dlg.btLoadKEGG.clicked.connect(self.load_model_kegg_data)
        self.dlg.btLoadKEGG.hide()
        
        self.dlg.btKEGG.clicked.connect(self.connect_KEGG)
        self.dlg.btKEGGOrganism.clicked.connect(self.select_organism_KEGG)
        #self.dlg.btKEGGPathway.clicked.connect(self.select_pathway_KEGG)
        
        self.dlg.cbPathway.activated.connect(self.select_pathway)
        
        self.dlg.cbSourceType.activated.connect(self.select_source_type)
        self.dlg.btLoadModelSrc.clicked.connect(self.openModelSrcFileNameDialog)
        
        self.dlg.gbSourcePathway.show()
        self.dlg.gbSourceModel.hide()
        
        
        self.dlg.btVisualiseAlignment.clicked.connect(self.visualise_alignment)
        self.dlg.btAddAlignmentReaction.clicked.connect(self.add_alignment_reaction)
        self.dlg.btSaveListReactions.clicked.connect(self.save_aligned_reactions)
        self.dlg.btLoadListReactions.clicked.connect(self.load_aligned_reactions)
        self.dlg.btExportAlignment.clicked.connect(self.export_aligned_reactions)
        
        
        #self.dlg.tabWidget.currentChanged.connect(self.tab_changed)
        
        ##################
    
        #self.dlg.figure = plt.figure()
        #self.dlg.canvas = FigureCanvas(self.dlg.figure)        
        #self.dlg.grid.addWidget(self.dlg.canvas)    
        
        #########
        #########
        # PATHS #
        #########
        #########
        
        self.dlg.btShortestPath.clicked.connect(self.shortest_path)
        self.dlg.btVisualiseSelectedShortestPath.clicked.connect(self.visualise_selected_path)
        self.dlg.btVisualiseAllShortestPaths.clicked.connect(self.visualise_all_paths)
        self.dlg.lePathMaxMetOccurrance.setText("50")
        
        self.dlg.btVisualiseSteiner.clicked.connect(self.visualise_path_multi)
        self.dlg.teSteiner.setText("EX_glc_e_\npyr[c]\nlac_L[c]")
        
        self.dlg.btDrawPerturbations_2.clicked.connect(self.shortest_paths_perturbation)

                
        ##############
        ##############
        # EDIT MODEL #
        ##############
        ##############
        
        
        self.dlg.btAddReaction.clicked.connect(self.add_reaction)
        self.dlg.cbEditReactions.activated.connect(self.select_edit_reaction)
        self.dlg.btDeleteReaction.clicked.connect(self.delete_reaction)
        self.dlg.btEditReactionsCommit.clicked.connect(self.edit_reactions_commit)
        
        self.dlg.btEditReactionsAddReactant.clicked.connect(self.add_reactant)
        self.dlg.btEditReactionsAddProduct.clicked.connect(self.add_product)
        self.dlg.btRemoveReactant.clicked.connect(self.remove_reactant)
        self.dlg.btRemoveProduct.clicked.connect(self.remove_product)
        
        self.dlg.btAddMetabolite_2.clicked.connect(self.add_metabolite)
        self.dlg.btDeleteMetabolite.clicked.connect(self.delete_metabolite)
        self.dlg.btEditMetaboliteCommit.clicked.connect(self.edit_metabolite)  
        self.dlg.cbEditMetabolitesMetabolites.activated.connect(self.select_edit_metabolite)
        
        
        self.dlg.btAddCompartment_2.clicked.connect(self.add_compartment)
        self.dlg.btDeleteCompartment.clicked.connect(self.delete_compartment)
        
        self.dlg.btAddObjective_2.clicked.connect(self.add_objective_2)
      
        
        ##################
        ##################
        # OTHER SETTINGS #
        ##################
        ##################
        self.dlg.actionPFBA.triggered.connect(self.check_pFBA)
        self.dlg.actionMOMA.triggered.connect(self.check_MOMA)
        self.dlg.actionFVA.triggered.connect(self.check_FVA)
          
        
        
        
        
    
        # varibales
        self.model = None            
        self.model_meta_data = None
        self.kegg = None
        self.alignment_loaded = False
        self.model_src = None
        
        self.alignment = None
        self.pathway_data = None
    
        
        
        ###############
        self.dlg.show()
   
        #self.app.exec()
        
    #def tab_changed(self):
    #    if self.dlg.tabWidget.currentIndex() == 1 and not self.alignment_initialised:
        
        
    def check_pFBA(self):             
        if self.model:
            self.model.pFBA = self.dlg.actionPFBA.isChecked()
            
            
    def check_MOMA(self):             
        if self.model:
            self.model.MOMA = self.dlg.actionMOMA.isChecked()
            
    def check_FVA(self):             
        if self.model:
            self.model.FVA = self.dlg.actionFVA.isChecked()
    
    
    def openFileNameDialog(self):   

        #options = QFileDialog.Options()
        #fileName = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        
        fname = QFileDialog.getOpenFileName(self.dlg, 'Open model', '', 'Model files (*.mat *.xml);; Matlab files (*.mat);; SBML files (*.xml);; All Files (*)')[0]
        
        if fname:
            self.model = md.model_data()
            model_name, ext = splitext(fname)
            self.dlg.setWindowTitle("Grohar - " + model_name + ext)
            if ext == ".mat":
                self.model.load_model_mat(fname)
            elif ext == ".xml":
                self.model.load_model_cobra_sbml(fname)                
            self.load_model()
            
            self.dlg.cbModels.clear()
            self.dlg.cbModels.setEnabled(False)
            self.dlg.btLoadModel.setEnabled(False)
            
            
    def openFileSaveAsDialog(self):    
        fname = QFileDialog.getSaveFileName(self.dlg, 'Save model', '', 'SBML files (*.xml);; Matlab files (*.mat);; All Files (*)')[0]
        
        if fname:
            model_name, ext = splitext(fname)
            self.dlg.setWindowTitle("Grohar - " + model_name + ext)
          
            if ext == ".mat":
                self.model.save_model_cobra_mat(fname)
            elif ext == ".xml":
                self.model.save_model_cobra_sbml(fname)       
    
    def filter_reactions(self):
        reactions = self.dlg.cbFilterReactions.currentText()
        
        self.dlg.cbReactions.clear()
        if reactions == 'Exchange':
            reaction_names = list(self.model.exchange_reactions)
        elif reactions == 'Uptake':
            reaction_names = list(self.model.uptake_reactions)
        else:
            reaction_names = list(self.model.reactions)
       
        reaction_names.sort()
        self.dlg.cbReactions.addItems(reaction_names)               
    
    def update_reaction_cbs(self):
        reaction_names = list(self.model.reactions)
        reaction_names.sort()
        self.dlg.cbReactions.clear()
        self.dlg.cbReactions.addItems(reaction_names)
        self.dlg.cbObjective.clear()
        self.dlg.cbObjective.addItems(reaction_names)
        
        self.dlg.cbEditReactions.clear()
        self.dlg.cbEditReactions.addItems(reaction_names)
        
        self.dlg.cbEditObjective.clear()
        self.dlg.cbEditObjective.addItems(reaction_names)
        
        
        
    
    
    def update_metabolite_cbs(self):
        met_names = list(self.model.metabolites)
        met_names.sort()
        self.dlg.cbMetabolites.clear()
        self.dlg.cbMetabolites.addItems(met_names)
        self.dlg.cbEditReactionsMetabolites.clear()
        self.dlg.cbEditReactionsMetabolites.addItems(met_names)
        self.dlg.cbEditMetabolitesMetabolites.clear()
        self.dlg.cbEditMetabolitesMetabolites.addItems(met_names)
    
    def update_compartments_cbs(self):
        compartments = list(self.model.compartments)
        compartments.sort()
        
        self.dlg.cbCompartments.clear()
        self.dlg.cbCompartments.addItems(compartments)
        
        self.dlg.cbEditMetabolitesCompartments.clear()
        self.dlg.cbEditMetabolitesCompartments.addItems(compartments)
        
        self.dlg.cbDeleteCompartments.clear()
        self.dlg.cbDeleteCompartments.addItems(compartments)
                
        
    
    def load_model(self):
        
        self.check_pFBA()
        self.check_FVA()
        self.check_MOMA()
        
                
        self.init_edit_reactions()
        
        self.update_compartments_cbs()
        
        
        self.update_metabolite_cbs()
        
        self.update_reaction_cbs()
        
        output = ""
        if self.model.objective:
            objectives = self.model.objective
            for i,r in enumerate(objectives):
                if i > 0:
                    output += " + "
                output += r.id + " * " + str(objectives[r])
            
            objective = list(self.model.objective.keys())[0].id
            index = self.dlg.cbObjective.findText(objective, QtCore.Qt.MatchFixedString)
            if index >= 0:
                self.dlg.cbObjective.setCurrentIndex(index)
                
        self.dlg.leObjective.setText(output)     
        # also in edit tab
        self.dlg.leEditObjective.setText(output)     
        
    
        index = self.dlg.cbMetabolites.findText("asn_L[c]", QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.dlg.cbMetabolites.setCurrentIndex(index)
         
   
            
        #index = self.dlg.cbReactions.findText("EX_o2_e_", QtCore.Qt.MatchFixedString)
        index = self.dlg.cbReactions.findText("EX_glc_e_", QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.dlg.cbReactions.setCurrentIndex(index)
            self.select_reaction()
            
        
     
        self.dlg.leMaxMetOccurance.setText("150")
        self.dlg.leMaxReactSize.setText("25")
        self.dlg.leMinFlux.setText("0.001")
        
        self.dlg.cbReactionTypes.clear()
        self.dlg.cbReactionTypes.addItems(["All", "Producing", "Consuming", "Both"])
        
        self.dlg.cbComparisonType.clear()
        self.dlg.cbComparisonType.addItems([ "Modified reactions", "Added reactions", "Removed reactions", "Remained reactions"])
       
     
         
        
        
        self.model.run_FBA()
        
        #if not model.sol:      
        #    QMessageBox.warning(self.dlg, "Warning","Solution not found!")
           
        
        self.init_alignment_tab()
    
        
        """
        GLC = 'EX_glc_e_'
        AA = ['EX_asn_L_e_',         
              'EX_gln_L_e_',
              'EX_asp_L_e_',
              'EX_arg_L_e_',
              'EX_cys_L_e_',
              'EX_his_L_e_',
              'EX_ile_L_e_',
              'EX_leu_L_e_',
              'EX_lys_L_e_',
              'EX_met_L_e_',
              'EX_phe_L_e_',
              'EX_pro_L_e_',
              'EX_ser_L_e_',
              'EX_thr_L_e_',
              'EX_trp_L_e_',
              'EX_tyr_L_e_',
              'EX_val_L_e_']
        
        for i in AA:
            self.model.objective_dependency(GLC, i, 1)
            self.model.objective_dependency(GLC, i, 0.1)
            self.model.objective_dependency(GLC, i, 0.01)
        """
        
        #self.select_reaction() 
        self.alignment = None
        self.pathway_data = None
    
    def load_model_fast(self):
        self.model = md.model_data()
        #self.model.load_model_sbml("iCHOv1_final.xml")
        #self.model.load_model_cobra_sbml("iCHOv1_final.xml")
        #self.model.load_model_mat("gimmeDG44.mat")
        #self.model.load_model_mat("gimmeS.mat")
        #self.model.load_model_mat("gimmeK1.mat")
             
        
        model_name = "".join([self.dlg.cbModels.currentText(), ".mat"])
        self.dlg.setWindowTitle("Grohar - " + model_name)
        
        self.model.load_model_mat(model_name)
        self.load_model()
        
        
    def load_model_BiGG(self):
        if not self.dlg.cbModelsBiGG.currentText():
            model_ids = [i['bigg_id'] for i in cobrababel.get_bigg_model_list()]
            model_ids.sort()
            self.dlg.cbModelsBiGG.addItems(model_ids)
            self.BiGG_models = model_ids
        else:
            model_id = self.dlg.cbModelsBiGG.currentText()
            self.dlg.setWindowTitle("Grohar - " + model_id)
            model = cobrababel.create_cobra_model_from_bigg_model(model_id)
            self.model = md.model_data()
            self.model.import_model(model)
            self.load_model()    
            
    def add_to_teCompartments(self):
        compartments = set(self.dlg.teCompartments.toPlainText().strip().split("\n")) - {""}
        compartments.add(self.dlg.cbCompartments.currentText())
        compartments = list(compartments)
        compartments.sort()
        compartments_disp = "\n".join(compartments)
        self.dlg.teCompartments.setText(compartments_disp)
  
        # zaenkrat dela samo dodajanje!
        
        
        self.dlg.cbMetabolites.clear()
        
        met_names = set(i for c in compartments for i in self.model.compartment_metabolites[c])      
        met_names = list(met_names)
        met_names.sort()
        self.dlg.cbMetabolites.addItems(met_names)
        
        index = self.dlg.cbMetabolites.findText("asn_L[c]", QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.dlg.cbMetabolites.setCurrentIndex(index)
    
        
    def add_to_teMetabolites(self):
        met_names = set(self.dlg.teMetabolites.toPlainText().strip().split("\n"))
        met_names.add(self.dlg.cbMetabolites.currentText())
        met_names = list(met_names)
        met_names.sort()
        met_names = "\n".join(met_names)
        self.dlg.teMetabolites.setText(met_names)
     
    def select_compartment(self):
        model = self.model
        compartment = self.dlg.cbCompartments.currentText()
        self.dlg.leDebug.setText(model.compartments[compartment])

          
    def select_metabolite(self):
        model = self.model
        metabolite = self.dlg.cbMetabolites.currentText()
        self.dlg.leDebug.setText(model.metabolites[metabolite])
        
    def select_reaction(self):
        reaction = self.dlg.cbReactions.currentText()
        model = self.model
        
        #i = model.reaction_position[reaction]
        #self.dlg.leLowerBound.setText(str(model.model.reactions[i].lower_bound))
        #self.dlg.leUpperBound.setText(str(model.model.reactions[i].upper_bound))
        self.dlg.leLowerBound.setText(str(model.reaction_bounds[reaction][0]))
        self.dlg.leUpperBound.setText(str(model.reaction_bounds[reaction][1]))
        
        
        if self.model.fluxes and reaction in self.model.fluxes:
            self.dlg.leDebug.setText(model.reactions[reaction] + ": "  + str(self.model.fluxes[reaction]))
        else:
            self.dlg.leDebug.setText(model.reactions[reaction])
           
        
    def get_constraints(self):
        constraints = self.dlg.teConstraints.toPlainText().split("\n")
        reactions = []
        lower = []
        upper = []


        # (constraints)
        if constraints[0]: 
            for constr in constraints:
                if (constr.strip()):
                    c = constr.split()
                    reactions.append(c[0])
                    lower.append(float(c[1]))
                    upper.append(float(c[2]))

        #print(reactions)

                   
        return (reactions, lower, upper)
        
        
    def add_constraint(self):
        (reactions, lower, upper) = self.get_constraints()
        reaction_name = self.dlg.cbReactions.currentText()
      
        bounds = (float(self.dlg.leLowerBound.text()), float(self.dlg.leUpperBound.text()))
        if reaction_name in reactions:
            lower[reactions.index(reaction_name)] = bounds[0]
            upper[reactions.index(reaction_name)] = bounds[1]
        else:
            reactions.append(reaction_name)
            lower.append(bounds[0])
            upper.append(bounds[1])
         
        self.dlg.teConstraints.clear()
        tt = ""        
        for (r, l, u) in zip(reactions, lower, upper):
            tt = tt + r + "\t" + str(l) + "\t" + str(u) + "\n"
             
        print(tt)
        self.dlg.teConstraints.setText(tt)
    
    def parse_objective(self, option = 0):
        objectives = defaultdict(int)
        if option == 2:
            inpt = self.dlg.leEditObjective.text()
        else:
            inpt = self.dlg.leObjective.text()
        if "*" in inpt:
            inpt = inpt.split(" + ")
            for x in inpt:
                i, v = x.split(" * ")
                v = float(v)
                objectives[i] = v
        return objectives
    
    def add_objective(self):
        try:   
            objectives = self.parse_objective()               
            r_id  = self.dlg.cbObjective.currentText()
            objectives[r_id] += 1
                
            output = ""
            for i,r_id in enumerate(objectives):
                if i > 0:
                    output += " + "
                output += r_id + " * " + str(objectives[r_id])
            self.dlg.leObjective.setText(output)                                 
        except Exception as inst:
            print(inst)
            
            
    def get_subnetwork_params(self):
        if self.model:
            model = self.model
        else:
            return
        
        model.run_FBA()
        
        
              
        
        met_names = set(self.dlg.teMetabolites.toPlainText().split("\n"))
        
     
        reaction_types = 2 if self.dlg.cbReactionTypes.currentText() == "Both" else 1 if self.dlg.cbReactionTypes.currentText() == "Producing" else -1 if self.dlg.cbReactionTypes.currentText() == "Consuming" else 0
        distance = int(self.dlg.sbDistance.value())
        
      
        if self.dlg.leMinFlux.text():
            min_flux = float(self.dlg.leMinFlux.text())
        else:
            min_flux = 0
            
       
        if self.dlg.leMaxReactSize.text():
            max_react_size = int(self.dlg.leMaxReactSize.text())
        else:
            max_react_size  = 0
            
      
        if self.dlg.leMaxMetOccurance.text():
            max_metabolite_occurance = int(self.dlg.leMaxMetOccurance.text())
        else:
            max_metabolite_occurance = 0 
        
        
        compartments = set(self.dlg.teCompartments.toPlainText().strip().split("\n")) - {""}
        compartments = list(compartments)
        compartments.sort()
        compartments = " ".join(compartments)
  
        return (reaction_types, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux)
    
    
    def draw_network(self):
        if self.model:
            model = self.model
        else:
            return
        
        
        if not model.sol:      
            QMessageBox.warning(self.dlg, "Warning","Solution not found!")
            return
        
        
        (reaction_types, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux) = self.get_subnetwork_params()
                               
        (pairs, fluxes) = make_pairs.make_pairs(reaction_types, model, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux = min_flux)
      
      
        
        self.draw(met_names, pairs, fluxes)
   
    def export_network(self):
        if self.model:
            model = self.model
        else:
            return
        
        (reaction_types, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux) = self.get_subnetwork_params()
                               
        return_pairs = False                
        reactions, fluxes = make_pairs.make_pairs(reaction_types, model, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux = min_flux, return_pairs = return_pairs)
        
        self.export_reactions(reactions, fluxes)
        
        
        
    # reactions: set of reaction IDs
    # fluxes: dictionary with reaction IDs as keys        
    def export_reactions(self, reactions, fluxes = []):
        if self.model:
            model = self.model
        else:
            return
        
        
        if not fluxes:
            fluxes = model.fluxes
        
        
        submodel = cobra.Model(model.model.id)
        subfluxes = {}
        
        for r in reactions:
            r_idx = model.reaction_position[r]
            reaction = model.model.reactions[r_idx]
            
            submodel.add_reactions([reaction])
            subfluxes[r] = fluxes[r]
        
        cobra.io.write_sbml_model(submodel, "submodel.xml")
        #
        #cobra.io.to_json(submodel)
        #print(subfluxes)
        #
        
        """submodel = model.model.copy()
        reactions_to_remove = []
        metabolites_to_keep = set()
        
        print(1)
        
        for r in submodel.reactions:            
            metabolites_to_keep |= set(r.metabolites.keys())
            if r.id not in reactions:
                reactions_to_remove.append(r)

        print(2)        
        submodel.remove_reactions(reactions_to_remove)
        
        print(3)
        meta = cobra.core.DictList(metabolites_to_keep)
        print(4)
        submodel.metabolites = meta
        print(5)
        
        
        #metabolites_to_remove = []
        #for m in submodel.metabolites:
        #    if m not in metabolites_to_keep:
        #        metabolites_to_remove.append(m)
        
        
        
        
        
        
        #cobra.io.save_json_model(submodel, "test.json")
        
        cobra.io.write_sbml_model(submodel, "test2.xml")
        
        print(6)
        """
        
        
        """ TODO
            iti čez reakcije v reactions in to shraniti v model       
            
            fluxes: lahko so v obliki
            - id_reakcije: flux
            - id_reakcije: (flux_pred, flux_po)
        """
        
        
    def draw_reactions(self, reacts, max_metabolite_occurance = 0, observed = None):
        if self.model:
            model = self.model
        else:
            return
        
           
        model.run_FBA()
        
        
        if not model.sol:      
            QMessageBox.warning(self.dlg, "Warning","Solution not found!")
            return
        
        #if self.dlg.leMaxMetOccurance.text():
        #    max_metabolite_occurance = int(self.dlg.leMaxMetOccurance.text())
        #else:
        #    max_metabolite_occurance = 0 
        
        #print(reacts)
        
        (pairs, fluxes) = make_pairs.make_pairs_reactions(model, reacts, max_metabolite_occurance)
        
        #print(pairs)
        
        if observed:
            met_names = observed
        else:
            met_names = ""
        
        self.draw(met_names, pairs, fluxes, path = reacts)
   
    def get_perturbation_params(self):
        if self.model:
            model = self.model
        else:
            return     
    
    
        (reactions, lower, upper) = self.get_constraints()  
        boundaries = []
        for (r, l, u) in zip(reactions, lower, upper):
            boundaries.append([r, 'l', l]) 
            boundaries.append([r, 'u', u]) 
        
        #objective_id = self.dlg.cbObjective.currentText() #"biomass_cho_producing"
        objectives = self.parse_objective()
                    
        #model.run_FBA()
        model.run_perturbation(boundaries, objectives = objectives)            
        
              
        if not model.perturbed_sol:      
            QMessageBox.warning(self.dlg, "Warning","Solution not found!")
            return
            
        print("Original solution {:10.10f}".format(model.sol.f))
        print("Perturbed solution {:10.10f}".format(model.perturbed_sol.f))
        
        
        #comparison_type = 1 # new reactions = 1, reaction removed = -1, reactions in both model= 0
        #comparison_type = 1 if self.dlg.cbComparisonType.currentText() == "Added reactions" else -1 if self.dlg.cbComparisonType.currentText() == "Removed reactions" else 2 if self.dlg.cbComparisonType.currentText() == "Modified reactions" else 0        
        comparison_type = 1 if self.dlg.cbComparisonType.currentText() == "Added reactions" else -1 if self.dlg.cbComparisonType.currentText() == "Removed reactions" else 3 if self.dlg.cbComparisonType.currentText() == "Modified reactions" else 0        
        
        
        
        met_names = set(self.dlg.teMetabolites.toPlainText().split("\n"))
        
        reaction_types = 2 if self.dlg.cbReactionTypes.currentText() == "Both" else 1 if self.dlg.cbReactionTypes.currentText() == "Producing" else -1 if self.dlg.cbReactionTypes.currentText() == "Consuming" else 0
        
        distance = int(self.dlg.sbDistance.value())
        if self.dlg.leMinFlux.text():
            min_flux = float(self.dlg.leMinFlux.text())
        else:
            min_flux = 0
        if self.dlg.leMaxReactSize.text():
            max_react_size = int(self.dlg.leMaxReactSize.text())
        else:
            max_react_size  = 0
        if self.dlg.leMaxMetOccurance.text():
            max_metabolite_occurance = int(self.dlg.leMaxMetOccurance.text())
        else:
            max_metabolite_occurance = 0 
       
            
        compartments = set(self.dlg.teCompartments.toPlainText().strip().split("\n")) - {""}
        compartments = list(compartments)
        compartments.sort()
        compartments = " ".join(compartments)
    
        return (comparison_type, reaction_types, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux)
    

    def export_perturbations(self):
        if self.model:
            model = self.model
        else:
            return 

        (comparison_type, reaction_types, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux) = self.get_perturbation_params()
        
        return_pairs = False
        (pairs, reactions) = make_pairs.make_pairs_perturbation(comparison_type, reaction_types, model, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux = min_flux, return_pairs = return_pairs)

        self.export_reactions(pairs, reactions)


    def draw_perturbations(self):
        if self.model:
            model = self.model
        else:
            return     
            
        (comparison_type, reaction_types, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux) = self.get_perturbation_params()
      
        (pairs, fluxes) = make_pairs.make_pairs_perturbation(comparison_type, reaction_types, model, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux = min_flux)
        #(pairs, fluxes) = make_pairs.make_pairs_perturbation(comparison_type, reaction_types, model, met_names, distance, compartments, max_react_size, max_metabolite_occurance, min_flux = min_flux)
        
        self.draw(met_names, pairs, fluxes)
        
       
    def on_plot_hover(self, event):
        pos = self.pos
        #annotation = self.annotation
        label = self.label
        
        nodes = [i for i in pos if abs(pos[i][0] - event.xdata) < 50 and abs(pos[i][1] - event.ydata) < 50]
                  
        if nodes:             
            #labels = [a + ":" + "{:.5f}".format(fluxes[a]) if a in self.model.reactions else a for a in nodes]
            if self.model.FVA:
                labels = self.label_FVA(nodes)
            else:
                labels = [self.model.reactions[a] if a in self.model.reactions else self.model.metabolites[a] for a in nodes]
            label.set_text("\n".join(labels))
            label.set_position((event.xdata, event.ydata))
            
            bb = label.get_bbox_patch()
            bb.set_boxstyle("round", pad=0.6)
            label.set_visible(True)
        else:
            label.set_visible(False)
       
        plt.draw()
       
    def label_FVA(self, nodes):
        labels = []
        for n in nodes:
            if n in self.model.metabolites:
                labels.append(self.model.metabolites[n])
            elif n in self.model.reactions:
                if n not in self.FVA_dict:
                    self.FVA_dict[n] = self.model.run_FVA(n)                                        
                fva = "[" + str(round(self.FVA_dict[n][0],3)) + ", " + str(round(self.FVA_dict[n][1],3)) + "]"
                
                
                labels.append(self.model.reactions[n] + " " + fva)                                    
            else:
                labels.append(n)
        return labels
            
                
                


    def on_key_press(self, event):
        if event.key == "l":
                self.lgd_visible = not self.lgd_visible
                self.lgd.set_visible(self.lgd_visible)
                plt.draw()
               
        
    def on_plot_click(self, event):
        pos = self.pos
        nodes = [i for i in pos if abs(pos[i][0] - event.xdata) < 50 and abs(pos[i][1] - event.ydata) < 50]
        
        if event.button == 1:   # move the selected node
            if nodes:
                self.selected = nodes[0]
                
            else:
                self.selected = None
        elif event.button == 3: # delete the selected node
            self.selected = None   
            node_name = nodes[0]
            node_id = self.node_ids[nodes[0]]
            
            self.G.remove_node(node_name)
            del self.node_labels[node_name]                       
            del self.node_colors[node_id]
            del self.node_sizes[node_id]
            for i in self.node_ids:
                if self.node_ids[i] > node_id:
                    self.node_ids[i] -= 1
                                
            self.redraw()
            
            #self.ax.clear()
            #nx.draw_networkx_labels(self.G, self.pos , self.node_labels, font_size=8, font_weight = "bold", font_color = "#1c2833", ax = self.ax)
            #nx.draw(self.G, self.pos, node_color = self.node_colors, node_size=self.node_sizes, edge_color = "#abb2b9", dge_cmap=matplotlib.cm.Reds, alpha = 0.5, ax = self.ax)
            
        
    def on_plot_click_release(self, event):  
       
        if event.button == 1:   # move the selected node     
            if self.selected:
                self.pos[self.selected] = (event.xdata, event.ydata)
                
                self.redraw()
                
                #self.ax.clear()
            
                       
                #nx.draw_networkx_labels(self.G, self.pos , self.node_labels, font_size=8, font_weight = "bold", font_color = "#1c2833", ax = self.ax)
                #nx.draw(self.G, self.pos, node_color = self.node_colors, node_size=self.node_sizes, edge_color = "#abb2b9", dge_cmap=matplotlib.cm.Reds, alpha = 0.5, ax = self.ax)
                

    
    def redraw(self):
        self.ax.clear()
        
        bbox_props = dict(boxstyle="round,pad=0.3", ec="0.5", lw=2)
        self.label = self.ax.text(100, 100, "text", color = "white", bbox=bbox_props, zorder=10)
        self.label.set_visible(False)
           
        
        
        if self.dlg.actionCBFriendly.isChecked():
            
            patches = [mpatches.Patch(color='#f0e442', label='observed metabolites/path'),                    
                       mpatches.Patch(color='#0072b2', label='metabolites'),
                       mpatches.Patch(color='#009e73', label='producing reactions'),
                       mpatches.Patch(color='#d55e00', label='consuming reactions'),                                 
                       mpatches.Patch(color='#7570b3', label='producing and consuming reactions'),                                 
                       mpatches.Patch(color='#cc79a7', label='other reactions')]
        
        else:            
        
            patches = [mpatches.Patch(color='m', label='observed metabolites/path'),                    
                       mpatches.Patch(color='c', label='metabolites'),
                       mpatches.Patch(color='g', label='producing reactions'),
                       mpatches.Patch(color='r', label='consuming reactions'),                                 
                       mpatches.Patch(color='y', label='producing and consuming reactions'),                                 
                       mpatches.Patch(color='b', label='other reactions')]
        
        
        self.lgd = plt.legend(handles = patches)
        self.lgd.set_visible(self.lgd_visible)
        
        
        nx.draw_networkx_labels(self.G, self.pos , self.node_labels, font_size=8, font_weight = "bold", font_color = "#1c2833", ax = self.ax)
        nx.draw(self.G, self.pos, node_color = self.node_colors, node_size=self.node_sizes, edge_color = "#abb2b9", dge_cmap=matplotlib.cm.Reds, alpha = 0.5, ax = self.ax)
        
        
        #self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)  
        #self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        #self.fig.canvas.mpl_connect('button_press_event', self.on_plot_click) 
        #self.fig.canvas.mpl_connect('button_release_event', self.on_plot_click_release) 
    
    
   
    
    
    def draw(self, met_names, pairs, fluxes, path = None):        
            
        #draw_plotly.draw_plotly(self.model, met_names, pairs, fluxes)
        
        model = self.model
                
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
       
        self.FVA_dict = {}
        
        #self.annotation = self.ax.annotate("Mouseover point",
        #xy=(100, 100), xycoords='data',
        #xytext=(100 + 1, 100), textcoords='data',
        #horizontalalignment="left",
        #arrowprops=dict(arrowstyle="simple",
        #                connectionstyle="arc3,rad=-0.2"),
        #bbox=dict(boxstyle="round", facecolor="w", 
        #          edgecolor="0.5", alpha=0.9)
        #)
        
        #self.annotation.set_visible(False)
        
        bbox_props = dict(boxstyle="round,pad=0.3", ec="0.5", lw=2)
        self.label = self.ax.text(100, 100, "text", color = "white", bbox=bbox_props, zorder=10)
        self.label.set_visible(False)
        
        
     
        pairs_mod = []
        for i,(a,b) in enumerate(pairs):
            pairs_mod.append((a.replace(":","---"), b.replace(":","---"))) # networkx uses ":" as a special symbol!
      
            
        #G = nx.DiGraph(splines = True)
        G = nx.DiGraph()
        
        #nodes = list(set(i for pair in pairs_mod for i in pair))
        #G.add_nodes_from(nodes)
        
        
        
        G.add_edges_from(pairs_mod)
        
        if fluxes and type(list(fluxes.values())[0]) == tuple:           
            fl = {}
            for i in fluxes:
                fl[i] = fluxes[i][1]
            producing = set(make_pairs.direct_producers(model, met_names, fl))
            consuming = set(make_pairs.direct_consumers(model, met_names, fl))            
        else:
            producing = set(make_pairs.direct_producers(model, met_names, fluxes))
            consuming = set(make_pairs.direct_consumers(model, met_names, fluxes))
            
        
        # https://www.color-hex.com/color-palette/49436
        
        
        
        #print(self.dlg.actionCBFriendly.isChecked())
        
        if self.dlg.actionCBFriendly.isChecked():
            
            val_map = {'reaction': '#cc79a7',
                       'metabolite': '#0072b2',
                       'producing': '#009e73',
                       'both': '#7570b3',
                       'consuming': '#d55e00',
                       'observed': '#f0e442',
                       'path': '#7570b3'}
            
            patches = [mpatches.Patch(color='#f0e442', label='observed metabolites/path'),                    
                       mpatches.Patch(color='#0072b2', label='metabolites'),
                       mpatches.Patch(color='#009e73', label='producing reactions'),
                       mpatches.Patch(color='#d55e00', label='consuming reactions'),                                 
                       mpatches.Patch(color='#7570b3', label='producing and consuming reactions'),                                 
                       mpatches.Patch(color='#cc79a7', label='other reactions')]
        
        else:
            
            val_map = {'reaction': 'b',
                       'metabolite': 'c',
                       'producing': 'g',
                       'both': 'y',
                       'consuming': 'r',
                       'observed': 'm',
                       'path': 'm'}
            
            patches = [mpatches.Patch(color='m', label='observed metabolites/path'),                    
                       mpatches.Patch(color='c', label='metabolites'),
                       mpatches.Patch(color='g', label='producing reactions'),
                       mpatches.Patch(color='r', label='consuming reactions'),                                 
                       mpatches.Patch(color='y', label='producing and consuming reactions'),                                 
                       mpatches.Patch(color='b', label='other reactions')]
    
        
        
        self.lgd = plt.legend(handles = patches)
        self.lgd_visible = False
        self.lgd.set_visible(self.lgd_visible)
        
        
        #node_colors = [val_map['observed'] if node in met_names else val_map['both'] if node in producing & consuming else val_map['producing'] if node in producing else val_map['consuming'] if node in consuming else val_map['reaction'] if node in model.reactions else val_map['metabolite'] for node in G.nodes()]
        node_colors = []
        for node in G.nodes():
          
            node_orig = node.replace("---",":")
            
            if node_orig in met_names or (path and node_orig in path):
                color = val_map['observed']
            elif node_orig in producing & consuming:
                color = val_map['both']
            elif node_orig in producing:
                color = val_map['producing']
            elif node_orig in consuming:
                color = val_map['consuming']
            elif node_orig in model.reactions:
                color = val_map['reaction'] 
            else:
                color = val_map['metabolite']             
            node_colors.append(color)            
        

        self.node_ids = {node:i for i, node in enumerate(G.nodes())}

        #node_labels = {a:a + ":" + "{:.5f}".format(fluxes[a]) if a in model.reactions else a for a in G.nodes()}
        node_labels = {}
        for node in G.nodes():                             
            node_orig = node.replace("---",":")
            if node_orig in model.reactions:
                f = fluxes[node_orig]
                if type(f) == tuple:
                    node_labels[node] = node_orig + " :" + "{:.3f}\n{:.3f}".format(round(f[0],3),round(f[1],3))
                else:
                    node_labels[node] = node_orig + " :" + "{:.3f}".format(round(fluxes[node_orig],3))
            else:
                node_labels[node] = node_orig
        
        #node_sizes = [500 if a in model.reactions else 250 for a in G.nodes()]
        node_sizes = []
        for node in G.nodes():
            node_orig = node.replace("---",":")
            if node_orig in model.reactions:
                node_sizes.append(500)
            else:
                node_sizes.append(250)
        
        
        
        force_layout = 0             
     
        if force_layout:
        
            if self.pathway_data.reaction_coordinates:         
                react_coords = {}
               
                
                x_coords = [i for i,_ in self.pathway_data.reaction_coordinates.values()]
                y_coords = [j for _,j in self.pathway_data.reaction_coordinates.values()]
                #min_x = min(x_coords)
                max_x = max(x_coords)
                #min_y = min(y_coords)
                max_y = max(y_coords)
                scale_x = 1800/max_x
                scale_y = 1200/max_y
         
            
                for i in self.alignment.last_alignment:
                    src_id = self.alignment.reactions_src[i]
                    dst_id = self.alignment.reactions_dst[self.alignment.last_alignment[i]]
                    
                    pos_x = float(self.pathway_data.reaction_coordinates[src_id][0] * scale_x)# + 1800 * min_x/max_x)
                    pos_y = float(self.pathway_data.reaction_coordinates[src_id][1] * scale_y)# + 1200 * min_y/max_y)
                    
                    react_coords[dst_id] = (pos_x, pos_y)
                    
                
                for node in G.nodes():                                             
                    node_orig = node.replace("---",":")
                    if node_orig in react_coords:
                        #G.node[node]['pos'] = react_coords[node_orig]
                        #G.node[node]['pin'] = True
                        
                        
                        G.node[node]['pos'] = str(react_coords[node_orig]).replace("(","").replace(")","").replace(" ","")+"!"
                        #print(node)
                        #print(G.node[node]['pos'])
               
                
            
            
            #pos = nx.nx_pydot.graphviz_layout(G, prog="twopi")              
            #pos = nx.nx_pydot.graphviz_layout(G, prog="circo")              
            #pos = nx.nx_pydot.graphviz_layout(G, prog="sfdp")              
            pos = nx.nx_pydot.graphviz_layout(G, prog="fdp")
            
            x_coords = [i for i,_ in pos.values()]
            y_coords = [j for _,j in pos.values()]
            #min_x = min(x_coords)
            max_x = max(x_coords)
            #min_y = min(y_coords)
            max_y = max(y_coords)
            scale_x = 1800/max_x
            scale_y = 1200/max_y
         
            for i in pos:
                pos[i] = (pos[i][0]*scale_x, pos[i][1]*scale_y)
            #pos = nx.spring_layout(G)
            #pos = nx.spectral_layout(G)
            #plt.clf()
            
            # set the positions!
            """
            print(G.node)
            print("****************")
            print("****************")
            print("****************")
            print("****************")
            print(pos)
            print("****************")
            print("****************")
            print("****************")
            print("****************")
            print(nx.get_node_attributes(G,'pos'))
            """
        
        else:
            pos = nx.nx_pydot.graphviz_layout(G, prog="dot")
        
        
    
    
        nx.draw_networkx_labels(G, pos , node_labels, font_size=8, font_weight = "bold", font_color = "#1c2833", ax = self.ax)
        #nx.draw_networkx_labels(G,labels=node_labels)
        
       
        nx.draw(G, pos, node_color = node_colors, node_size=node_sizes, edge_color = "#abb2b9", dge_cmap=matplotlib.cm.Reds, alpha = 0.5, ax = self.ax)
      
        
        ##nx.draw_networkx_nodes(G, pos, nodelist = [a for a in G if a in model.reactions], node_color = node_colors, node_size=250,  alpha=0.5)
        ##nx.draw_networkx_nodes(G, pos, nodelist = [a for a in G if a in model.metabolites], node_color = node_colors, node_size=100,  alpha=0.5)
        ##nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
        
        ##pylab.show()
        
        
        #self.dlg.canvas.draw_idle()
        
        #handles, labels = self.ax.get_legend_handles_labels()
        
        
        self.G = G
        self.pos = pos
        #self.pos = nx.get_node_attributes(G,'pos')
        self.node_labels = node_labels
        self.node_colors = node_colors
        self.node_sizes = node_sizes
        
      
        
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)  
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.fig.canvas.mpl_connect('button_press_event', self.on_plot_click) 
        self.fig.canvas.mpl_connect('button_release_event', self.on_plot_click_release) 
    
       
                 
        #nodes = G.nodes()
        #edges = G.edges()
        #nodes_dict = [{"id":n} for n in nodes]
        #node_map = dict(zip(nodes,range(len(nodes)))) # map to indices for source/target in edges
        #edges_dict = [{"source":node_map[edges[i][0]], "target":node_map[edges[i][1]],"title":'test'} for i in range(len(edges))]
        #visJS.visjs_network(nodes_dict, edges_dict, time_stamp=0)
        
    def init_alignment_tab(self):
        
        if not self.alignment_loaded:
        
            self.dlg.cbSourceType.addItems(['KEGG pathway', 'MAT or SBML'])
            
            self.dlg.cbBiGGModels.addItem('None')
            
            try:
                if not self.BiGG_models:
                    model_ids = [i['bigg_id'] for i in cobrababel.get_bigg_model_list()]
                    model_ids.sort()
                    self.dlg.cbModelsBiGG.addItems(model_ids)
                    self.BiGG_models = model_ids    
                                        
                self.dlg.cbBiGGModels.addItems(self.BiGG_models)
                """        
                bigg_models = requests.get('http://bigg.ucsd.edu/api/v2/models/').json()
                for b in bigg_models['results']:
                    self.dlg.cbBiGGModels.addItem(b['bigg_id'])
                """
            except:
                pass
        
        #print(self.model.model.id)    
        
        index = self.dlg.cbBiGGModels.findText(self.model.model.id, QtCore.Qt.MatchFixedString)
        #print(index)
        if index >= 0:
            self.dlg.cbBiGGModels.setCurrentIndex(index)
        else:
            self.dlg.cbBiGGModels.setCurrentIndex(0)
        
        
        
        if self.model.model.id == "e_coli_core":
            self.dlg.lePathwayId.setText("eco00020")
        if self.model.model.id == "test_model":
            self.dlg.lePathwayId.setText("test_pathway")
    
        #self.dlg.leBiGGID.setText("iCHOv1")#self.model.bigg_id)
          
        self.dlg.cbAlignmentCompartments.clear()
        compartments = list(self.model.compartments)
        compartments.sort()
        self.dlg.cbAlignmentCompartments.addItem("None")
        self.dlg.cbAlignmentCompartments.addItems(compartments)
        index = self.dlg.cbAlignmentCompartments.findText('c', QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.dlg.cbAlignmentCompartments.setCurrentIndex(index)
            
            
        if not self.alignment_loaded:
            
            self.dlg.leAlignmentMaxMetOccurance.setText("10")
                
            self.dlg.cbReversibility.clear()
            self.dlg.cbReversibility.addItems(["Yes","No"])
            self.dlg.cbReversibility.setCurrentIndex(0)
            
            self.dlg.cbKm.clear()
            self.dlg.cbKm.addItems(list(map(str,range(1,11))))
            self.dlg.cbKm.setCurrentIndex(1)
            
            self.dlg.cb1ToMany.clear()
            self.dlg.cb1ToMany.addItems(["Yes","No"])
            self.dlg.cb1ToMany.setCurrentIndex(1)
            
            self.dlg.cbManyTo1.clear()
            self.dlg.cbManyTo1.addItems(["Yes","No"])
            self.dlg.cbManyTo1.setCurrentIndex(1)
            
            self.dlg.cbRemoveBorders.clear()
            self.dlg.cbRemoveBorders.addItems(["Yes","No"])
            self.dlg.cbRemoveBorders.setCurrentIndex(0)
            
            self.dlg.cbTopologicalWeight.clear()
            self.dlg.cbTopologicalWeight.addItems(list(map(str,range(0,11))))
            self.dlg.cbTopologicalWeight.setCurrentIndex(1)
            
            self.dlg.cbECWeight.clear()
            self.dlg.cbECWeight.addItems(list(map(str,range(0,11))))
            self.dlg.cbECWeight.setCurrentIndex(1)
            
            
            self.dlg.cbKEGGWeight.clear()
            self.dlg.cbKEGGWeight.addItems(list(map(str,range(0,11))))
            self.dlg.cbKEGGWeight.setCurrentIndex(1)
            
            self.dlg.cbAdditive.clear()
            self.dlg.cbAdditive.addItems(["Additive", "Multiplicative"])
            self.dlg.cbAdditive.setCurrentIndex(1)
            
            self.dlg.cbAlignmentType.clear()
            self.dlg.cbAlignmentType.addItems(["Simple", "Greedy", "Probabilistic"])
            self.dlg.cbAlignmentType.setCurrentIndex(1)
            
            self.dlg.cbRemoveTooDistant.clear()
            self.dlg.cbRemoveTooDistant.addItems(["Yes","No"])
            self.dlg.cbRemoveTooDistant.setCurrentIndex(0)
            
            self.dlg.cbRemoveTooClose.clear()
            self.dlg.cbRemoveTooClose.addItems(["Yes","No"])
            self.dlg.cbRemoveTooClose.setCurrentIndex(1)
            
            self.dlg.cbDistanceToRemove.clear()
            self.dlg.cbDistanceToRemove.addItems(list(map(str,range(1,11))))
            self.dlg.cbDistanceToRemove.setCurrentIndex(2)
        
            self.dlg.cbRippled.clear()
            self.dlg.cbRippled.addItems(["Yes","No"])
            self.dlg.cbRippled.setCurrentIndex(0)
        
            self.alignment_loaded = True 
        
        
        
        
        self.dlg.teAlignedReactions.clear()
        
        self.dlg.cbAlignmentReactions.clear()
        reaction_names = list(self.model.reactions)
        reaction_names.sort()
        self.dlg.cbAlignmentReactions.addItems(reaction_names) 
        
        
     
    def select_source_type(self):
        source_type = self.dlg.cbSourceType.currentText()
        if source_type == 'KEGG pathway':
            self.dlg.gbSourcePathway.show()
            self.dlg.gbSourceModel.hide()
        else: #'SBML'
            self.dlg.gbSourcePathway.hide()
            self.dlg.gbSourceModel.show()
    
    
    def openModelSrcFileNameDialog(self):   

        #options = QFileDialog.Options()
        #fileName = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        
        fname = QFileDialog.getOpenFileName(self.dlg, 'Open model', '', 'Model files (*.mat *.xml);; Matlab files (*.mat);; SBML files (*.xml);; All Files (*)')[0]
        
        if fname:
            self.model_src = md.model_data()
            _ , ext = splitext(fname)
            #self.dlg.setWindowTitle("Grohar - " + model_name + ext)
            if ext == ".mat":
                self.model_src.load_model_mat(fname)
            elif ext == ".xml":
                self.model_src.load_model_cobra_sbml(fname)                
                
            self.dlg.cbAlignmentCompartmentsSrc.clear()
            self.dlg.cbAlignmentCompartmentsSrc.addItem("None")
            self.dlg.cbAlignmentCompartmentsSrc.addItems(self.model_src.compartments)
            index = self.dlg.cbAlignmentCompartmentsSrc.findText('c', QtCore.Qt.MatchFixedString)
            if index >= 0:
                self.dlg.cbAlignmentCompartmentsSrc.setCurrentIndex(index)
                
        
    
    def load_bigg_data(self):
       
        bigg_id = self.dlg.cbBiGGModels.currentText() 
        if bigg_id != 'None':
            self.model_meta_data = model_meta_data(bigg_id, from_bigg = 1)       
            self.model_meta_data.load_metabolite_kegg_ids()
            self.model_meta_data.load_reaction_ECs()
            print("BiGG data loaded")
        else:
            self.model_meta_data = None
            
            ret = QMessageBox.question(self.dlg, '', "BiGG id is not selected. Try to load the data from KEGG?", QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if ret == QMessageBox.Yes:
                self.load_model_kegg_data()
            else:
                print("BiGG data not loaded")
    
    def load_model_kegg_data(self):                
        model_id = self.model.model.id                
        self.model_meta_data = model_meta_data(model_id, model = self.model, from_bigg = 0)
        self.model_meta_data.load_reaction_ECs()

    
    
    def align(self):
        
        
        #reload = 0        
        
        #pathway_id = "cge00250"
        #pathway = pd.pathway_data(pathway_id)

       

        #if reload:
        #    pathway.load_kegg("cge00250")
        #    pathway.save_ECs("cge00250_ec.txt")
        #    pathway.save_kegg("cge00250.xml")
        #else:
        #    pathway.load_kgml("cge00250.xml", "cge00250_ec.txt")
       
        
        
    
    
    
    
        # load parameters from the GUI
        compartment = self.dlg.cbAlignmentCompartments.currentText()
        if compartment == "None":
            compartment = ""
        
        max_metabolite_occurence = int(self.dlg.leAlignmentMaxMetOccurance.text())
            
        
        if self.dlg.cbReversibility.currentText() == "Yes":
            reversibility = 1
        else:
            reversibility = 0
        
        K_m = int(self.dlg.cbKm.currentText())
        
        
        
        if self.dlg.cb1ToMany.currentText() == "Yes":
            to_many = 1
        else:
            to_many = 0
        
        if self.dlg.cbManyTo1.currentText() == "Yes":
            from_many = 1
        else:
            from_many = 0
        
        if self.dlg.cbRemoveBorders.currentText() == "Yes":
            remove_borders = 1
        else:
            remove_borders = 0
        
        topological_weight = 10 * int(self.dlg.cbTopologicalWeight.currentText())
        homological_weight_EC = 10 * int(self.dlg.cbECWeight.currentText())
        homological_weight_KEGG = 10 * int(self.dlg.cbKEGGWeight.currentText())
        
        if self.dlg.cbAdditive.currentText() == "Additive":
            homological_additive = 1
        else:
            homological_additive = 0
        
        alignment_type = self.dlg.cbAlignmentType.currentText()
        # ["Simple", "Greedy", "Probabilistic"]
        
        if self.dlg.cbRemoveTooDistant.currentText() == "Yes":
            remove_too_distant = 1
        else:
            remove_too_distant = 0
            
        if self.dlg.cbRemoveTooClose.currentText() == "Yes":
            remove_too_close = 1
        else:
            remove_too_close = 0
        
        if self.dlg.cbRippled.currentText() == "Yes":
            max_from_neighbours = 1
        else:
            max_from_neighbours = 0
        
        
        distance_to_remove = int(self.dlg.cbDistanceToRemove.currentText())
        
        
        align = alignment()
        print("alignment object created")
        
        source_type = self.dlg.cbSourceType.currentText()
        
        can_proceed = False
        if source_type == 'KEGG pathway':
            #pathway = pd.pathway_data("cge00100")
            pathway = pd.pathway_data(self.dlg.lePathwayId.text())
            pathway.load_kgml()
            align.load_model_pathway(pathway, self.model, model_meta_data = self.model_meta_data, max_metabolite_occurence = max_metabolite_occurence, compartment = compartment, reversibility = reversibility )
            print("model and pathway loaded")
            can_proceed = True
        else: # 'sbml'
            if self.model_src:
               
                compartment_src = self.dlg.cbAlignmentCompartmentsSrc.currentText()
                if compartment_src == "None":
                    compartment_src = ""
                
                align.load_model_model(self.model, self.model_src, max_metabolite_occurence = max_metabolite_occurence, compartment_dst = compartment, compartment_src = compartment_src, reversibility = reversibility)
                print("both models loaded")
                can_proceed = True
        if can_proceed:
            
        
            a, score = align.make_alignment(alignment_type, K_m = K_m, to_many = to_many, from_many = from_many, remove_borders = remove_borders, topological_weight = topological_weight, homological_weight_EC = homological_weight_EC, homological_weight_KEGG = homological_weight_KEGG, homological_additive = homological_additive, remove_too_distant = remove_too_distant, remove_too_close = remove_too_close, distance_to_remove = distance_to_remove, max_from_neighbours = max_from_neighbours)
            
            #a, score = align.make_alignment(3, to_many = 0, from_many = 0, topological = 0, homological_weight_EC = 10, homological_weight_KEGG = 10, homological_additive = 0)
            align.print_alignment(a)
            #align.print_alignment_full(a)
            reacts = []
            
            """
            f = open("debug_KEGG.txt","w")
            for i in align.KEGG_metabolites_dst:
                for j in i[0]:
                    f.write(j + " ")
                f.write("\n")
                for j in i[1]:
                    f.write(j + " ")
                f.write("\n")
                f.write("\n")
            f.close()
                
            f = open("debug_EC.txt","w")
            for i in align.EC_dst:
                for j in i:
                    f.write(j + " ")
                f.write("\n")            
            f.close()
            """ 
            
            
            
            for i in a:
                reacts.append(align.reactions_dst[a[i]])
                
            """
            compartments = set(self.dlg.teCompartments.toPlainText().strip().split("\n")) - {""}
            compartments.add(self.dlg.cbCompartments.currentText())
            compartments = list(compartments)
            compartments.sort()
            compartments_disp = "\n".join(compartments)
            self.dlg.teCompartments.setText(compartments_disp)
            """
            
            reacts0 = set(self.dlg.teAlignedReactions.toPlainText().strip().split("\n")) - {""}
            reacts = reacts0 | set(reacts)      
            
            self.dlg.teAlignedReactions.setText("\n".join(reacts))

            self.alignment = align
            self.pathway_data = pathway
            
            
            #self.draw_reactions(reacts)          
            #self.draw_reactions(["ASNN","r1580","r1637"])
      
    def add_alignment_reaction(self):
        reaction = self.dlg.cbAlignmentReactions.currentText()
        reacts = set(self.dlg.teAlignedReactions.toPlainText().strip().split("\n")) - {""}
        reacts.add(reaction)
        reacts = list(reacts)
        self.dlg.teAlignedReactions.setText("\n".join(reacts))
        
    
    def save_aligned_reactions(self):
        fname = QFileDialog.getSaveFileName(self.dlg, 'Save file', '', 'All Files (*)')[0]
        
        f = open(fname, 'w')
        reacts = set(self.dlg.teAlignedReactions.toPlainText().strip().split("\n")) - {""}
        for i in reacts:
            f.write(i+"\n")
        f.close()
    
    
    def load_aligned_reactions(self):
        fname = QFileDialog.getOpenFileName(self.dlg, 'Open file', '', 'All Files (*)')[0]
        
        if fname:
            try:
                f = open(fname, 'r')
                reacts = set()
                for i in f:
                    reacts.add(i.strip())
                reacts = reacts - {""}
                self.dlg.teAlignedReactions.setText("\n".join(reacts))
                f.close()               
            except:
                print("File read error!")
    
    
    def visualise_alignment(self):
        reacts = set(self.dlg.teAlignedReactions.toPlainText().strip().split("\n")) - {""}
                
        self.draw_reactions(reacts, max_metabolite_occurance = int(self.dlg.leAlignmentMaxMetOccurance.text()))          
    
    def export_aligned_reactions(self):
         reacts = set(self.dlg.teAlignedReactions.toPlainText().strip().split("\n")) - {""}
         if reacts:
             self.export_reactions(reacts)
    

    
    def connect_KEGG(self):
        if not self.kegg:
            self.kegg = KEGG()
            
        organisms =  self.kegg.organismIds
        self.dlg.cbOrganism.clear()
        self.dlg.cbOrganism.addItems(organisms)
        
    def select_organism_KEGG(self):
        self.kegg.organism = self.dlg.cbOrganism.currentText()
        pathways = self.kegg.pathwayIds
        pathways = [p.replace("path:","") for p in pathways]
        self.dlg.cbPathway.clear()
        self.dlg.cbPathway.addItems(pathways)
            
        
    def select_pathway(self):
        self.dlg.lePathwayId.setText(self.dlg.cbPathway.currentText())
    #def select_pathway_KEGG(self):
    #    pass
   
    #########
    #########
    # PATHS #
    #########
    #########
        
    def shortest_path(self):
        node1 = self.dlg.lePathSourceNode.text()
        node2 = self.dlg.lePathDstNode.text()

        ignore_weights = self.dlg.cbIgnoreWeights.checkState()
        

        if (node1 in self.model.reactions or node1 in self.model.metabolites)  and (node2 in self.model.reactions or node2 in self.model.metabolites):  
            
            max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text())
            
            
            P = paths(self.model, ignore_weights, max_metabolite_occurance = max_metabolite_occurance)
            try:
                all_paths = P.get_all_shortest_paths(node1, node2)
                all_paths = list(all_paths)
               
                #self.draw_reactions([i for i in all_paths[0] if i in self.model.reactions], observed = [node1, node2]) #only reactions
                self.draw_reactions([i for i in all_paths[0] if i in self.model.reactions], max_metabolite_occurance = max_metabolite_occurance) #only reactions
                                   
                all_paths = [" --> ".join([i for i in p if i in self.model.reactions]) for p in all_paths]
                self.dlg.cbAllShortestPaths.clear()
                self.dlg.cbAllShortestPaths.addItems(all_paths)           
            except:
                self.dlg.cbAllShortestPaths.clear()
                self.dlg.cbAllShortestPaths.addItems(["Path does not exist! Increase maximal metabolite occurrence!?"])           
                
        else:
            self.dlg.cbAllShortestPaths.clear()
            self.dlg.cbAllShortestPaths.addItems(["Wrong node name(s)!"])
            
            
    def shortest_paths_perturbation(self):
        plot_all = not self.dlg.cbPathPerturbationDifferences.checkState()
                
        node1 = self.dlg.lePathSourceNode.text()
        node2 = self.dlg.lePathDstNode.text()
        
        if (node1 in self.model.reactions or node1 in self.model.metabolites)  and (node2 in self.model.reactions or node2 in self.model.metabolites):  
            
            max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text())
            
            
            (reactions, lower, upper) = self.get_constraints()  
            boundaries = []
            for (r, l, u) in zip(reactions, lower, upper):
                boundaries.append([r, 'l', l]) 
                boundaries.append([r, 'u', u]) 
                     
            objectives = self.parse_objective()
                    
            self.model.run_FBA()              
            self.model.run_perturbation(boundaries, objectives = objectives)      
            
            if not self.model.sol or not self.model.perturbed_sol:      
                QMessageBox.warning(self.dlg, "Warning","Solution not found!")
                return
            
            
            
            
            try:
                self.dlg.cbAllShortestPaths.clear()
                
                                            
                P = paths(self.model, ignore_weights = 0, max_metabolite_occurance = max_metabolite_occurance)
                all_paths = P.get_all_shortest_paths(node1, node2)
                all_paths = list(all_paths)
                #print(all_paths)
                reactions1 = [r for path in all_paths for r in path if r in self.model.reactions]
                #print(reactions1)
                
                                
                P_perturbed = paths(self.model, ignore_weights = 0, max_metabolite_occurance = max_metabolite_occurance, perturbed = 1)
                all_paths_perturbation = P_perturbed.get_all_shortest_paths(node1, node2)
                all_paths_perturbation = list(all_paths_perturbation)
                #print(all_paths_perturbation)
                reactions2 = [r for path in all_paths_perturbation for r in path if r in self.model.reactions]
                #print(reactions2)
                
                                
                if plot_all:
                    comparison_type = 4
                    (pairs, fluxes) = make_pairs.compare_perturbation_reactions(comparison_type, self.model, reactions1, reactions2, max_metabolite_occurance = max_metabolite_occurance)                                        
                    #print(reactions1)
                    #print(reactions2)
                    self.draw_pathway_perturbations(pairs, fluxes, reactions1, reactions2)                     
                else:
                    comparison_type = 3
                    (pairs, fluxes) = make_pairs.compare_perturbation_reactions(comparison_type, self.model, reactions1, reactions2, max_metabolite_occurance = max_metabolite_occurance)
                    met_names = []#[node1, node2]
                    reacts = list(set(reactions1 + reactions2))
                    self.draw(met_names, pairs, fluxes, path = reacts)
                
               
            except Exception as inst:
                self.dlg.cbAllShortestPaths.clear()
                self.dlg.cbAllShortestPaths.addItems(["Path does not exist!"])           
                print(inst)
                
        else:
            self.dlg.cbAllShortestPaths.clear()
            self.dlg.cbAllShortestPaths.addItems(["Wrong node name(s)!"])
   

    def draw_pathway_perturbations(self, pairs, fluxes, reactions1, reactions2):
        model = self.model
        
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
    
        bbox_props = dict(boxstyle="round,pad=0.3", ec="0.5", lw=2)
        self.label = self.ax.text(100, 100, "text", color = "white", bbox=bbox_props, zorder=10)
        self.label.set_visible(False)
        
        
     
        pairs_mod = []
        for i,(a,b) in enumerate(pairs):
            pairs_mod.append((a.replace(":","---"), b.replace(":","---"))) # networkx uses ":" as a special symbol!
      
            
        #G = nx.DiGraph(splines = True)
        G = nx.DiGraph()
        
        #nodes = list(set(i for pair in pairs_mod for i in pair))
        #G.add_nodes_from(nodes)
        
        
        
        G.add_edges_from(pairs_mod)
        
        if self.dlg.actionCBFriendly.isChecked():
            
 
            val_map = {'metabolite': '#f0e442',
                   'pathway1': '#009e73',
                   'pathway2': '#d55e00',
                   'both': '#7570b3'}
        
            patches = [mpatches.Patch(color='#f0e442', label='metabolites'),
                   mpatches.Patch(color='#009e73', label='before'),
                   mpatches.Patch(color='#d55e00', label='after'),                                 
                   mpatches.Patch(color='#7570b3', label='both')]
                     
        else:
            
           val_map = {'metabolite': 'c',
                   'pathway1': 'g',
                   'pathway2': 'r',
                   'both': 'y'}
        
           patches = [mpatches.Patch(color='c', label='metabolites'),
                   mpatches.Patch(color='g', label='before'),
                   mpatches.Patch(color='r', label='after'),                                 
                   mpatches.Patch(color='y', label='both')]
        
        
        
        
        
        
        
        
        self.lgd = plt.legend(handles = patches)
        self.lgd_visible = False
        self.lgd.set_visible(self.lgd_visible)
        
        
        #node_colors = [val_map['observed'] if node in met_names else val_map['both'] if node in producing & consuming else val_map['producing'] if node in producing else val_map['consuming'] if node in consuming else val_map['reaction'] if node in model.reactions else val_map['metabolite'] for node in G.nodes()]
        node_colors = []
        for node in G.nodes():
          
            node_orig = node.replace("---",":")
            
            if node_orig in set(reactions1) & set(reactions2):
                color = val_map['both']
            elif node_orig in reactions1:
                color = val_map['pathway1']
            elif node_orig in reactions2:
                color = val_map['pathway2']
            else:
                color = val_map['metabolite']             
            node_colors.append(color)            
        

        self.node_ids = {node:i for i, node in enumerate(G.nodes())}

        #node_labels = {a:a + ":" + "{:.5f}".format(fluxes[a]) if a in model.reactions else a for a in G.nodes()}
        node_labels = {}
        for node in G.nodes():                             
            node_orig = node.replace("---",":")
            if node_orig in model.reactions:
                f = fluxes[node_orig]
                if type(f) == tuple:
                    node_labels[node] = node_orig + " :" + "{:.3f}\n{:.3f}".format(round(f[0],3),round(f[1],3))
                else:
                    node_labels[node] = node_orig + " :" + "{:.3f}".format(round(fluxes[node_orig],3))
            else:
                node_labels[node] = node_orig
        
        #node_sizes = [500 if a in model.reactions else 250 for a in G.nodes()]
        node_sizes = []
        for node in G.nodes():
            node_orig = node.replace("---",":")
            if node_orig in model.reactions:
                node_sizes.append(500)
            else:
                node_sizes.append(250)
        
        
        
       
        pos = nx.nx_pydot.graphviz_layout(G, prog="dot")
        
        
    
    
        nx.draw_networkx_labels(G, pos , node_labels, font_size=8, font_weight = "bold", font_color = "#1c2833", ax = self.ax)
        #nx.draw_networkx_labels(G,labels=node_labels)
        
       
        nx.draw(G, pos, node_color = node_colors, node_size=node_sizes, edge_color = "#abb2b9", dge_cmap=matplotlib.cm.Reds, alpha = 0.5, ax = self.ax)
      
        
        ##nx.draw_networkx_nodes(G, pos, nodelist = [a for a in G if a in model.reactions], node_color = node_colors, node_size=250,  alpha=0.5)
        ##nx.draw_networkx_nodes(G, pos, nodelist = [a for a in G if a in model.metabolites], node_color = node_colors, node_size=100,  alpha=0.5)
        ##nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
        
        ##pylab.show()
        
        
        #self.dlg.canvas.draw_idle()
        
        #handles, labels = self.ax.get_legend_handles_labels()
        
        
        self.G = G
        self.pos = pos
        #self.pos = nx.get_node_attributes(G,'pos')
        self.node_labels = node_labels
        self.node_colors = node_colors
        self.node_sizes = node_sizes
        
      
        
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_plot_hover)  
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.fig.canvas.mpl_connect('button_press_event', self.on_plot_click) 
        self.fig.canvas.mpl_connect('button_release_event', self.on_plot_click_release) 

     
    def visualise_selected_path(self):
        #node1 = self.dlg.lePathSourceNode.text()
        #node2 = self.dlg.lePathDstNode.text()
        
        path = self.dlg.cbAllShortestPaths.currentText().split(" --> ")
        #self.draw_reactions(path, observed = [node1, node2])
        self.draw_reactions(path, max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text()))


    def visualise_all_paths(self):
        #node1 = self.dlg.lePathSourceNode.text()
        #node2 = self.dlg.lePathDstNode.text()
        
        all_paths = [self.dlg.cbAllShortestPaths.itemText(i) for i in range(self.dlg.cbAllShortestPaths.count())]
        #print(all_paths)
        all_paths = list(set(r for path in all_paths for r in path.split(" --> ")))
        #print(all_paths)
        #self.draw_reactions(all_paths, observed = [node1, node2])
        self.draw_reactions(all_paths, max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text()))
        
    def draw_steiner(self):
        nodes = list(set(self.dlg.teSteiner.toPlainText().strip().split("\n")) - {""})
        
        if nodes:       
            ignore_weights = 1
            max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text())
            P = paths(self.model, ignore_weights = ignore_weights, max_metabolite_occurance = max_metabolite_occurance)         
            tree = P.get_sub_graph(nodes)
            #print(tree)
            reactions = [t for t in tree if t in self.model.reactions]
            #print(reactions)
            self.draw_reactions(reactions, observed = tree, max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text()))
            
            
    def visualise_path_multi(self):
        nodes = self.dlg.teSteiner.toPlainText().strip().split("\n")
        #print("Nodes: ", nodes)
        
        if nodes:       
            try:
                ignore_weights = self.dlg.cbIgnoreWeights.checkState()
                max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text())
                P = paths(self.model, ignore_weights = ignore_weights, max_metabolite_occurance = max_metabolite_occurance)
                                
                path = []
                for i in range(len(nodes) - 1):
                    node1 = nodes[i]
                    node2 = nodes[i+1]
                    path += P.get_shortest_path(node1, node2)
                
                path = list(set(path))
                reactions = [r for r in path if r in self.model.reactions]
            
                self.draw_reactions(reactions, observed = path, max_metabolite_occurance = int(self.dlg.lePathMaxMetOccurrance.text()))
            except Exception as inst:
                print(inst)
        
    ##############
    ##############
    # EDIT MODEL #
    ##############
    ##############
    
    def init_edit_reactions(self):
   
        self.dlg.leAddReactionID.setText("")
        self.dlg.leAddReactionName.setText("")
        self.dlg.cbEditReactions.clear()
        self.dlg.cbEditReactionsMetabolites.clear()
        self.dlg.cbEditReactionsReactants.clear()
        self.dlg.cbEditReactionsReactants.clear()
        self.dlg.leEditReactionsLB.setText("")
        self.dlg.leEditReactionsUB.setText("")
  
        
        
    def add_reaction(self):
        r_id = self.dlg.leAddReactionID.text()
        r_name = self.dlg.leAddReactionName.text()
        
        if self.model.add_reaction_basic(r_id, r_name):
            self.update_reaction_cbs()
        
      
    def select_edit_reaction(self):
            
        r_id = self.dlg.cbEditReactions.currentText()
        #print(r_id)

        r_name = self.model.reactions[r_id]        
        #print(r_name)
        self.dlg.leEditReactionsName.setText(r_name)
        
        reactants = self.model.reaction_reactants[r_id]
        #print(reactants)
        self.dlg.cbEditReactionsReactants.clear()
        self.dlg.cbEditReactionsReactants.addItems(reactants)


        products = self.model.reaction_products[r_id]
        #print(products)
        self.dlg.cbEditReactionsProducts.clear()
        self.dlg.cbEditReactionsProducts.addItems(products)


        lb, ub = self.model.reaction_bounds[r_id]
        #print(lb, ub)
        self.dlg.leEditReactionsLB.setText(str(lb))
        self.dlg.leEditReactionsUB.setText(str(ub))
            
            
    def edit_reactions_commit(self):
        try:
            r_id = self.dlg.cbEditReactions.currentText()
            r_name = self.dlg.leEditReactionsName.text()
            reactants = [self.dlg.cbEditReactionsReactants.itemText(i) for i in range(self.dlg.cbEditReactionsReactants.count())]
            products = [self.dlg.cbEditReactionsProducts.itemText(i) for i in range(self.dlg.cbEditReactionsProducts.count())]
            lb = float(self.dlg.leEditReactionsLB.text())
            ub = float(self.dlg.leEditReactionsUB.text())
            self.model.edit_reaction(r_id, r_name, reactants, products, [lb, ub])
        except Exception as inst:
            print(inst)                

        
    def delete_reaction(self):
        r_id = self.dlg.cbEditReactions.currentText()
        if self.model.delete_reaction(r_id):            
            self.update_reaction_cbs()
        
    def add_reactant(self):
        m_id = self.dlg.cbEditReactionsMetabolites.currentText()
        reactants = [self.dlg.cbEditReactionsReactants.itemText(i) for i in range(self.dlg.cbEditReactionsReactants.count())]
        reactants.append(m_id)
        reactants = list(set(reactants))
        reactants.sort()
        self.dlg.cbEditReactionsReactants.clear()
        self.dlg.cbEditReactionsReactants.addItems(reactants)
        
        
    def add_product(self):
        m_id = self.dlg.cbEditReactionsMetabolites.currentText()
        products = [self.dlg.cbEditReactionsProducts.itemText(i) for i in range(self.dlg.cbEditReactionsProducts.count())]
        products.append(m_id)
        products = list(set(products))
        products.sort()
        self.dlg.cbEditReactionsProducts.clear()
        self.dlg.cbEditReactionsProducts.addItems(products)
        
        
    def remove_reactant(self):
        m_id = self.dlg.cbEditReactionsReactants.currentText()
        reactants = [self.dlg.cbEditReactionsReactants.itemText(i) for i in range(self.dlg.cbEditReactionsReactants.count()) if self.dlg.cbEditReactionsReactants.itemText(i) != m_id]
        reactants = list(set(reactants))
        reactants.sort()
        self.dlg.cbEditReactionsReactants.clear()
        self.dlg.cbEditReactionsReactants.addItems(reactants)
        
    def remove_product(self):
        m_id = self.dlg.cbEditReactionsProducts.currentText()
        products = [self.dlg.cbEditReactionsProducts.itemText(i) for i in range(self.dlg.cbEditReactionsProducts.count()) if self.dlg.cbEditReactionsProducts.itemText(i) != m_id]
        products = list(set(products))
        products.sort()
        self.dlg.cbEditReactionsProducts.clear()
        self.dlg.cbEditReactionsProducts.addItems(products)
        
    def add_metabolite(self):
        m_id = self.dlg.leAddMetaboliteID.text()
        m_name = self.dlg.leAddMetaboliteName.text()
        compartment_id = self.dlg.cbEditMetabolitesCompartments.currentText()
        
        if self.model.add_metabolite(m_id, m_name, compartment_id):
            self.update_metabolite_cbs()
   
    def delete_metabolite(self):
        m_id = self.dlg.cbEditMetabolitesMetabolites.currentText()
        
        if self.model.delete_metabolite(m_id):
            self.update_metabolite_cbs()
    
    def select_edit_metabolite(self):
        m_id = self.dlg.cbEditMetabolitesMetabolites.currentText()
        self.dlg.leEditMetaboliteName.setText(self.model.metabolites[m_id])
     
    
    def edit_metabolite(self):
        m_id =  self.dlg.cbEditMetabolitesMetabolites.currentText()
        m_name = self.dlg.leEditMetaboliteName.text()
        
        self.model.edit_metabolite(m_id, m_name)
    
    
    def add_compartment(self):
        c_id = self.dlg.leAddCompartmentID.text()
        c_name = self.dlg.leAddCompartmentName.text() 
        
        if self.model.add_compartment(c_id, c_name):
            self.update_compartments_cbs()
        
    def delete_compartment(self):
        c_id = self.dlg.cbDeleteCompartments.currentText()
        
        if self.model.delete_compartemnt(c_id):
            self.update_compartments_cbs()
        
    def add_objective_2(self):
        try:   
            objectives = self.parse_objective(option = 2)               
            r_id  = self.dlg.cbEditObjective.currentText()
            objectives[r_id] += 1
                
            output = ""
            for i,r_id in enumerate(objectives):
                if i > 0:
                    output += " + "
                output += r_id + " * " + str(objectives[r_id])
            self.dlg.leEditObjective.setText(output)                                 
        except Exception as inst:
            print(inst)
    
    def objective_commit(self):
        try:   
            objectives = self.parse_objective(option = 2)               
            objective = {self.model.model.reactions[self.model.reaction_position[r_id]]:int(objectives[r_id]) for r_id in objectives}
            self.model.update_objective(objective)
            
        except Exception as inst:
            print(inst)
    
        
app = QApplication([])        
t = top()     
app.exec()
