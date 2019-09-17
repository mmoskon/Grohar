class FBA_data:
    
    def __init__(self, model):
        self.model = model
        self.run_FBA()
        
    def run_FBA(self):
        model = self.model
        #if model.type == 'cobra':
        self.sol = model.model.optimize()
        self.fluxes = self.sol.x_dict
        
        