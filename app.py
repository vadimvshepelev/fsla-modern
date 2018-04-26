class CApp:
    """Implements the whole calculation process"""
    def __init__(self, problem, eos, field, solver):
        print("Class CApp: Initializing numerical experimnet process...", end="")        
        self.problem = problem
        self.eos = eos
        self.field = field
        self.solver = solver	
        self.t = self.problem.t_min
        print("done!")

    def run(self):
        print("Starting computational process...")	
        self.field.set(self.problem, self.eos)