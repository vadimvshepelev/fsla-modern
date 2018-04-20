class EOSIdeal:
    """Ideal equation of state"""
    def __init__(self, GAMMA):
        print("Class eos.CEOSIdeal: Initializing ideal equation of state...", end="")
        self.GAMMA = GAMMA  
        print("done!")