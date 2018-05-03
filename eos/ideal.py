from math import sqrt

class EOSIdeal:
    """Ideal equation of state"""
    def __init__(self, GAMMA):
        print("Class eos.CEOSIdeal: Initializing ideal equation of state...", end="")
        self.GAMMA = GAMMA  
        print("done!")
        
    def getp(self, ro, e):
        return (self.GAMMA-1.)*ro*e
        
    def gete(self, ro, p):
        return p/(self.GAMMA-1.)/ro
        
    def getc(self, ro, e):
        p = self.getp(ro, e)
        return sqrt(self.GAMMA*p/ro)
        