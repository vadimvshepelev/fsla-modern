import numpy as np

class CRiemannSolver:
    """Класс реализует различные римановские солверы"""
    def __init__(self, problem):
        self.GAMMA = problem.gamma
		
def calc_flux(const, U_l, U_r):
    """Calculates intercell flux based on exact Riemann problem solution"""
    F = np.zeros(5)
    return F


def solve_RP(ro_l, u_l, p_l, ro_r, u_r, p_r):	
    """Function finds resulting ro, v, p of Riemann task solving nonlinear equation by tangentials method of Newton"""	
    res = CRPSolutionPrimitive(0., 0., 0., 0., "")
    tol = 1.e-6
    c_l = 0
    c_r = 0
    if ro_l != 0:
        c_l = sqrt(GAMMA*p_l/ro_l)
    if ro_r != 0.:
        c_r = sqrt(GAMMA*p_r/ro_r)
    # Trivial case UL = UR
    if ro_l == ro_r and u_l == u_r and p_l == p_r:
        res.ro_l = ro_l
        res.ro_r = ro_r
        res.p = p_l
        res.u = u_l
        return res
    # Vacuum at the left side
    if ro_l == 0.:
        res.s_type = "VacRW"
        res.ro_l = 0.
        res.ro_r = ro_r
        res.p = 0.
        res.u = u_r - 2.*c_r/(GAMMA-1.)
        return res
    # Vacuum at the right side
    if ro_r == 0.:
        res.s_type = "RWVac"
        res.ro_l  = ro_l
        res.ro_r  = 0.
        res.p = 0.
        res.u = u_l + 2.*c_l/(GAMMA-1.)
        return res
    # Two rarefaction waves generating vacuum
    if 2.*c_l/(GAMMA-1.) + 2.*c_r/(GAMMA-1.) < u_r-u_l:
        res.type = "RWVacRW"
        res.ro_l = 0.
        res.ro_r = 0.
        res.u = 0.
        res.p = 0.
        return res
    # General case -- tangentials
    # fl_min = fl(p_l, ro_l, u_l, p_l) + fr(p_l, ro_r, u_r, p_r) + u_r-u_l
    # fr_max = fl(p_r, ro_l, u_l, p_l) + fr(p_r, ro_r, u_r, p_r) + u_r-u_l
    # Initial approximation
    p = p_l/2.
    it_c = 0
    while True:
        p_prev = p
        p = p_prev - \
        (fl(p_prev, ro_l, u_l, p_l) + fr(p_prev, ro_r, u_r, p_r) + u_r - u_l) / \
        (dfldp(p_prev, ro_l, u_l, p_l) + dfrdp(p_prev, ro_r, u_r, p_r))
        if p <= 0.:
            p = tol
        it_c += 1
        if fabs(2*(p-p_prev)/(p+p_prev))<=tol:
            break
    res.p = p
    res.u = 0.5*(u_l + u_r) + 0.5*(fr(p, ro_r, u_r, p_r) - fl(p, ro_l, u_l, p_l))
    if p < p_l and p > p_r:
        res.s_type = "RWSW"
        res.ro_l  = ro_l*pow(res.p/p_l, 1./GAMMA)
        res.ro_r  = ro_r*(res.p/p_r + (GAMMA-1.)/(GAMMA+1.))/((GAMMA-1.)/(GAMMA+1.)*res.p/p_r + 1.)
    elif p <= p_l and p <= p_r:
        res.s_type = "RWRW"
        res.ro_l = ro_l*pow(res.p/p_l, 1./GAMMA)
        res.ro_r = ro_r*pow(res.p/p_r, 1./GAMMA)
    elif p > p_l and p < p_r:
        res.s_type = "SWRW"
        res.ro_l = ro_l*(res.p/p_l + (GAMMA-1.)/(GAMMA+1.))/((GAMMA-1.)/(GAMMA+1.)*res.p/p_l + 1.)
        res.ro_r = ro_r*pow(res.p/p_r, 1./GAMMA)
    else:
        res.s_type = "SWSW"
        res.ro_l = ro_l*(res.p/p_l + (GAMMA-1.)/(GAMMA+1.))/((GAMMA-1.)/(GAMMA+1.)*res.p/p_l + 1.)
        res.ro_r = ro_r*(res.p/p_r + (GAMMA-1.)/(GAMMA+1.))/((GAMMA-1.)/(GAMMA+1.)*res.p/p_r + 1.)
#	//Testing (no vacuum)
#	/*double L =  - fL(p, roL, vL, pL);
#	double R = fR(p, roR, vR, pR) + vR - vL;
#	double delta = fabs((L-R)/0.5/(L+R));
#	double LToro = - fL(res.p, roL, vL, pL);
#	double RToro = fR(res.p, roR, vR, pR) + vR - vL;
#	double deltaToro = fabs((LToro-RToro)/0.5/(LToro+RToro));*/
    return res

def fl(p, ro_l, u_l, p_l):
    if p>p_l:
        al = 2./(GAMMA+1)/ro_l
        bl = (GAMMA-1.)/(GAMMA+1.)*p_l
        f = (p-p_l) * sqrt(al/(p+bl))
        return f
    else:
        c_l = sqrt(GAMMA*p_l/ro_l)
        f = 2.*c_l/(GAMMA-1.) * ( (pow(p/p_l, (GAMMA-1.)/2./GAMMA)) - 1. )
        return f

def dfldp(p, ro_l, u_l, p_l):
    if p>p_l:
        al = 2./(GAMMA+1)/ro_l
        bl = (GAMMA-1.)/(GAMMA+1.)*p_l
        dfdp = sqrt(al/(p+bl)) * (1. - (p-p_l)/2./(p+bl))
        return dfdp
    else:
        c_l = sqrt(GAMMA*p_l/ro_l)
        dfdp = c_l/p_l/GAMMA*pow(p/p_l, -(GAMMA+1)/2./GAMMA)
        return dfdp

def fr(p, ro_r, u_r, p_r):
    if p > p_r:
        ar = 2./(GAMMA+1)/ro_r
        br = (GAMMA-1.)/(GAMMA+1.)*p_r
        f = (p-p_r) * sqrt(ar/(p+br))
        return f
    else:
        c_r = sqrt(GAMMA*p_r/ro_r)
        f = 2.*c_r/(GAMMA-1.) * ((pow(p/p_r, (GAMMA-1.)/2./GAMMA)) - 1.)
        return f

def dfrdp(p, ro_r, u_r, p_r):
    if p > p_r:
        AR = 2./(GAMMA + 1) / ro_r
        BR = (GAMMA - 1.) / (GAMMA + 1.) * p_r
        dfdp = sqrt(AR/(p+BR)) * (1. - (p-p_r)/2./(p+BR))
        return dfdp
    else:
        cR = sqrt(GAMMA * p_r / ro_r)
        dfdp = cR/p_r / GAMMA * pow(p / p_r, -(GAMMA + 1) / 2. / GAMMA)
        return dfdp