# Import package and module
import numpy as np
import utility as util
import tools


def EGM(sol,t,par):
    sol = EGM_loop(sol,t,par) 
    #sol = EGM_vec(sol,t,par) 
    return sol

def EGM_loop (sol,t,par):
    for i_a,a in enumerate(par.grid_a[t,:]):
        if t+1 <= par.Tr: # no retirement in next period
            fac = par.G*par.L[t]*par.psi_vec
            w = par.w
            xi = par.xi_vec
            inv_fac = 1/fac

            # Future m and c
            m_next = inv_fac*par.R*a+xi
            c_next = tools.interp_linear_1d(sol.m[t+1,:],sol.c[t+1,:], m_next)
        else: 
            fac = par.G*par.L[t]
            w = 1
            xi = 1
            inv_fac = 1/fac

            # Future m and c
            m_next = inv_fac*par.R*a+xi
            c_next = tools.interp_linear_1d_scalar(sol.m[t+1,:],sol.c[t+1,:], m_next)

        # Future marginal utility
        marg_u_next = util.marg_util(fac*c_next,par)
        avg_marg_u_next = np.sum(w*marg_u_next)

        # Correct C and m
        sol.c[t,i_a+1] = util.inv_marg_util(par.beta*par.R*avg_marg_u_next,par)
        sol.m[t,i_a+1] = a+sol.c[t,i_a+1]

    return sol

def EGM_vec (sol,t,par):

    

    return sol
