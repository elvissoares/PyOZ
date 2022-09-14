import numpy as np
from scipy.fftpack import dst, idst

def uhs(r,sigma):
    return np.piecewise(r,[r<=sigma,r>sigma],[np.inf,0.0])

def closure_relation(beta,u,gamma,model='PY'):
    if model == 'PY':
        b_r = np.log(1+gamma)-gamma
    elif model == 'HNC':
        b_r = 0.0
    return np.exp(-beta*u+gamma+b_r)-gamma-1

def PyOZ(r,rho,kBT,u=uhs,params=1.0,model='PY'):
    alpha = 0.001
    finc = 1.1
    fdec = 0.5
    atol = 1.e-5
    max_iter = 9999
    n_iter = 0
    # Fourier variables
    n_pts = r.size
    dr = r[1]-r[0]
    dk = np.pi/(n_pts*dr)
    k = np.arange(0,n_pts * dk, dk)+0.5*dk
    # Calculate c(r)
    c_r = np.empty_like(r)
    gamma_r = np.zeros_like(r)
    c_k = np.empty_like(r)
    gamma_k = np.empty_like(r)
    gamma_r_old = np.empty_like(r)
    errorlast = np.inf
    while n_iter < max_iter:
        n_iter +=1
        # Using closure relation
        c_r[:] = closure_relation(1.0/kBT,u(r,params),gamma_r,model)
        # First, make the FT of c(r)
        constant = 2 * np.pi * dr / k
        transform = dst( c_r * r, type=1)
        c_k[:] = constant * transform
        gamma_r_old[:] = gamma_r
        # Calculate gamma(k)
        gamma_k[:] = rho*c_k**2/(1 - rho*c_k)
        # Calculate IFT of gamma(k)
        constant = n_pts * dk / (4 * np.pi**2 * (n_pts + 1) * r)
        transform = idst(gamma_k * k, type=1)
        gamma_r[:] = constant * transform
        # Calculate the error
        error = np.linalg.norm(gamma_r-gamma_r_old)
        #updating alpha value
        if errorlast > error:  alpha = min(0.9,alpha*finc)
        else: alpha = max(1.0e-3,alpha*fdec)
        if error < atol:
            break
        gamma_r[:] = (1-alpha)*gamma_r_old + alpha*gamma_r
        errorlast = error
    # return the g_r and the c_r
    return [gamma_r+c_r+1,c_r]

def rdf(r,rho,kBT=1.0,u=uhs,params=1.0,model='PY'):
    [g_r,c_r] = PyOZ(r,rho,kBT,u,params,model)
    return g_r

def dcf(r,rho,kBT=1.0,u=uhs,params=1.0,model='PY'):
    [g_r,c_r] = PyOZ(r,rho,kBT,u,params,model)
    return c_r