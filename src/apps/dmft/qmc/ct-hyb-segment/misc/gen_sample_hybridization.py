"""
 A short script to generate a hybridization function Delta(tau) and Delta(omega) out of the given DOS(omega).
 A. E. Antipov (2015)
"""

import numpy as np
from scipy import integrate

def dos_function(omega):
    """ put here your dos"""
    nu = 10.0
    w0 = 10.0
    dos_f = lambda x : 1./(1 + np.exp(nu*(x-w0))) * 1./(1 + np.exp(-(nu*(x+w0)))) 
    return dos_f(omega)
    

def generate_delta(params):
    print params
    D = params["D"]
    beta = params["beta"]
    npts = params["npts"]
    print "DOS half-bandwidth : ", D
    Gamma = params["gamma"]
    dos_prec = 0.1 

    omega_grid = np.arange(-2*D,2*D,dos_prec)
    dos_vals = dos_function(omega_grid) 
    data = np.vstack([omega_grid,dos_vals]).transpose()
    np.savetxt("dos.dat",data)

    fermi = lambda w : 1. / (1.+np.exp(beta*w))
    delta_wt = lambda tau, w : -fermi(w) * np.exp(tau*w) * dos_function(w) * Gamma
    delta_f = lambda tau : integrate.quad(lambda w: -fermi(w) * np.exp(tau*w) * dos_function(w) * Gamma, -2*D, 2*D) 

    tau_grid = np.linspace(0,beta,npts+1)
    delta_vals = np.array([delta_f(x)[0] for x in tau_grid])
    data_out = np.vstack([range(npts+1), delta_vals, delta_vals]) 
    fname = "delta_tau.dat"
    np.savetxt("delta_tau.dat", data_out.transpose())
    print "Saved", fname

    kramers_kronig_imag = lambda z : integrate.quad(lambda w: np.imag(dos_function(w) / (1j*z - w)), -2*D, 2*D)
    kramers_kronig_real = lambda z : integrate.quad(lambda w: np.real(dos_function(w) / (1j*z - w)), -2*D, 2*D)
    matsubara_grid = (2*np.arange(0, npts, 1, dtype=np.float) + 1)*np.pi/beta
    delta_iw = np.array([[x, kramers_kronig_real(x)[0], kramers_kronig_imag(x)[0]] for x in matsubara_grid])
    fname = "delta_iw.dat"
    np.savetxt(fname, delta_iw)
    print "Saved", fname
    
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='gen imag delta')
    parser.add_argument('--beta', help='Inverse temperature', type=float, default = 5)
    parser.add_argument('--npts', help='Number of points on tau and Matsubara grid', type=int, default = 100)
    parser.add_argument('--D', help='half bandwidth parameter (grid is between [-2D; 2D])', type=float, default = 5)
    parser.add_argument('--gamma', help='coupling strength', type=float, default = 1)
    args = parser.parse_args()
    generate_delta(vars(args))


