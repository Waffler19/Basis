import numpy as np
from integrator import Integrator
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import math

class Basis():
    def __init__(self, orbital, Z=1, verbose=1, rmax=5, r_inc = 0.05, p_thresh = 0.2, plot_type='Radial', orbital_type='STO'):
        '''Create A Basis Set Object
        Kwargs:
            orbital : String
                Orbital name currently supporting [1s, 2s, 2pz, 3s, 3pz, 3dz, 4s]

            rmax : float
                Maximum Radius to use when constructing orbitals

            r_inc : float
                Step size for moving along the r/theta axis 

            p_thresh : float
                probabillity threshold to include points on output graph in polar coordinates

            plot_type : string
                Radial - generate plot of p vs r
                Polar - generate plot of r vs theta with p as a colorscale
            '''
        self.orbital = orbital
        self.verbose = verbose
        self.Z = Z
        self.QN = orbital[0]
        self.rmax =rmax
        self.r_inc = r_inc
        self.p_thresh = p_thresh
        self.plot_type = plot_type
        self.log = []
        self.orbital_type = orbital_type
        self.KE = 0
        self.PE = 0
        self.E = 0
        self.coeffs = []
        self.alphas = []
        #self.r, self.psi, self.psi_2, self.psi_rad = self.generate_correct_function(self.orbital) #generate H atom function'

    def fit(self, number_of_gaussians=3, stats=False, numerical=False, plot=True, plot_ind=False, opt='psi'):
        '''Kernel ==> Calculate STO-NG approximation for Hydrogen Atom Orbitals'''
        self.opt = opt
        self.number_of_gaussians = number_of_gaussians
        self.stats = stats
        self.numerical = numerical
        self.plot_Ind = plot_ind
        r, psi, psi_2, psi_rad = self.generate_correct_function(self.orbital) #generate H atom function'
        #generate gaussian linear combination using 
        if self.number_of_gaussians == 0:
            self.plot(r, p='None', psi=psi)
            return
        if self.orbital == '1s':
            guess = np.array([3.42525, 0.154, 0.62391, 0.53, 0.16886, 0.44])*self.Z**(1.5)
            guess = list(guess)
        guess = np.array([0.1, 0.4, 0.5, -0.1, 1, 0.75])*self.Z**(1.5)
        guess = list(guess)
        if self.number_of_gaussians == 1: 
            if self.opt == 'ovlp': 
                min = minimize(self.overlap_loss, guess[0:2], args=(r, psi))
                popt = min.x
            elif self.opt == 'psi':
                popt, pcov = curve_fit(self.model_function_1g, r, psi, p0=[0.01, 0.4])
            p = self.model_function_1g(r, *popt)
        if self.number_of_gaussians == 2:
            if self.opt == 'ovlp': 
                min = minimize(self.overlap_loss, guess[0:4], args=(r, psi), bounds=((0.01,10*self.Z**2), (-100,100), (0.0001,10*self.Z**2), (-100,100)))
                popt = min.x
            elif self.opt == 'psi':
                popt, pcov = curve_fit(self.model_function_2g, r, psi, p0=[0.1, 0.4, 0.5, 0.33], bounds=[[0, -100, 0, -100],[10*self.Z**2, 100, 10*self.Z**2, 100]], maxfev=25000)
            p = self.model_function_2g(r, *popt)
        if self.number_of_gaussians == 3:
            if self.opt == 'ovlp':
                min = minimize(self.overlap_loss, guess[0:6], args=(r, psi), bounds=((0.01,10*self.Z**2), (-100,100), (0.0001,10*self.Z**2), (-100,100), (0.0001,10*self.Z**2), (-100,100)), options={'maxiter':25000})
                popt = min.x
            elif self.opt == 'psi': 
                popt, pcov = curve_fit(self.nmodel_function_3g,r, psi, p0=[0.1, 0.154, 0.05, 0.53, 0.2, 0.44],bounds=[[0, -100, 0, -100, 0, -100],[10*self.Z**2, 100, 10*self.Z**2, 100, 10*self.Z**2, 100]], maxfev=25000)
            p = self.nmodel_function_3g(r, *popt)

        if self.number_of_gaussians == 6:
            if self.opt == 'ovlp':
                min = minimize(self.overlap_loss, [0.0483, 1,0.1677,1,0.4505,1,0.9614,1,1.7361,1,3.3243,1], bounds=((0.01,10*self.Z**2), (-100,100), (0.0001,25*self.Z**2), (-100,100), (0.0001,10*self.Z**2), (-100,100), (0.01,10*self.Z**2), (-100,100), (0.0001,25), (-100,100), (0.0001,10*self.Z**2), (-100,100)), args=(r, psi, True))
                min2 = minimize(self.overlap_loss, [0.0483, 1,0.1677,1,0.4505,1,0.9614,1,1.7361,1,3.3243,1], bounds=((0.01,10*self.Z**2), (-100,100), (0.0001,25*self.Z**2), (-100,100), (0.0001,10*self.Z**2), (-100,100), (0.01,10*self.Z**2), (-100,100), (0.0001,25), (-100,100), (0.0001,10*self.Z**2), (-100,100)), args=(r, psi, None, min.x[0::2]))
                #min = minimize(self.overlap_loss, [0.0483, 1,0.1677,1,0.4505,1,0.9614,1,1.7361,1,3.3243,1], args=(r, psi), fixed=True)
                popt = []
                for ai, ci in zip(min.x[0::2], min2.x[1::2]):
                    popt.append(ai)
                    popt.append(ci)
                popt = np.array(popt)

            elif self.opt == 'psi': 
                popt, pcov = curve_fit(self.model_function_6g_full, r, psi, p0=[0.0483, 1,0.1677,1,0.4505,1,0.9614,1,1.7361,1,3.3243,1], bounds=[[0, -100, 0, -100, 0, -100, 0, -100, 0, -100, 0, -100],[10*self.Z**2, 100, 10*self.Z**2, 100, 10*self.Z**2, 100, 10*self.Z**2, 100, 10*self.Z**2, 100, 10*self.Z**2, 100]], maxfev=25000)
            
            p = self.model_function_6g_full(r, *popt)

        if self.numerical == True:
            p, psi, Np, Npsi = self.normalize(r, p, psi) #Normalize p, psi
            overlap_int = self.overlap(r, p, p) #Compute int(pxpsi dr)
            overlap_int2 = self.overlap(r, psi, psi) #Compute int(pxpsi dr)
            overlap_int = self.overlap(r, p, psi) #Compute int(pxpsi dr)
            
        if self.numerical == False:
            pnon, psi, Np, Npsi = self.normalize(r, p, psi) 
            k = Integrator(self.number_of_gaussians, popt, self.rmax, 'ana')
            k.kernel_normalize()
            Np = k.normalized
            k.coeffs = k.coeffs * math.sqrt(Np)
            p = self.model(r, popt[1::2]*math.sqrt(Np), popt[0::2])
            self.KE, self.PE, self.E = k.kernel_KE_PE()
            if self.verbose > 1:
                print('1e kinetic energy eV', self.KE*27.2114)
                print('1e potential energy eV', self.PE*27.2114)
                print('1e total energy eV', self.E*27.2114)


            
            self_ovlp = self.overlap(r, p, p) #Check normalization
            self_ovlp2 = self.overlap(r, psi, psi)
            overlap_int = self.overlap(r, p, psi)
        
        if self.stats:
            expect_psi, expect_p = self.expectation_r(r,p,psi)
            f = self.convert_to_functions(popt)
        
        if self.plot_Ind:
            count = 1
            plots = []
            for i in popt[0::2]:
                p_new = self.model_function_1g(r, i, math.sqrt(Np)*popt[count])
                count += 2
                plots.append(p_new)
        if plot:
            if self.plot_type in ['Radial', 'Polar']:
                if self.plot_type == 'Polar':
                    self.to_polar(r, p, psi)
                if self.plot_type=='Radial':
                    if self.stats:
                        self.plot(r, p, psi, overlap_int, expect_psi, expect_p)
                    elif self.plot_Ind:
                        self.plot(r, p, psi, overlap_int, inds = [plots])
                    else:
                        self.plot(r, p, psi, overlap_int)
                else:
                    return
        self.r = r
        return
        
    def overlap_loss(self, params, r, psi_sto, ci_fixed=None, alpha_fixed=None):            
        count = 1
        counter = 0
        psi_gto = 0
        for i in params[0::2]:
             alphai = i
             ci = params[count]
             if ci_fixed is not None:
                ci = 1
             elif alpha_fixed is not None:
                alphai = alpha_fixed[counter]
             counter += 1
             count += 2
             psi_gto += ci*((2*alphai/math.pi)**(0.75))*math.e**(-alphai*(r**2))
        top = np.trapz(psi_sto*psi_gto*r**2, r)
        bottom = math.sqrt(np.trapz(psi_gto*psi_gto*r**2)*np.trapz(psi_sto*psi_sto*r**2))
        S = top/bottom
        return -S

    def model_function_1g(self, r, alphai, ci):
        '''Model STO-1G'''
        return ci*((2*alphai/math.pi)**(0.75))*math.e**(-alphai*(r**2))
    def model_function_2g(self, r, alphai, ci, alphaj, cj):
        '''Model STO-2G'''
        return ci*((2*alphai/math.pi)**(0.75))*math.e**(-alphai*(r**2)) + cj*((2*alphaj/math.pi)**(0.75))*math.e**(-alphaj*(r**2))
    def nmodel_function_3g(self, r, alphai, ci, alphaj, cj, alphak, ck):
        '''Model STO-3G'''
        return ci*((2*alphai/math.pi)**(0.75))*math.e**(-alphai*(r**2)) + cj*((2*alphaj/math.pi)**(0.75))*math.e**(-alphaj*(r**2)) + ck*((2*alphak/math.pi)**(0.75))*math.e**(-alphak*(r**2))        
    def model_function_3g(self, r, alphai, ci, alphaj, cj, alphak, ck):
        '''Model STO-3G'''
        return ci*math.e**(-alphai*(r**2)) + cj*math.e**(-alphaj*(r**2)) + ck*math.e**(-alphak*(r**2))
    def model_function_6g(self, r, ci, cj, ck, cl, cn, cm):
        '''Model STO-6G'''
        # Exponents from literature b/c otherwize to complicated to do with least squares. 
        return ci*((2*0.0483/math.pi)**(0.75))*math.e**(-0.0483*(r**2)) + cj*((2*0.1677/math.pi)**(0.75))*math.e**(-0.1677*(r**2)) + ck*((2*0.4505/math.pi)**(0.75))*math.e**(-0.4505*(r**2))+ cl*((2*0.9614/math.pi)**(0.75))*math.e**(-0.9614*(r**2)) + cm*((2*1.7361/math.pi)**(0.75))*math.e**(-1.7361*(r**2)) + cn*((2*3.3243/math.pi)**(0.75))*math.e**(-3.3243*(r**2))
    def model_function_6g_full(self, r, ai, ci, aj, cj, ak, ck, al, cl, an, cn, am, cm):
        '''Model STO-6G'''
        # Exponents from literature b/c otherwize to complicated to do with least squares. 
        return ci*((2*ai/math.pi)**(0.75))*math.e**(-ai*(r**2)) + cj*((2*aj/math.pi)**(0.75))*math.e**(-aj*(r**2)) + ck*((2*ak/math.pi)**(0.75))*math.e**(-ak*(r**2))+ cl*((2*al/math.pi)**(0.75))*math.e**(-al*(r**2)) + cm*((2*am/math.pi)**(0.75))*math.e**(-am*(r**2)) + cn*((2*an/math.pi)**(0.75))*math.e**(-an*(r**2))
    def model(self, r, coeffs, alphas=None):
        self.coeffs = coeffs
        self.alphas = alphas
        if self.number_of_gaussians == 6:
            p = self.model_function_6g_full(r, alphas[0], coeffs[0], alphas[1], coeffs[1], alphas[2], coeffs[2], alphas[3], coeffs[3], alphas[4], coeffs[4], alphas[5], coeffs[5])
        if self.number_of_gaussians == 1:
            p = self.model_function_1g(r, alphas[0], coeffs[0])
        if self.number_of_gaussians == 2:
            p = self.model_function_2g(r, alphas[0], coeffs[0], alphas[1], coeffs[1])
        if self.number_of_gaussians == 3:
            p = self.nmodel_function_3g(r, alphas[0], coeffs[0], alphas[1], coeffs[1], alphas[2], coeffs[2])
        if self.verbose > 1:
            print('coefficients: ', coeffs)
            print('alphas: ', alphas)
        return p
    def to_polar(self, r, p, psi):
            '''Convert to Polar Coordinates to Handle p like orbitals'''
            angles = np.linspace(0,2*math.pi, int(self.rmax/self.r_inc + 1))
            if self.orbital in ['2pz', '3pz']:
                theta = np.cos(angles)
            elif self.orbital in ['3dz']:
                theta = 3*(np.cos(angles)**2)-1
            else:
                theta = angles
            R_sto = psi
            R_g = p
            p = []
            psi = []
            for i,j in zip(R_g,R_sto):
                for k in theta:
                    if self.orbital in ['2pz', '3pz', '3dz']:
                        p.append(i*k)
                        psi.append(j*k)
                    else:
                        p.append(i)
                        psi.append(j)
            p = np.array(p)**2
            psi = np.array(psi)**2
            p_index = np.where(p>self.p_thresh)
            psi_index = np.where(psi>self.p_thresh)
            p_thresh = []
            psi_thresh = []
            r_thresh_p = []
            r_thresh_psi = []
            theta_thresh = []
            theta_thresh_psi = []
            for i in p_index:
                p_thresh.append(i)
                r_thresh_p.append(r[i//(int(self.rmax/self.r_inc + 1))])
                theta_thresh.append(angles[i%(int(self.rmax/self.r_inc + 1))])
            for i in psi_index:
                psi_thresh.append(i)
                r_thresh_psi.append(r[i//(int(self.rmax/self.r_inc + 1))])
                theta_thresh_psi.append(angles[i%(int(self.rmax/self.r_inc + 1))])
            self.plot_polar(fitted=True, r=r_thresh_psi, theta = theta_thresh_psi, r_p=r_thresh_p, theta_p=theta_thresh, p=p_thresh, psi=psi_thresh)

    def generate_correct_function(self, orbital, r=0):
        '''Generate Correct Hydrogen Atom Functions (Spherical Harmonics for Different Orbitals)'''
        r_list = [] #list of radius for graphing
        psi_list = [] #psi for graphing
        psi_2_list = [] #psi**2 for prob
        psi_2_r_2 = [] #radial probability 
        a0 = 0.529 #bohr radius in atomic units = 1
        #self.orbital_type = 'STO'
        while r < self.rmax:
            if self.orbital_type == 'Spherical Harmonics':
                if orbital == '1s':
                    N = (1/math.sqrt(math.pi)) * a0**(-3/2)
                    psi = N*math.e**(-r/a0)
                if orbital == '2s':
                    N = (1/(4*math.sqrt(2*math.pi))) * a0**(-3/2)
                    psi = N*(2-(r/a0))*math.e**(-r/(2*a0))
                if orbital == '2pz':
                    N = (1/(4*math.sqrt(2*math.pi))) * a0**(-3/2)    
                    psi = N*(r/a0)*math.e**(-r/(2*a0))
                if orbital == '3s':
                    N = (1/(81*math.sqrt(3*math.pi))) * a0**(-3/2)
                    psi = N*(27-(18*r/a0) + 2*r**2/a0**2)*math.e**(-r/(2*a0))
                if orbital == '3pz':
                    N = (1/(81*math.sqrt(3*math.pi)))* a0**(-3/2)
                    psi = N*(6-(r/a0))*(r/a0)*math.e**(-r/(3*a0))
                if orbital == '3dz':
                    N = (1/(81*math.sqrt(3*math.pi))) * a0**(-3/2)
                    psi = N*(r**2/a0**2)*math.e**(-r/(3*a0))
                if orbital == '4s':
                    N = (1/(81*math.sqrt(3*math.pi))) * a0**(-3/2)
                    psi = N*(r**2/a0**2)*math.e**(-r/(3*a0))
            if self.orbital_type == 'STO':
                N = int(orbital[0])
                Norm = (2*self.Z)**N * ((2*self.Z)/(math.factorial(2*N)))**(1/2)
                #Norm =1 
                psi = Norm*((r**(N-1))*math.e**(-self.Z*r))*(1/2*math.sqrt(math.pi)) #for Angular Part
            r_list.append(r)
            psi_prob = psi**2
            rad_prob = (psi**2) * (r**2)
            psi_list.append(psi)
            psi_2_list.append(psi_prob)
            psi_2_r_2.append(rad_prob)
            r += self.r_inc
        return np.array(r_list), np.array(psi_list), np.array(psi_2_list), np.array(psi_2_r_2)

    def convert_to_functions(self, coefficients):
        '''Convert STO-3G Generated Function to a Readable String'''
        coefficients = coefficients.tolist()
        f = r"$\Psi_{\text{GTO}} = "   
        for i in range(0, len(coefficients), 2):
            c = str(round(coefficients[i],3))
            alpha = str(round(coefficients[i+1], 3))
            print(alpha)
            f+=f"{c}e^{{-{alpha}(r^2)}}"
            if i < len(coefficients)-2:
                f += ' + '
        f += '$'
        print(f)
        return f
    
#Numerical Integrals:
    def normalize(self,r,  p, psi):
        r_s = r**2
        normalize_psi = (1/(4*math.pi*np.trapz(r_s*psi*psi, r)))
        psi = psi*math.sqrt(normalize_psi)
        normalize_p =  (1/((4*math.pi*np.trapz(r_s*p*p, r))))
        p = p*math.sqrt(normalize_p)
        return p, psi, normalize_p, normalize_psi
    
    def overlap(self, r, p, psi):
        r_s = r**2
        #('Overlap', 4*math.pi*np.trapz(r_s*p*psi, r))
        return 4*math.pi*np.trapz(r_s*p*psi, r)
    
    def expectation_r(self, r, p, psi):
        return 4*math.pi*np.trapz((r**3)*psi**2, r), 4*math.pi*np.trapz((r**3)*p**2)

        
    def plot(self, r, p, psi, overlap=0, psi_e=0, p_e=0, inds=[]):
            '''Plot in R vs Psi Space'''
            if type(p) == str:
                plt.plot(r, psi, label='Hydrogen Atom Wavefunction: ' + self.orbital)
                plt.xlabel('r (Angstrom)', fontsize=12)
                plt.ylabel(rf'$\psi(r)$', fontsize=15)
                plt.show()
                return
            plt.plot(r,psi, label='Hydrogen Atom Wavefunction: ' + self.orbital)
            plt.plot(r, p, 'g--', label='STO-'+str(self.number_of_gaussians)+'g')
            plt.text(0.6 * max(r),  0.8 * max(psi), rf'$\langle \psi_{{\mathrm{{GTO}}}} | \psi_{{\mathrm{{STO}}}} \rangle = {overlap:.3f}$', fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            if self.stats:
                plt.axvline(x=psi_e, linestyle='-')
                plt.axvline(x=p_e, linestyle='--', color='green')
            if self.plot_Ind:
                count = 1
                for i in inds[0]:
                    plt.plot(r, i, label='Gaussian'+str(count), alpha=0.75, linestyle=':')
                    count += 1
            plt.xlabel('r (Angstrom)', fontsize=12)
            plt.ylabel(rf'$\psi(r)$', fontsize=15)
            plt.legend()
            plt.show()

    def plot_polar(self, fitted=True, r=None, theta=None, theta_p=None, r_p=None, p=None, psi=None):
        '''Plot in R vs Theta Space with Psi as Shading'''
        if fitted:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.scatter(theta, r, c=psi, label='Hydrogen Atom Wavefunction: '+self.orbital, alpha=0.5)
            plt.legend()
            fig1, ax1 = plt.subplots(subplot_kw={'projection': 'polar'})
            ax1.scatter(theta_p, r_p, c=p, label='STO-'+str(self.number_of_gaussians)+'g', alpha=0.5)
        else:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.scatter(theta, r, c=psi, label='Hydrogen Atom Wavefunction: '+self.orbital, alpha=0.5)
        plt.legend()
        plt.show()

if __name__ == '__main__':
    orbital = Basis('1s', Z=1, verbose=3, rmax=5, r_inc=0.05, p_thresh=0.1, plot_type='Radial', orbital_type='STO')
    orbital.fit(3, plot_ind=True)