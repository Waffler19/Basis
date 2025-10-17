import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math

class Basis():
    def __init__(self, orbital, number_of_gaussians=3, rmax=5, r_inc = 0.05, p_thresh = 0.2, plot_type='radial', stats=False):
        '''Create A Basis Set Object
        Kwargs:
            orbital : String
                Orbital name currently supporting [1s, 2s, 2pz, 3s, 3pz, 3dz, 4s]
            number_of_gaussians : integer
                Number of Gaussians to Include in LC : [1, 2, 3, 6]
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
        self.QN = orbital[0]
        self.number_of_gaussians = number_of_gaussians
        self.rmax =rmax
        self.r_inc = r_inc
        self.p_thresh = p_thresh
        self.plot_type = plot_type
        self.stats = stats
        self.log = []
    
    def kernel(self):
        '''Kernel ==> Calculate STO-NG approximation for Hydrogen Atom Orbitals'''
        r, psi, psi_2, psi_rad = self.generate_correct_function(self.orbital) #generate H atom function'
        #generate gaussian linear combination using 
        if self.number_of_gaussians == 1: 
            popt, pcov = curve_fit(self.model_function_1g, r, psi, p0=[1.2, 0.4])
            p = self.model_function_1g(r, *popt)
        if self.number_of_gaussians == 2:
            popt, pcov = curve_fit(self.model_function_2g, r, psi, p0=[1.2, 0.4, 1.3, 0.33])
            p = self.model_function_2g(r, *popt)
        if self.number_of_gaussians == 3:
            popt, pcov = curve_fit(self.model_function_3g, r, psi, p0=[1.2, 0.4, 1.3, 0.33, 0.33, 0.33])
            p = self.model_function_3g(r, *popt)

        if self.number_of_gaussians == 6:
            popt, pcov = curve_fit(self.model_function_6g, r, psi, p0=[1,1,1,1,1,1])
            p = self.model_function_6g(r, *popt)
        p, psi = self.normalize(r, p, psi) #Normalize p, psi
        overlap_int = self.overlap(r, p, psi) #Compute int(pxpsi dr)
        if self.stats:
            expect_psi, expect_p = self.expectation_r(r,p,psi)
            f = self.convert_to_functions(popt)
            print(f)
        self.psi = psi
        self.p = p
        self.r = r
        if self.plot_type == 'Polar':
            self.to_polar(r, p, psi)
        elif self.plot_type=='Radial':
            if self.stats:
                self.plot(r, p, psi, overlap_int, expect_psi, expect_p)
            else:
                self.plot(r, p, psi, overlap_int)
        else:
            return
    
    def model_function_1g(self, r, alphai, ci):
        '''Model STO-1G'''
        return ci*math.e**(-alphai*(r**2))
    def model_function_2g(self, r, alphai, ci, cj, alphaj):
        '''Model STO-2G'''
        return ci*math.e**(-alphai*(r**2)) + cj*math.e**(-alphaj*(r**2))
    def model_function_3g(self, r, alphai, ci, alphaj, cj, alphak, ck):
        '''Model STO-3G'''
        return ci*math.e**(-alphai*(r**2)) + cj*math.e**(-alphaj*(r**2)) + ck*math.e**(-alphak*(r**2))
    def model_function_6g(self, r, ci, cj, ck, cl, cn, cm):
        '''Model STO-6G'''
        # Exponents from literature b/c otherwize to complicated to do with least squares. 
        return ci*math.e**(-0.0483*(r**2)) + cj*math.e**(-0.1677*(r**2)) + ck*math.e**(-0.4505*(r**2))+ cl*math.e**(-0.9614*(r**2)) + cm*math.e**(-1.7361*(r**2)) + cn*math.e**(-3.3243*(r**2))
    
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
            self.plot_polar(r=r_thresh_psi, theta = theta_thresh_psi, r_p=r_thresh_p, theta_p=theta_thresh, p=p_thresh, psi=psi_thresh)

    def generate_correct_function(self, orbital, r=0):
        '''Generate Correct Hydrogen Atom Functions (Spherical Harmonics for Different Orbitals)'''
        r_list = [] #list of radius for graphing
        psi_list = [] #psi for graphing
        psi_2_list = [] #psi**2 for prob
        psi_2_r_2 = [] #radial probability 
        a0 = 0.529 #bohr radius in A
        while r < self.rmax:
            if orbital == '1s':
                N = (1/math.sqrt(math.pi)) * a0**(-3/2)
                psi = N*math.e**(-r/(a0))
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

    def normalize(self,r,  p, psi):
        normalize_psi = math.sqrt(1/np.trapz(psi*psi, r))
        psi = psi*normalize_psi
        normalize_p = math.sqrt(1/np.trapz(p*p, r))
        p = p*normalize_p
        return p, psi
    
    def overlap(self, r, p, psi):
        return np.trapz(p*psi, r)
    
    def expectation_r(self, r, p, psi):
        return np.trapz(r*psi**2, r), np.trapz(r*p**2, r)
    
    def plot(self, r, p, psi, overlap, psi_e=0, p_e=0):
            '''Plot in R vs Psi Space'''
            plt.plot(r,psi, label='Hydrogen Atom Wavefunction: ' + self.orbital)
            plt.plot(r, p, 'g--', label='STO-'+str(self.number_of_gaussians)+'g')
            plt.text(0.6 * max(r),  0.8 * max(psi), rf'$\langle \psi_{{\mathrm{{GTO}}}} | \psi_{{\mathrm{{STO}}}} \rangle = {overlap:.3f}$', fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
            if self.stats:
                plt.axvline(x=psi_e, linestyle='-')
                plt.axvline(x=p_e, linestyle='--', color='green')
            plt.xlabel('r (Angstrom)', fontsize=12)
            plt.ylabel(rf'$\psi(r)$', fontsize=15)
            plt.legend()
            plt.show()

    def plot_polar(self, r=None, theta=None, theta_p=None, r_p=None, p=None, psi=None):
        '''Plot in R vs Theta Space with Psi as Shading'''
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.scatter(theta, r, c=psi, label='Hydrogen Atom Wavefunction: '+self.orbital, alpha=0.5)
        plt.legend()
        fig1, ax1 = plt.subplots(subplot_kw={'projection': 'polar'})
        ax1.scatter(theta_p, r_p, c=p, label='STO-'+str(self.number_of_gaussians)+'g', alpha=0.5)
        plt.legend()
        plt.show()

k = Basis('3pz', number_of_gaussians=1, p_thresh=0.01, plot_type='Radial', r_inc=0.05, rmax=6, stats=True)
k.kernel()
