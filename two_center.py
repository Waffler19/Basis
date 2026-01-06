from Basis import Basis
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import math
import numpy as np
class TwoCenter:
    def __init__(self, orb1='1s', orb2='1s', x1=0, x2=2, Z1=1, Z2=1, bond_type='Bonding', slide=True):
        self.x1 = x1
        self.x2 = x2
        self.dx = x2-x1
        self.orb1 = orb1
        self.orb2 = orb2
        self.Z1 = Z1
        self.Z2 = Z2
        self.bond_type = bond_type
        self.slide = slide
    def fit(self, n_gaussians, rmax=5):
        fit1 = Basis(orbital = self.orb1, Z=self.Z1, rmax=rmax, r_inc=0.001, orbital_type='STO', plot_type='None')
        fit2 = Basis(orbital = self.orb2, Z=self.Z2, rmax=rmax, r_inc=0.001, orbital_type='STO', plot_type='None')
        fit1.fit(number_of_gaussians=n_gaussians)
        fit2.fit(number_of_gaussians=n_gaussians)
        self.c1 = fit1.coeffs
        self.a1 = fit1.alphas
        self.c2 = fit2.coeffs
        self.a2 = fit2.alphas

        self.r = np.linspace(-rmax, rmax, num=int(2*rmax/0.001))
        self.p1 = self.model_2c(self.r, fit1.coeffs, fit1.alphas, self.x1)
        self.p2 = self.model_2c(self.r, fit2.coeffs, fit2.alphas, self.x2)
        ovlp = self.ovlp(self.r, self.p1, self.p2)
        self.plot(self.r, self.p1, self.p2, ovlp)
    def update(self, val):
        self.x2 = self.slider.val
        self.p1 = self.model_2c(self.r, self.c1, self.a1, self.x1)
        self.p2 = self.model_2c(self.r, self.c2, self.a2, self.x2)
        ovlp = self.ovlp(self.r, self.p1,self.p2)
        if self.bond_type == 'Bonding':
            self.line_sum.set_ydata(0.5*self.p1 + 0.5*self.p2)
        elif self.bond_type == 'Anti-Bonding':
            self.line_sum.set_ydata(0.5*self.p1 - 0.5*self.p2)
        self.line1.set_ydata(self.p1)
        self.line2.set_ydata(self.p2)
        self.text.set_text(rf'$\langle \psi_{{\mathrm{{i}}}} | \psi_{{\mathrm{{j}}}} \rangle = {ovlp:.3f}$')
        self.ax.relim()        
        self.ax.autoscale_view()
        self.fig.canvas.draw_idle()

    def model_2c(self, r, coeffs, alphas, x):
        tot = 0
        for ci, ai in zip(coeffs, alphas):
            tot += ci*((2*ai/math.pi)**(0.75))*math.e**(-ai*((r-x)**2))
        return tot
    def ovlp(self, r, p1, p2):
        N = np.trapz(p1**2, r)
        M = np.trapz(p2**2, r)
        p1 = p1/np.sqrt(N)
        p2 = p2/np.sqrt(M)
        ovlp = np.trapz(p1*p2, r)
        return abs(ovlp)
    def plot(self, r, p1, p2, ovlp):
        self.fig, self.ax = plt.subplots()
        self.line1, = self.ax.plot(r, p1, label='Orbital 1')
        self.line2, = self.ax.plot(r, p2, label='Orbital 2')
        self.text = self.ax.text(0.4 * max(r),  0.6 * max(p1), rf'$\langle \psi_{{\mathrm{{i}}}} | \psi_{{\mathrm{{j}}}} \rangle = {ovlp:.3f}$', fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        if self.bond_type == 'Bonding':
            self.line_sum, = self.ax.plot(r, (0.5*p1+0.5*p2), label='Bonding', linestyle='--')
        if self.bond_type == 'Anti-Bonding':
            self.line_sum, = self.ax.plot(r, (0.5*p1-0.5*p2), label='Anti-Bonding', linestyle='--')
        plt.legend()
        if self.slide == True:
            ax_slider = plt.axes([0.925, 0.1, 0.03, 0.8])  
            self.slider = Slider(ax_slider, 'Bond Length', -10.0, 10.0, valinit=self.x2, orientation='vertical')
            self.slider.on_changed(self.update)
        plt.xlabel('Distance (Angstroms)')
        plt.show()
    def nuc_rep(self):
        return self.Z1*self.Z2/(self.dx)
    def kinetic_1e(self):
        return 1
    def potential_1e(self):
        return 2
    def coulomb_2e(self):
        return 2
    def exchange_2e(self):
        return 2

test = TwoCenter('1s', '1s', x1=0, x2=0, Z1=1, Z2=1, bond_type='Anti-Bonding')
test.fit(3, rmax=10)
