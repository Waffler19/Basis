import numpy as np
import math 
class Integrator():

    def __init__(self, g_num, pconv, rmax, mode='erf'):
        self.g_num = g_num
        self.pconv = pconv
        self.rmax = rmax
        self.normalized = 1
        self.mode = mode
        self.coeffs = []
    def kernel_normalize(self):
        #Build Matrix of All Possible Combinations of Functions
        coeffs = np.array(self.pconv[1::2])
        alphas = np.array(self.pconv[0::2])
        for c, a in zip(coeffs, alphas):
           self.coeffs.append(self.init_norm(c, a))
        self.coeffs = np.array(self.coeffs)
        coeff_arr1 = (self.coeffs.reshape(self.g_num, 1))@(self.coeffs.reshape(1, self.g_num))
        alpha_arr = np.add(alphas.reshape(self.g_num, 1), alphas.reshape(1, self.g_num))
        flat_c = coeff_arr1.flatten()
        flat_a = alpha_arr.flatten()
        integral = 0
        for i, j in zip(flat_c, flat_a):
            if self.mode == 'erf':
                integral += self.erf_normalization(i, j)
            if self.mode =='ana':
                integral += self.full_analytical(i, j)
        N = 1/(integral)
        self.normalized = N
        return N
    def init_norm(self, ci, alphai):
        ci = ci*(2*alphai/math.pi)**0.75
        return ci
    def kernel_KE_PE(self):
        alphas = np.array(self.pconv[0::2])
        coeff_arr1 = (self.coeffs.reshape(self.g_num, 1))@(self.coeffs.reshape(1, self.g_num))
        suma_arr = np.add(alphas.reshape(self.g_num, 1), alphas.reshape(1, self.g_num))
        proda_arr = (alphas.reshape(self.g_num, 1))@(alphas.reshape(1, self.g_num))
        flat_c = coeff_arr1.flatten()
        sum_a = suma_arr.flatten()
        prod_a = proda_arr.flatten()
        integral_KE = 0
        integral_PE = 0
        s = 0
        for i, j,k in zip(flat_c, sum_a, prod_a):
            s += self.full_analytical(i, j)
            integral_KE += self.ke(i, j, k)
            integral_PE += self.pe(i, j)
        E = integral_KE/s + integral_PE/s
        return integral_KE/s, integral_PE/s, E
    def erf_normalization(self, ci_cj, alpha_ij):
        '''See ReadMe for Integral Derivation'''
        ai_aj = math.sqrt(alpha_ij)
        term1 = (self.rmax/(-2*alpha_ij))*math.e**(-alpha_ij*(self.rmax**2))
        #Define Constants for ERF Approximation
        const2 = math.sqrt(math.pi)/(4*((alpha_ij)**(1.5)))
        a1 = 0.278393
        a2 = 0.230389
        a3 = 0.000972
        a4 = 0.078108
        erf_denom = 1 + a1*math.sqrt(ai_aj)*self.rmax + a2*ai_aj*self.rmax**2 + a3*(ai_aj**(1.5))*self.rmax**3 + a4*(ai_aj**2)*self.rmax**4
        erf = (1 - 1/(erf_denom**4))
        term2 = const2*erf
        integral = (term1 + term2)*ci_cj
        integral = integral
        return integral
    
    def full_analytical(self, ci_cj, alpha_ij):
        '''See Derivation in ReadMe'''
        integral = ci_cj * (math.pi / alpha_ij)**1.5
        return integral

    def ke(self, ci_cj, sum_a, prod_a):
        h_bar = 1 # in atomic units = 1
        m = 1
        top = 3*ci_cj*prod_a*(math.pi**(1.5))*h_bar**2
        bottom = (sum_a**(5/2))*m
        integral = top/bottom
        return integral
    
    def pe(self, ci_cj, sum_a):
        integral = -(math.pi**(1.5))*ci_cj/(sum_a**0.5)
        return integral
        
##k = Integrator(3, [3.42525, 0.154, 0.62391, 0.53, 0.16886, 0.44], 5, mode='ana')
#N = k.kernel_normalize()
#k.coeffs = k.coeffs * math.sqrt(N)
#print(k.kernel_KE())

#3.42525, 0.154, 0.62391, 0.53, 0.16886, 0.44