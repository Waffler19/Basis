# TWO CENTER INTEGRALS 
Normalized 1D Gaussians of the form:
$\psi_i = c_ie^{-\alpha_i(x-a)^{2}}$
$\psi_j = c_je^{-\alpha_j(x-b)^{2}}$

Center 1 at coordinate x=u1, Center 2 at coordinate x=u2
## 1. Overlap Integral 

$\langle \psi_i \psi_j \rangle = \int c_ic_je^{-\alpha_i(x-u_1)^{2}}e^{-\alpha_j(x-u_2)^{2}}dx$

$\langle \psi_i \psi_j \rangle = \int c_ic_je^{-\alpha_i(x-u_1)^{2}-\alpha_j(x-u_2)^{2}}dx$

$\langle \psi_i \psi_j \rangle = \int c_ic_je^{-\alpha_i(x^{2}-2xu_1+u_1^{2})-\alpha_j(x^{2}-2xu_2+u_2^{2})}dx$

$\langle \psi_i \psi_j \rangle = \int c_ic_je^{-x^{2}(\alpha_i+\alpha_j) + 2x(\alpha_iu_1+\alpha_ju_2) - u_1^{2}\alpha_i-u_2^{2}\alpha_j}dx$

We can simplify by completing the square in the exponent: $-x^{2}(\alpha_i+\alpha_j) + 2x(\alpha_iu_1+\alpha_ju_2) - u_1^{2}\alpha_i-u_2^{2}\alpha_j$

$u_1^{2}\alpha_i + u_2^{2}\alpha_j = -x^{2}(\alpha_i+\alpha_j) + 2x(\alpha_iu_1+\alpha_ju_2)$

$\frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{-(\alpha_i+\alpha_j)} = x^{2} + \frac{2x(\alpha_iu_1+\alpha_ju_2)}{-(\alpha_i+\alpha_j)}$

$\frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{-(\alpha_i+\alpha_j)} + (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}= x^{2} + \frac{2x(\alpha_iu_1+\alpha_ju_2)}{-(\alpha_i+\alpha_j)} + (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}$


$\frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{-(\alpha_i+\alpha_j)} + (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}= (x + \frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}$

$0 = (x + \frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2} + \frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{-(\alpha_i+\alpha_j)} + (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}$

Now plug in our expression to the integral: 

$\langle \psi_i \psi_j \rangle = \int c_ic_je^{(x + \frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2} + \frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{-(\alpha_i+\alpha_j)} + (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}}dx$

$\langle \psi_i \psi_j \rangle = c_ic_j\int e^{(x + \frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}}e^{\frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{-(\alpha_i+\alpha_j)} + (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}}dx$

$\langle \psi_i \psi_j \rangle = c_ic_je^{\frac{u_1^{2}\alpha_i + u_2^{2}\alpha_j}{(\alpha_i+\alpha_j)} - (\frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}} \int e^{-(x + \frac{\alpha_iu_1+\alpha_ju_2}{-2(\alpha_i+\alpha_j)})^{2}}dx$

Now our integral is solvable!
