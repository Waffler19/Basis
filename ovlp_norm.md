# Derivations

## I) Deriving the Overlap Integrals $\langle \psi_i | \psi_j* \rangle$

Overlap Integrals can be approximated numerically using various riemann sums like techniques. However, analytical solutions are more computationally viable and can be easily modified to incorporate different types of functions. Here I will briefly present a short derivation for the overlap integral of two entirely radial gaussian functions ($\psi_i$, $\psi_j$) with coefficients ($c_i$, $c_j$) and exponential paramters ($\alpha_i$, $\alpha_j$) on the same center. 

### 1. Setup the Integral in 3D Polar Coordinates with Respect to radius R

$S_{ij} = \langle \psi_i | \psi_j* \rangle = \int_{0}^{2\pi}\,\int_{0}^{\pi}\,\int_{0}^{\infty} \psi_i\psi_j \, r^{2}sin(\theta)dr d\theta d\phi$

### 1.5 Let's Deal with Those Pesky Angle Dependencies first

$S_{ij} =  \int_{0}^{2\pi}\ d\phi \int_{0}^{\pi}\ sin(\theta)d\theta \int_{0}^{\infty} \psi_i\psi_j \, r^{2}dr$

$S_{ij} = \left. \phi \right|_0^{2\pi} \left. * -cos(\theta)\right|_0^{\pi}  \int_{0}^{\infty} \psi_i\psi_j \, r^{2}dr$

$S_{ij} = 2\pi* 2\int_{0}^{\infty} \psi_i\psi_j \, r^{2}dr = 4\pi \int_{0}^{\infty} \psi_i\psi_j \, r^{2}dr$

Now let's plug in our Gaussians:
$\psi_i = c_ie^{-\alpha_ir^{2}}$
$\psi_j = c_je^{-\alpha_jr^{2}}$

### 2. Plugging in and Rearrange:

$S_{ij} = 4\pi\int_{0}^{\infty} c_ie^{-\alpha_ir^{2}} * c_je^{-\alpha_jr^{2}} \, r^{2}dr$

$S_{ij} = 4\pi c_ic_j \int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}} \, r^{2}dr$

### 3. Apply Integration by parts where 

$dv = re^{-(\alpha_i+\alpha_j)r^{2}}$ 

$v = \int \,dv = -e^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j)$ 

$u = r$

$du = 1$

$S_{ij} = 4\pi c_ic_j\int_{0}^{\infty} u \, dv = 4\pi c_ic_ju*v - \int_{0}^{\infty} v du = 4\pi c_ic_j(-re^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j) + \int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j))\, dr$

### 4. Take the Limit of the Function as R approaches infinity and Zero

$\lim_{r \to \infty} (S_{ij}) = 4\pi c_ic_j(-\infty e^{-(\alpha_i+\alpha_j)\infty^{2}}/2(\alpha_i+\alpha_j) + \int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j))\, dr$ 

$\lim_{r \to 0} (S_{ij}) = 4\pi c_ic_j(-0 e^{-(\alpha_i+\alpha_j)0^{2}}/2(\alpha_i+\alpha_j) + \int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j))\, dr$ 

First Term goes to 0 in both cases so we are left with: 

$S_{ij} = 0 + 4\pi c_ic_j\int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j) \, dr$ 

### 5. Rearrange and Simplify More to Apply an Identity

$S_{ij} = 4\pi c_ic_j/2(\alpha_i+\alpha_j)\int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}} \, dr$

u-sub with 

$u = \sqrt{\alpha_i+\alpha_j}r$

$du = \sqrt{\alpha_i+\alpha_j}$

$S_{ij} = 4\pi c_ic_j/2(\alpha_i+\alpha_j)^{3/2}\int_{0}^{\infty} e^{-u^{2}} \, dr$

### 6. Apply the Euler-Poissen Integral: $\int_{0}^{\infty} e^{-u^{2}}\, dr = \sqrt{\pi}/2$

$S_{ij} = 4\pi c_ic_j/2(\alpha_i+\alpha_j)^{3/2} * \sqrt{\pi}/2 = c_ic_j\pi^{3/2}/(\alpha_i+\alpha_j)^{3/2}$ 

$S_{ij} =  c_ic_j (\frac {\pi}{(\alpha_i+\alpha_j)})^{3/2}$

### 7. Applying Our Result to Normalization: 

$S_{ij} =  c_i^{2}(\frac {\pi}{(2\alpha_i)})^{3/2}$

## We can normalize and Solve for N:

$S_{ij} =  c_i^{2}(\frac {\pi}{(2\alpha_i)})^{3/2}$
We want wavefunctions of the form: $\Psi_{gto} = Nc_ie^{-\alpha_i r^{2}}$

Where $\langle \psi_i | \psi_i* \rangle = S_{ii} = 1$

Using Our Previous Result and Adding an N constant to $\psi_i$: 

$S_{ij} = 1 = N^{2}c_i^{2}(\frac {\pi}{(2\alpha_i)})^{3/2}$

$N^{2} =\frac{(\frac{2\alpha_i}{\pi})^{3/2}}{c_i^{2}}$
$N = \frac{(\frac{2\alpha_i}{\pi})^{3/4}}{c_i}$

### 8. Therefore to be normalized, each specific gaussian in a contraction should have the form: 

## $\psi_i = {(\frac{2\alpha_i}{\pi})^{3/4}}c_ie^{-\alpha_i r^{2}}$

## II) Normalizing the STO's
https://en.wikipedia.org/wiki/Slater-type_orbital

Slater Type Orbitals Have the General Form: 

$R(r) = Nr^{n-1}e^{-\gamma r}$

where n is the principle quantum number, r is the radial distance, and $\gamma$ is the nuclear charge

### 1. Taking the General Overlap Integral in Spherical Coordiantes:

$\langle \psi_{sto} | \psi_{sto}^{*} \rangle = \int_{0}^{2\pi}\,\int_{0}^{\pi}\,\int_{0}^{\infty} \psi_{sto}\psi_{sto}^{*} \, r^{2}sin(\theta)dr d\theta d\phi = 1$

$\langle \psi_{sto} | \psi_{sto}^{*} \rangle = \int_{0}^{2\pi}\,\int_{0}^{\pi}\,\int_{0}^{\infty} N^{2}r^{2(n-1)}e^{-2\gamma r}\, r^{2}sin(\theta)dr d\theta d\phi$


$\langle \psi_{sto} | \psi_{sto}^{*} \rangle = N^{2}\int_{0}^{2\pi}\ d\phi,\int_{0}^{\pi}\ sin(\theta) d\theta,\int_{0}^{\infty} r^{2n-2+2}e^{-2\gamma r}\, dr $

$\langle \psi_{sto} | \psi_{sto}^{*} \rangle = 4\pi N^{2}\int_{0}^{\infty} r^{2n}e^{-2\gamma r}\, dr $

### All Spherical Harmonics Additions to R(r) Have a 1/2 sqrt(pi) Piece so we can simplify to:

$\langle \psi_{sto} | \psi_{sto}^{*} \rangle = N^{2}\int_{0}^{\infty} r^{2n}e^{-2\gamma r}\, dr $


### 2. Using a Known Integral 

$\int_{0}^{\infty} x^{n}e^{-ax} \, dx = \frac{n!}{a^{n+1}}$

$n = 2n$ $a = 2\gamma $

### 3. Simplify
$\langle \psi_{sto} | \psi_{sto}^{*} \rangle = \frac{N^{2} (2n)!}{(2\gamma)^{2n+1}} = 1$

$N = \sqrt{\frac{(2\gamma)^{2n+1}}{4\pi(2n)!}} = \sqrt{\frac{(2\gamma)^{2n}(2\gamma)^1}{(2n)!}} = (2\gamma)^n \sqrt{\frac{(2\gamma)}{(2n)!}}$
