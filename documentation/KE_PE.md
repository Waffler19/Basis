## 1e Kinetic and Potential Energy Derivations for Radial Wavefunction on 1-Center of the Form:
## $\psi_i = {(\frac{2\alpha_i}{\pi})^{3/4}}c_ie^{-\alpha_i r^{2}}$

## I) Deriving an Analytical Expression for 1e K.E. 

#

### 1. Setting Everything Up 

$\langle KE \rangle = \frac{-\hbar^{2}}{2m} \langle \psi_i |\nabla^{2}| \psi_j* \rangle$

$\nabla_r^{2} = \frac{1}{r^{2}}\frac{\partial}{\partial r} r^{2} \frac{\partial}{\partial r} + \frac{1}{r^{2}sin\theta} \frac{\partial}{\partial \theta}sin(\theta)\frac{\partial}{\partial \theta} + \frac{1}{r^{2}sin^{2}\theta}\frac{\partial^{2}}{\partial \theta^{2}}$


$\psi_i = c_ie^{-\alpha_ir^{2}}$, 
$\psi_j = c_je^{-\alpha_jr^{2}}$

Note $c_i = N*c_i$ as the wavefunctions have already been normalized!

### 2. Applying $\nabla_r^{2}$ to $\psi_i$ and $\psi_j$ 

Both wavefunctions are entirely radial (for now) so we all non-radial partial derivative terms of $\nabla_r^{2}$ should go to 0 when applied to $\psi_i$ or $\psi_j$ leaving us with 

$\langle KE \rangle = \frac{-\hbar^{2}}{2m} \langle \psi_i |\nabla^{2}| \psi_j* \rangle = \frac{-\hbar^{2}}{2m} \langle \psi_i | \frac{1}{r^{2}}\frac{\partial}{\partial r} r^{2} \frac{\partial}{\partial r} | \psi_j* \rangle = \frac{-\hbar^{2}}{2m} \langle c_ie^{-\alpha_ir^{2}} | \frac{1}{r^{2}}\frac{\partial}{\partial r} r^{2} \frac{\partial}{\partial r} | c_je^{-\alpha_jr^{2}} \rangle$

### 3. Let's Compute Some Derivatives 

$| \frac{1}{r^{2}}\frac{\partial}{\partial r} r^{2} \frac{\partial}{\partial r} | c_je^{-\alpha_jr^{2}} \rangle$

#### 3.1 : $\frac{\partial}{\partial r} c_je^{-\alpha_jr^{2}} = -2\alpha_j rc_je^{-\alpha_jr^{2}}$

#### 3.2 : $r^{2}\frac{\partial}{\partial r} c_je^{-\alpha_jr^{2}} = -2\alpha_j r^{3}c_je^{-\alpha_jr^{2}}$

#### 3.3 : $\frac{\partial}{\partial r}r^{2}\frac{\partial}{\partial r} c_je^{-\alpha_jr^{2}} = -6c_j\alpha_jr^{2}e^{-\alpha_jr^{2}} + 4c_j\alpha_j^{2}r^{4}e^{-\alpha_jr^{2}}$

#### 3.4 : $\frac{1}{r^{2}}\frac{\partial}{\partial r}r^{2}\frac{\partial}{\partial r} c_je^{-\alpha_jr^{2}} = -6c_j\alpha_je^{-\alpha_jr^{2}} + 4c_j\alpha_j^{2}r^{2}e^{-\alpha_jr^{2}}$

### 4. Plug our operated expression into the 3D Spherical Coordiantes Integral 

$\langle KE \rangle = \frac{-\hbar^{2}}{2m}\int_{0}^{2\pi}\,\int_{0}^{\pi}\,\int_{0}^{\infty} \psi_i|\nabla^{2}|\psi_j \, r^{2}sin(\theta)dr d\theta d\phi$

Using Same Thing from Previous Derivation This Can Simplify to: 

$\langle KE \rangle = \frac{-4\pi\hbar^{2}}{2m}\int_{0}^{\infty} \psi_i|\nabla^{2}|\psi_j \, r^{2}dr$

And then Applying $\nabla^{2}$ to $\psi_j$ 

$\langle KE \rangle = \frac{-4\pi\hbar^{2}}{2m}\int_{0}^{\infty} c_ie^{-\alpha_ir^{2}} (-6c_j\alpha_je^{-\alpha_jr^{2}} + 4c_j\alpha_j^{2}r^{2}e^{-\alpha_jr^{2}}) \, r^{2}dr$

$\langle KE \rangle = \frac{-4\pi\hbar^{2}}{2m} \int_{0}^{\infty} -6c_ic_j\alpha_jr^{2}e^{-(\alpha_i+\alpha_j)r^{2}} + 4c_ic_j\alpha_j^{2}r^{4}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr$

### 5. Getting a Free Integral: 

$\langle KE \rangle = \frac{-4\pi\hbar^{2}}{2m} (\int_{0}^{\infty} -6c_ic_jr^{2}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr  + \int_{0}^{\infty}4c_ic_j\alpha_j^{2}r^{4}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr)$

Notice that integral one is functionally the same as the overlap integral we evaluated earlier when doing overlap integrals so it should reduce to : 

$\int_{0}^{\infty} -6c_ic_j\alpha_j r^{2}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr = -6 * c_ic_j\alpha_j/2(\alpha_i+\alpha_j)^{3/2} * \sqrt{\pi}/2$

### 6. Getting A Less Free Integral: 

Now Considering Integral 2: $c_ic_j\int_{0}^{\infty}4c_ic_j\alpha_j^{2}r^{4}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr$

We can use the same integration by parts approach: 

$dv = re^{-(\alpha_i+\alpha_j)r^{2}}$ 

$v = \int \,dv = -e^{-(\alpha_i+\alpha_j)r^{2}}/2(\alpha_i+\alpha_j)$ 

$u = 4r^{3}$

$du = 12r^{2}$

$\int_{0}^{\infty} u \, dv = u*v - \int_{0}^{\infty} v du = \frac{-4r^{3}e^{-(\alpha_i+\alpha_j)r^{2}}}{2(\alpha_i+\alpha_j)} + \int_{0}^{\infty} \frac{12r^2e^{-(\alpha_i+\alpha_j)r^{2}}}{2(\alpha_i+\alpha_j)}\, dr$

Now We Rinse and Repeat for $\int_{0}^{\infty} \frac{12r^2e^{-(\alpha_i+\alpha_j)r^{2}}}{2(\alpha_i+\alpha_j)}\, dr$

$dv = re^{-(\alpha_i+\alpha_j)r^{2}}$ 

$v = \int \,dv = \frac{-e^{-(\alpha_i+\alpha_j)r^{2}}}{2(\alpha_i+\alpha_j)}$ 

$u = \frac{12r}{2(\alpha_i+\alpha_j)}$

$du = \frac{12}{2(\alpha_i+\alpha_j)}$

$\int_{0}^{\infty} u \, dv = u*v - \int_{0}^{\infty} v du = \frac{3re^{-(\alpha_i+\alpha_j)r^{2}}}{(\alpha_i+\alpha_j)^{2}} - \int_{0}^{\infty} \frac{-3e^{-(\alpha_i+\alpha_j)r^{2}}}{(\alpha_i+\alpha_j)^{2}}$

Adding up the two integral by parts 

$\int_{0}^{\infty}4c_ic_j\alpha_j^{2}r^{4}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr = c_ic_j\alpha_j^{2}(\frac{-2r^{3}e^{-(\alpha_i+\alpha_j)r^{2}}}{(\alpha_i+\alpha_j)} + \frac{-3re^{-(\alpha_i+\alpha_j)r^{2}}}{(\alpha_i+\alpha_j)^{2}} + \int_{0}^{\infty} \frac{3e^{-(\alpha_i+\alpha_j)r^{2}}}{(\alpha_i+\alpha_j)^{2}} \, dr)$

Using the Same Limit Idea for r -> 0 and r -> infinity the leading terms go to 0 so:

$\int_{0}^{\infty}4c_ic_j\alpha_j^{2}r^{4}e^{-(\alpha_i+\alpha_j)r^{2}} \, dr = \int_{0}^{\infty} \frac{+3\alpha_j^{2}e^{-(\alpha_i+\alpha_j)r^{2}}}{(\alpha_i+\alpha_j)^{2}} \, dr = +3\alpha_j^{2}c_ic_j/(\alpha_i+\alpha_j)^{5/2} * \sqrt{\pi}/2$

### 7 Plugging Back in and Simplify

$\langle KE \rangle = \frac{-4\pi\hbar^{2}}{2m}\int_{0}^{\infty} \psi_i|\nabla^{2}|\psi_j \, r^{2}dr$

$\langle KE \rangle =\frac{-2\pi^{3/2}\hbar^{2}}{m}(-3 * c_ic_j\alpha_j/2(\alpha_i+\alpha_j)^{3/2} + 3c_ic_j\alpha_j^{2}/2(\alpha_i+\alpha_j)^{5/2})$

$\langle KE \rangle =\frac{c_ic_j2\pi^{3/2}\hbar^{2}}{(\alpha_i+\alpha_j)^{5/2}m}(\frac{3\alpha_j}{2(\alpha_i+\alpha_j)^{-1}} - \frac{3\alpha_j^{2}}{2})$

$\langle KE \rangle =\frac{c_ic_j2\pi^{3/2}\hbar^{2}}{(\alpha_i+\alpha_j)^{5/2}m}(\frac{3\alpha_j^{2}+3\alpha_j\alpha_i}{2} - \frac{3\alpha_j^{2}}{2})$

$\langle KE \rangle =\frac{c_ic_j2\pi^{3/2}\hbar^{2}}{(\alpha_i+\alpha_j)^{5/2}m}(\frac{3\alpha_j\alpha_i}{2})$

$\langle KE \rangle =\frac{c_ic_j3\alpha_j\alpha_i\pi^{3/2}\hbar^{2}}{(\alpha_i+\alpha_j)^{5/2}m}$
### 8. Unit Conversion

In Hartree $E_h = \frac{\hbar}{m_ea_0^{2}}$
So we can plug in 1 for $\hbar$ and $m_e$ which yields:

$\langle KE \rangle = E_{ha} = \frac{3c_ic_j\alpha_j\alpha_i\pi^{3/2}}{(\alpha_i+\alpha_j)^{5/2}}$

## II) 1e Potential Energy Derivation

### 1. The 1e Potential Energy Operator Comes from Coulumbs Law:

$V(r) = k_ee/r^{2}$

Where r is nuclear distance, e is the electrical charge (1 in Ha) and k_e is Cuolomb constant =  $8.99*10^{9} N * m^{2}/C^{2}$

### 2. Using The Same Expectation Value Approach:

$\langle PE \rangle = \langle \psi_i | V(x) | \psi_j \rangle = \langle c_ie^{-\alpha_ir^{2}} | k_ee/r^{2} | c_je^{-\alpha_jr^{2}} \rangle$ 

### 3. Apply the Same Spherical Coordinate Simplification (assuming no angular nodes (l=0)) 

$\langle PE \rangle = 4\pi c_ic_jk_ee \int_{0}^{\infty} \frac{r^{2}e^{-(\alpha_i+\alpha_j)r^{2}}}{r^{2}} \, dr$

$\langle PE \rangle = 4\pi c_ic_jk_ee \int_{0}^{\infty} e^{-(\alpha_i+\alpha_j)r^{2}}\, dr$

### 4. Using the Integral Rules We've Already Established in Normalization (let e =1)

$\langle PE \rangle = \pi c_ic_jk_e (\frac {\pi}{(\alpha_i+\alpha_j)})^{1/2}$

### 5. Converting to Ha where k_e = 1 and e = 1
$\langle PE \rangle = \frac{\pi^{3/2} c_ic_j}{(\alpha_i+\alpha_j)^{1/2}}$

