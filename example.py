from Basis import Basis #importing the Basis Class

#Calculate Overlap Integral B/W 2 different wave functions
k = Basis('2pz', number_of_gaussians=3, p_thresh=0.01, plot_type='False', r_inc=0.05, rmax=6)
k.kernel()
j=Basis('2s', number_of_gaussians=3, p_thresh=0.01, plot_type='False', r_inc=0.05, rmax=6)
j.kernel()
print('STO Overlap: ', j.overlap(j.r, j.psi, k.psi)) #STO Overlap 2s/2p
print('STO-3g Overlap: ', j.overlap(j.r, j.p, k.p)) #STO-nG Overlap 2s/2p

k=Basis('2pz', number_of_gaussians=3, p_thresh=0.01, plot_type='Radial', r_inc=0.05, rmax=6) #generate a basis object with desired params for 2pz orbital
k.kernel() #run basis calculations 
