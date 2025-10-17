from Basis import Basis #importing the Basis Class
k=Basis('2pz', number_of_gaussians=3, p_thresh=0.01, plot_type='Radial', r_inc=0.05, rmax=6) #generate a basis object with desired params for 2pz orbital
k.kernel() #run basis calculations 
