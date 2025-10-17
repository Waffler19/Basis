from Basis import Basis
k=Basis('2pz', number_of_gaussians=3, p_thresh=0.01, plot_type='Radial', r_inc=0.05, rmax=6)
k.kernel()
