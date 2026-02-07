import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import os

output_dir = "saved_images"
os.makedirs(output_dir, exist_ok=True)

theta = np.loadtxt('theta.out')
Nx = 128
Ny = 128
theta_matrix = np.zeros((Nx,Ny))
for k in range (0,201,1):
    ii=0
    theta_plot = theta[k*(Nx*Ny):(k+1)*(Nx*Ny)]
    for i in range (0,Nx):
        for j in range (0,Ny):
            theta_matrix[i][j]=theta_plot[ii]
            ii=ii+1
    
    plt.figure(figsize=(5,4))
    contour = plt.contourf(theta_matrix.T,cmap='plasma',levels=200)
    contour.set_clim(-0.6,0.15)
    cbar = plt.colorbar(contour)
    plt.title(r'$\theta$',size=24)
    plt.xlabel(r'$x$',size=22)
    plt.ylabel(r'$y$',size=22)
    
    plt.gca().xaxis.set_major_locator(MultipleLocator(16))
    plt.gca().xaxis.set_major_formatter('{x:.0f}')
    plt.gca().xaxis.set_minor_locator(MultipleLocator(8))
    plt.gca().yaxis.set_major_locator(MultipleLocator(16))
    plt.gca().yaxis.set_major_formatter('{x:.0f}')
    plt.gca().yaxis.set_minor_locator(MultipleLocator(8))
    
    plt.gca().tick_params(axis="both", which="both", top=True, right=True, bottom=True, left=True)
    plt.gca().tick_params(axis="x", which="both", direction="out", bottom=True)
    plt.gca().tick_params(axis="x", which="both", direction="in", top=True)
    plt.gca().tick_params(axis="y", which="both", direction="out", left=True)
    plt.gca().tick_params(axis="y", which="both", direction="in", right=True)
    
    plt.savefig(os.path.join(output_dir,f"poly_{k+1}.jpg"),dpi=300,bbox_inches='tight',transparent=True)
    plt.close()