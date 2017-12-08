# coding: utf-8

# Modules à utiliser
import scipy.sparse as sp
import numpy as np

def matrix_neumann2D(Omega,Nx,Ny):
    """Returns the matrix that discretizes the mapping u --> -Laplacian(u)
on the domain Omega = [xmin,xmax,ymin,ymax] split into Nx parts in x
and Ny parts in y. The final matrix is a scipy.sparse CSR matrix. It
is of size (Nx+1)*(Ny+1)."""
    
    hx = (Omega[1]-Omega[0])/Nx
    hy = (Omega[3]-Omega[2])/Ny
    hx2 = hx*hx
    hy2 = hy*hy

    # Les inconnues sont numérotés de 0 à Nx suivant x et 0 à Ny
    # suivant y. La taille du problème est donc (Nx+1)*(Ny+1).

    # Pour -Laplacien(u), la matrice est constituée de (Ny+1)x(Ny+1)
    # blocs de taille (Nx+1)x(Nx+1), de la forme
    #
    # A = [ A0 B              ]
    #     [ B  A1 B           ]
    #     [    B  A1 B        ]
    #     [       .  .  .     ]
    #     [           B A1 B  ]
    #     [             B  A0 ]
    #
    # Au final, on peut commencer à remplir avec des diagonales
    N = (1+Nx)*(1+Ny)
    diags = np.zeros((5,N))
    # La diagonale est constante
    diags[2,:] = 2./hx2+2./hy2
    # Diagonale -1
    diags[1,:] = -1./hx2                      # en général
    diags[1,np.arange(Nx,N,Nx+1)] = 0.        # bord gauche
    diags[1,np.arange(Nx-1,N,Nx+1)] = -2./hx2 # bord droit
    # Diagonale +1
    diags[3,:] = -1./hx2                   # en général
    diags[3,np.arange(0,N,Nx+1)] = 0.      # bord droit
    diags[3,np.arange(1,N,Nx+1)] = -2./hx2 # bord gauche
    # Diagonale -(Nx+1)
    diags[0,:] = -1./hy2                       # en général
    diags[0,(Nx+1)*(Ny-1):(Nx+1)*Ny] = -2./hy2 # bord bas
    # Diagonale +(Nx+1)
    diags[4,:] = -1./hy2             # en général
    diags[4,Nx+1:2*(Nx+1)] = -2./hy2 # bord haut

    # Construction de la matrice creuse de u --> -Laplacien(u)
    A = sp.spdiags(diags,[-(Nx+1),-1,0,1,(Nx+1)], (Nx+1)*(Ny+1),
                   (Nx+1)*(Ny+1), format="csc")

    return A
