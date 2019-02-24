import torch
from math import sin, exp, sqrt, pi

def get_device():
    try:
        device = torch.device('cuda')
        assert device
    except:
        device = torch.device('cpu')
    return torch.device('cpu') #device

def zeros(m,n,k):
    return torch.zeros((m,n,k), device=get_device(), dtype=torch.double)

def diff(a,dim):
    if dim == 0:
        return a[1:,:,:] - a[:-1,:,:]
    elif dim == 1:
        return a[:,1:,:] - a[:,:-1,:]
    elif dim == 2:
        return a[:,:,1:] - a[:,:,:-1]
 
#physical constants
c    = 2.998e8
eta0 = 120*pi
mu0  = pi*4e-7
eps0 = 1e-9/(36*pi)
#environment parameters
width  = 1
height = 1
depth  = 1
f0     = 1e9
tw     = 1e-8/pi
t0     = 4*tw
#discretization parameters
dx = 0.01
dy = 0.01
dz = 0.01
nx = width  / dx
ny = height / dy
nz = depth  / dz
dt   = 0.95 / (c*sqrt(dx**-2+dy**-2+dz**-2))
#calculation parameters
n_iter = 100
#initalization
Hx = zeros(nx,ny,nz)
Hy = zeros(nx,ny,nz)
Hz = zeros(nx,ny,nz)
Ex = zeros(nx,ny,nz)
Ey = zeros(nx,ny,nz)
Ez = zeros(nx,ny,nz)
#iteration
i = 0

for n in range(n_iter):

    #derivatives
   
    Hxy = diff(Hx,1);
    Hxz = diff(Hx,2);
    Hzx = diff(Hz,0);
    Hzy = diff(Hz,1);
    Hyx = diff(Hy,0);
    Hyz = diff(Hy,2);

    #Maxwell Equations  

    Ex[:,1:-1,1:-1] += (dt/(eps0*dy))*Hzy[:,:-1,1:-1] - (dt/(eps0*dz))*Hyz[:,1:-1,:-1]
    Ey[1:-1,:,1:-1] += (dt/(eps0*dz))*Hxz[1:-1,:,:-1] - (dt/(eps0*dx))*Hzx[:-1,:,1:-1]
    Ez[1:-1,1:-1,:] += (dt/(eps0*dx))*Hyx[:-1,1:-1,:] - (dt/(eps0*dy))*Hxy[1:-1,:-1,:]

    #Gaussian Source

    Ez[int(nx/2),int(ny/2),int(nz/2)] += sin(2*pi*f0*n*dt)*exp(-(n*dt-t0)**2/(tw**2))/dy

    #derivatives

    Exy = diff(Ex,1);
    Exz = diff(Ex,2);############################!!!!!!!!!!!
    Ezx = diff(Ez,0);
    Ezy = diff(Ez,1);
    Eyx = diff(Ey,0);
    Eyz = diff(Ey,2);
    
    #Maxwell Equations 

    Hx[:,:-1,:-1] += -(dt/(mu0*dy))*Ezy[:,:,:-1] + (dt/(mu0*dz))*Eyz[:,:-1,:]
    Hy[:-1,:,:-1] += -(dt/(mu0*dz))*Exz[:-1,:,:] + (dt/(mu0*dx))*Ezx[:,:,:-1]
    Hz[:-1,:-1,:] += -(dt/(mu0*dx))*Eyx[:,:-1,:] + (dt/(mu0*dy))*Exy[:-1,:,:]

    ####################### I AM HERE #########################

    """
    #display
    if not i%5:
       pcolor(Ez[:,int(ny/2),:])
    drawnow
    """

    i += 1
    print(i)


