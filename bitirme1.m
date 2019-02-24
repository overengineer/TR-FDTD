%TMz Polarization
%physical constants
c    = 2.998e8;
eta0 = 120*pi;
mu0  = pi*4e-7;
eps0 = 1e-9/(36*pi);
%environment parameters
width  = 1;
height = 1;
tw     = 1e-9/pi;%?
t0     = 4*tw;
%discretization parameters
dx = 0.005;
dy = 0.005;
nx = width/dx;
ny = height/dy;
dt   = 0.95/(c*sqrt(dx^-2+dy^-2));
%calculation parameters
n_iter = 1000;
c1 = dt/(mu0*dy);
c2 = dt/(mu0*dx);
c3 = dt/(eps0*dx);
c4 = dt/(eps0*dy);
%initalization
Hx = zeros(nx,ny);
Hy = zeros(nx,ny);
Ez = zeros(nx,ny);
%iteration
for n=1:1:n_iter
    %Maxwell Equations (TMz)
    Ezx = diff(Ez,1,1);
    Ezy = diff(Ez,1,2);
    Hx(2:nx-1,2:ny) = Hx(2:nx-1,2:ny) - c1*Ezy(2:nx-1,:);
    Hy(2:nx,2:ny-1) = Hy(2:nx,2:ny-1) + c2*Ezx(:,2:ny-1);
    Hxy = diff(Hx,1,2);
    Hyx = diff(Hy,1,1);
    Ez(2:nx-1,2:ny-1) = Ez(2:nx-1,2:ny-1) + c3*Hyx(2:nx-1,2:ny-1) - c4*Hxy(2:nx-1,2:ny-1);
    %Gaussian Source
    f(n)= exp(-(n*dt-t0)^2/(tw^2))/dy;
    Ez(round(nx/2),round(ny/2)) = Ez(round(nx/2),round(ny/2)) + f(n);
    %Neuman Condition
    Ez(:,2)  = -Ez(:,1);
    Ez(2,:)  = -Ez(1,:);
    Ez(:,ny-1) = -Ez(:,ny);
    Ez(nx-1,:) = -Ez(nx,:);
    %display
    %n = n + 1;
    pcolor(Ez);
    shading interp;
    drawnow
end

