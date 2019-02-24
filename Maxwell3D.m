%physical constants
clear all;
c    = 2.998e8;
eta0 = 120*pi;
mu0  = pi*4e-7;
eps0 = 1e-9/(36*pi);
%environment parameters
width  = 1;
height = 1;
depth  = 1;
f0     =1e9;
tw     = 1e-8/pi;%?
t0     = 4*tw;
%discretization parameters
dx = 0.01;
dy = dx;
dz = dx;
nx = width/dx;
ny = height/dy;
nz = depth/dz;
dt   = 0.95/(c*sqrt(dx^-2+dy^-2+dz^-2));
%calculation parameters
n_iter = 2000;
%initalization
Hx = zeros(nx,ny,nz);
Hy = zeros(nx,ny,nz);
Hz = zeros(nx,ny,nz);
Ex = zeros(nx,ny,nz);
Ey = zeros(nx,ny,nz);
Ez = zeros(nx,ny,nz);
%iteration
i = 0
for n=1:1:n_iter
    %derivatives
   
    %H
    Hxy = diff(Hx,1,2);
    Hxz = diff(Hx,1,3);
    Hzx = diff(Hz,1,1);
    Hzy = diff(Hz,1,2);
    Hyx = diff(Hy,1,1);
    Hyz = diff(Hy,1,3);
    %Maxwell Equations
    % boyutlar e? uzunlukta olmas? için türev al?nmayan boyutu bir eksik
    % indeksliyoruz. E vektörlerinin türev boyutunu k?saltmadan al?p di?rev
    % vektörün türev boyutunda bir ileriden ba?lat?yoruz.
    
    %E
    Ex(:,2:end-1,2:end-1) = Ex(:,2:end-1,2:end-1) + (dt/(eps0*dy))*Hzy(:,1:end-1,2:end-1) - (dt/(eps0*dz))*Hyz(:,2:ny-1,1:end-1); %OK
    Ey(2:end-1,:,2:end-1) = Ey(2:end-1,:,2:end-1) + (dt/(eps0*dz))*Hxz(2:end-1,:,1:end-1) - (dt/(eps0*dx))*Hzx(1:end-1,:,2:end-1); %OK
    Ez(2:end-1,2:end-1,:) = Ez(2:end-1,2:end-1,:) + (dt/(eps0*dx))*Hyx(1:end-1,2:end-1,:) - (dt/(eps0*dy))*Hxy(2:end-1,1:end-1,:); %OK
     %Gaussian Source
      f(n)= sin(2*pi*f0*n*dt)*exp(-(n*dt-t0)^2/(tw^2))/dy;
    Ez(round(nx/2),round(ny/2),round(nz/2)) = Ez(round(nx/2),round(ny/2),round(nz/2)) + f(n);
    Ezn(n)=Ez(round(nx/2),round(ny/2)-4,round(nz/2));
    Exy = diff(Ex,1,2);
    Exz = diff(Ex,1,3);
    Ezx = diff(Ez,1,1);
    Ezy = diff(Ez,1,2);
    Eyx = diff(Ey,1,1);
    Eyz = diff(Ey,1,3);
    
    %H
    Hx(:,1:end-1,1:end-1) = Hx(:,1:end-1,1:end-1) - (dt/(mu0*dy))*Ezy(:,:,1:end-1) + (dt/(mu0*dz))*Eyz(:,1:end-1,:); %OK
    Hy(1:end-1,:,1:end-1) = Hy(1:end-1,:,1:end-1) - (dt/(mu0*dz))*Exz(1:end-1,:,:) + (dt/(mu0*dx))*Ezx(:,:,1:end-1); %OK
    Hz(1:end-1,1:end-1,:) = Hz(1:end-1,1:end-1,:) - (dt/(mu0*dx))*Eyx(:,1:end-1,:) + (dt/(mu0*dy))*Exy(1:end-1,:,:); %OK
    %display
    if (mod(i,5)==0)
        aa(:,:)= Ez(:,round(ny/2),:);
        a = log(aa.*(aa>0));
        a(1,1)= -60;
        a(1,2)=0;
        a(end,1)=-60;
        a(end,end)= -60;
        a(1,end)=-60;
       pcolor(a);colorbar
    end
    %shading interp;
    drawnow
    i = i+1
end

