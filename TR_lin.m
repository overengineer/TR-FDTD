%physical constants
clear all;
close all;
load 'withtumor.mat';
load 'withouttumor.mat';

c0    = 2.998e8;
eta0 = 120*pi;
mu0  = pi*4e-7;
eps0 = 1e-9/(36*pi);
%box dimensions
width  = 0.5; % 30cm 
height = 0.5;
length  = 0.5; % 1cm
%source parameters
f0     = 1e9; % GHz
band   = 2e9;
tw     = sqrt(-log(0.1)/(pi*band)^2);%1e-8/pi;
t0     = 4*tw;
%spatial discretization
adipose = 5;
tumor   = 10;
sigma   = 0;
epsr    = tumor;
w    = 2 * pi * band;
k    = (w/c0)*sqrt(epsr-1j*sigma/(w*eps0));
beta = real(k);
c    = w / beta;
lambda = c/f0;
dxmax  = lambda / 20;
dx = dxmax;
dy = dx;
dz = dx;
nx = round(width/dx);
ny = round(height/dy);
nz = round(length/dz);

%source position
srcx = round(nx / 2);
srcy = round( 3 * ny / 4);
srcz = round(nz / 2);

% material
eps = ones(nx,ny,nz) * eps0; %* adipose;
sigma = zeros(nx,ny,nz);% * f0 * 1e-9 * 0.5 - 0.5;
%temporal discretization
dt   = 0.99/(c0*sqrt(dx^-2+dy^-2+dz^-2));

rec1 = trec - rec;
tau = 20e-12;
[foo,tp] = max(abs(rec1),[],2); 
for k=1:1:nrec
    recn(k,:) = exp(-((dt*((1:1:n_iter)-tp(k)))/tau).^2) .* rec1(k,:);
end
% hold on
% plot(rec(15,:))
% plot(exp(-((dt*((1:1:n_iter)-tp(15)))/tau).^2))
% draw now
% 
% while 1
% end


%EM field dimensions
Hx = zeros(nx,ny,nz);
Hy = zeros(nx,ny,nz);
Hz = zeros(nx,ny,nz);
Ex = zeros(nx,ny,nz);
Ey = zeros(nx,ny,nz);
Ez = zeros(nx,ny,nz);
%iteration
i = 0;
for n=1:1:n_iter
    %magnetic field derivatives
    Hxy = diff(Hx,1,2);
    Hxz = diff(Hx,1,3);
    Hzx = diff(Hz,1,1);
    Hzy = diff(Hz,1,2);
    Hyx = diff(Hy,1,1);
    Hyz = diff(Hy,1,3);
    %electric field maxwell equations
    epsi = eps(:,2:end-1,2:nz-1);
    ksi = (dt * sigma(:,2:end-1,2:nz-1)) ./ ( 2 * epsi );
    c2 = (1./(1+ksi)).*(dt./epsi);
    c1 = (1-ksi)./(1+ksi);
    Ex(:,2:end-1,2:end-1) = c1.*Ex(:,2:end-1,2:nz-1) - c2.*((1/dy)*Hzy(:,1:end-1,2:end-1) - (1/dz)*Hyz(:,2:ny-1,1:end-1));
    
     
    epsi = eps(2:end-1,:,2:end-1);
    ksi = (dt * sigma(2:end-1,:,2:end-1)) ./ ( 2 * epsi );
    c2 = (1./(1+ksi)).*(dt./epsi);
    c1 = (1-ksi)./(1+ksi);
    Ey(2:end-1,:,2:end-1) = c1.*Ey(2:end-1,:,2:end-1) - c2.*((1/dz)*Hxz(2:end-1,:,1:end-1) - (1/dx)*Hzx(1:end-1,:,2:end-1));
    
    epsi = eps(2:end-1,2:end-1,:);
    ksi = (dt * sigma(2:end-1,2:end-1,:)) ./ ( 2 * epsi );
    c2 = (1./(1+ksi)).*(dt./epsi);
    c1 = (1-ksi)./(1+ksi);
    Ez(2:end-1,2:end-1,:) = c1.*Ez(2:end-1,2:end-1,:) - c2.*((1/dx)*Hyx(1:end-1,2:end-1,:) - (1/dy)*Hxy(2:end-1,1:end-1,:));
   
    %TR sources
    for k=1:nrec
        Ez(recdx * k, recy, recz) = Ez(recdx * k, recy, recz) + recn(k, n_iter-n+1);
    end
    %Ez(recx, recdy , recz)
    %rec(1,n_iter-n)
    %electric field derivatives
    Exy = diff(Ex,1,2);
    Exz = diff(Ex,1,3);
    Ezx = diff(Ez,1,1);
    Ezy = diff(Ez,1,2);
    Eyx = diff(Ey,1,1);
    Eyz = diff(Ey,1,3);
    
    %magnetic field maxwell equations
    Hx(:,1:end-1,1:end-1) = Hx(:,1:end-1,1:end-1) + (dt/(mu0*dy))*Ezy(:,:,1:end-1) - (dt/(mu0*dz))*Eyz(:,1:end-1,:);
    Hy(1:end-1,:,1:end-1) = Hy(1:end-1,:,1:end-1) + (dt/(mu0*dz))*Exz(1:end-1,:,:) - (dt/(mu0*dx))*Ezx(:,:,1:end-1);
    Hz(1:end-1,1:end-1,:) = Hz(1:end-1,1:end-1,:) + (dt/(mu0*dx))*Eyx(:,1:end-1,:) - (dt/(mu0*dy))*Exy(1:end-1,:,:);

    %display
    if 1 %n>120 && n<160)
        %slice(:,:)=Ez(30:60,round(ny/2)-20:round(ny/2)+3,srcz);
        slice(:,:)=Ez(35:55,35:55,srcz);
        pcolor(slice.');
        colorbar;
        shading interp
        drawnow
    end
    i = i+1;
    disp(i)
     
    %R(n) = varimax_norm(Ez(30:56,round(ny/2)-20:round(ny/2)+3,srcz));
    R(n) = varimax_norm(Ez(35:55,35:55,srcz));
end

figure;plot(R)

function R = varimax_norm(Ez)
    R = sum(sum(sum(Ez.^2)))^2 / sum(sum(sum(Ez.^4)));
end
