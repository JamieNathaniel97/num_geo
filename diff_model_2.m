%***** 2D DIFFUSION MODEL OF HEAT TRANSPORT *******************

%***** Initialise Model Setup

% Image data obtained

% create x-coordinate vectors
xc = h/2:h:W-h/2;               % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;               % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;                     % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;                     % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc);      % create 2D coordinate arrays

%Coordinate grid on the output / for code interpretation: use this as the
%model grid
%h = W/Nx;
%xc = h/2:h:W-h/2;
%zc = h/2:h:D-h/2;
%[Xc,Zc] = meshgrid(xc,zc);

%Insulating Sides
ix3 = [       1,1:Nx,Nx      ];
ix5 = [    1, 1,1:Nx,Nx,Nx   ];
ix7 = [ 1, 1, 1,1:Nx,Nx,Nx,Nx];
%Set up Insulating Top and Bottom
iz3 = [       1,1:Nz,Nz      ];
iz5 = [    1, 1,1:Nz,Nz,Nz   ];
iz7 = [ 1, 1, 1,1:Nz,Nz,Nz,Nz];

% create smooth random perturbation field
rng(15);
dr = randn(Nz,Nx);
for ii = 1:10
    dr = dr + (diff(dr(iz3,:),2,1) + diff(dr(:,ix3),2,2))/8;
end

% set initial condition for temperature at cell centres
T = Ttop + geotherm.*Zc + dr*0; % initialise T array on linear gradient
T(air) = Ttop;

% initialise density and mobility
rho = rho0.*(1 - aT.*(T-Ttop));
kT = kappa.*ones(Nz,Nx);

% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,0,yr)

%***** Solve Model Equations
dt = CFL * (h/2)^2/max(kT(:)); % initial time step [s]
t = 0; % initial time [s]
k = 0; % initial time step count

% loop through time steps until stopping time reached
while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

    % print time step header
    %fprintf(1,'\n\n***** step = %d; dt = %1.3e; time = %1.3e \n\n',k,dt,t)

    % store old temperature and rate
    dTdt = dTdto;
    To = T;

    dTdt1 = diffusion(T           ,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);
    dTdt2 = diffusion(T+dTdt1/2*dt,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);
    dTdt3 = diffusion(T+dTdt2/2*dt,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);
    dTdt4 = diffusion(T+dTdt3  *dt,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);

    T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt;

    Hs = Hr ./ (rho.*Cp);

    T(air) = Ttop;
    % plot model progress every 'nop' time steps
    
    if ~mod(k,nop)
        makefig(xc,zc,T,t,yr);
    end
    
end

%***** Utility Functions ************************************************

% Function to make output figure
function makefig(x,z,T,t,yr)

% plot temperature in subplot 1
imagesc(x,z,T); axis equal; c = colorbar; hold on
contour(x,z,T,[100,150,200],'k');
drawnow

[C, h] = contour(x, z, T, [150,150], 'r', 'Linewidth', 2; % adds contour line at 150m
ylabel(c,'[°C]','FontSize',15)
ylabel('Depth [m]','FontSize',15)
xlabel('Horizontal Distance [m]','FontSize',15)
title(['Temperature [C]; time = ',num2str(t/yr), 'years'],'FontSize',17)

end

% Function to calculate diffusion rate
function [dTdt] = diffusion(f,k,h,ix,iz,geotherm,Hr,rho,Cp)

% calculate heat flux by diffusion
kx = (k(:,ix(1:end-1)) + k(:,ix(2:end)))/2;
kz = (k(iz(1:end-1),:) + k(iz(2:end),:))/2;
qx = - kx .* diff(f(:,ix), 1, 2)/h;
qz = - kz .* diff(f(iz,:), 1, 1)/h;

% basal boundary
qz(end,:) = - kz(end,:) .* geotherm;

% calculate flux balance for rate of change
dTdt_diffuion = - (diff(qx,1,2)/h+diff(qz,1,1)/h);

% add heat source from variable table
heat_source = Hr./(rho .* Cp);

% add to diffusion
dTdt = dTdt_diffuion + heat_source;

end