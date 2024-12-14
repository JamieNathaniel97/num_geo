%***** 2D DIFFUSION MODEL OF HEAT TRANSPORT *******************

%***** Initialise Model Setup

% Image data obtained

% create x-coordinate vectors
xc = h/2:h:W-h/2;               % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;               % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;                     % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;                     % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc);      % create 2D coordinate arrays

%Insulating Sides
ix3 = [       1,1:Nx,Nx      ];
%Set up Insulating Top and Bottom
iz3 = [       1,1:Nz,Nz      ];

% set initial condition for temperature at cell centres
T = Ttop + geotherm.*Zc; % initialise T array on linear gradient
T(air) = Ttop;

% initialise density and mobility
rho = rho0.*(1 - aT.*(T-Ttop));
%kT = kappa.*ones(Nz,Nx);
kT = kappa;

% initialise output figure with initial condition
figure(1); clf
makefig(xc,zc,T,0,yr)

%***** Solve Model Equations
dt = CFL * (h/2)^2/max(kT(:)); % initial time step [s]
t = 0; % initial time [s]
k = 0; % initial time step count

if gauss_test % (ChatGPT helped with this whole part)
    % initialise gaussian T for error test
    T = T0 + T_p .* exp( -( ( (Xc - W/2).^2 + (Zc - D/2).^2 ) / (4*T_w^2) ) );
    k0 = kT(1,1);  
end

% loop through time steps until stopping time reached
while t <= tend

    % increment time and step count
    t = t+dt;
    k = k+1;

% RK4 scheme steps

    dTdt1 = diffusion(T           ,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);  %R1
    dTdt2 = diffusion(T+dTdt1*dt/2,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);  %R2
    dTdt3 = diffusion(T+dTdt2*dt/2,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);  %R3
    dTdt4 = diffusion(T+dTdt3*dt  ,kT,h,ix3,iz3,geotherm,Hr,rho,Cp);  %R4

    T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)*dt/6 ;

    Hs = Hr ./ (rho.*Cp);

    T(air) = Ttop;
    % plot model progress every 'nop' time steps
    
    if ~mod(k,nop)
        makefig(xc,zc,T,t,yr);
    end
    
end

if gauss_test

    % analytical solution
    
    % Effective Gaussian width 
    w_effective = sqrt(T_w^2 + 4 * k0 * t); 
    
    % Analytical solution
    Ta = T0 + (T_p / (w_effective.^2)) * ...
           exp(-((Xc - W/2).^2 + (Zc - D/2).^2) / w_effective.^2);

    % calculate error of final T
    Errorz = norm(T-Ta,1)./norm(Ta,1);
    Errorx = norm(T-Ta,2)./norm(Ta,2);
    
    
    disp(['Numerical error on x = ', num2str(Errorx)]);
    disp(['Numerical error on z = ', num2str(Errorz)]);
    

end

%***** Utility Functions ************************************************




% Function to make output figure
function makefig(x,z,T,t,yr)

% plot temperature in subplot 1
imagesc(x,z,T); axis equal; c = colorbar; hold on
[C, h] = contour(x,z,T,[50,100,150],'k');

[j, g] = contour(x, z, T, [150,150], 'r', 'Linewidth', 2); % adds contour line at 150m (ChatGPT helped here)
clabel(C, h, 'Fontsize',12,'Color', 'r')
clabel(j, g, 'Fontsize',12,'Color', 'r')
drawnow
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
