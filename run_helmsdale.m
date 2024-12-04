%***** RUN 1D ADVECTION DIFFUSION MODEL ***********************************

% clear workspace
clear all; close all; %clc;

% set model parameters

run('./run_from_img.m');

TGrad = 35;            % temperature gradient
Ttop  = 5;            % surface temperature
Tbot  = Ttop+TGrad*(D/1000); % top/base T-gradient

KD0   = 1e-8;         % Darcy mobility coefficient [m2/Pas]
cp    = 0e-14;        % KD p-dependence prefactor
mp    = 0;            % KD p-dependence powerlaw
aT    = 1e-4;         % thermal expansivity [1/C]
rho0  = 1000;         % reference density [kg/m3]
kT0   = 1e-7;         % heat diffusivity [m2/s]
cT    = 1e-9;         % kT T-dependence prefactor
mT    = 2;            % kT T-dependence powerlaw
g0    = 9.8;           % gravity [m/s2]

ADVN  = 'WENO5';      % advection scheme ('UPW1', 'CFD2', 'UPW3', 'WENO5')

yr    = 3600*24*365;  % seconds per year [s]
tend  = 1e4*yr;       % stopping time [s]
CFL   = 1/2;          % Time step limiter
nop   = 10;           % output figure produced every 'nop' steps
alpha = 0.99;         % iterative step size limiter
beta  = 0.95;         % iterative lag parameter
tol   = 1e-8;         % residual tolerance
nup   = 100;          % update T, check residual every nup iterations

%*****  RUN MODEL

run('./diff_model_helmsdale.m');