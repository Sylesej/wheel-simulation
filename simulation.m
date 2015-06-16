%-------------------------------------------------------------------------
% SPOKED WHEEL SIMULATION
%
% Simulation of a spoked wheel as forces are applied in vertical and
% horizontal direction.
% The simulation is made by constructing a stiffnessmatrix using beam
% theory to calculate the stiffness-values, and then solving for
% displacements due to the applied forces.
% ------------------------------------------------------------------------

close all
clear;
clc;

%-------------------------------------------------------------------------
%     Variables
%-------------------------------------------------------------------------

% Load
fWeight =0;
fTurn = 250;

% Spoke properties
nSpokes = 32; %[-] Number of spokes
ESpokes = 210e9; %[Pa] Youngs module for spokes (steel)
aSpokes = pi*(2e-3)^2; %[mm^2] Area of spoke
t0 = 1e-3; %[N] Pretension in spokes

% Rim porperties
Erim = 70e9; % [N/m^2] Youngs modulus for the rim
EA = Erim * 119.15;%[N/m]
EIz = Erim * 2964e-8; %E-modul for twill (g�t) gange inertimoment fra SW
EIy = Erim * 0.40e-8; %[Nm^2] bending stiffness of rim around y axis

% Wheel dimensions
rHub = 10e-3; %[m] Radius of hub
rRim = 40e-2; %[m] Radius of rim (where spokes are attached)
wHub = 80e-3; %[m] Width of hub (where spokes are attached)


kSpokes = ESpokes*aSpokes; %[N/m] Stiffness of spokes
aBracing = asin(wHub/(rRim-rHub));%[Rad] Bracing angle of spokes
fy = fWeight;
fz = fTurn;

% Wheeldata struct
wheeldata.t0 = t0;
wheeldata.EIy = EIy;
wheeldata.EIz = EIz;
wheeldata.EA = EA;
wheeldata.nSpokes = nSpokes;
wheeldata.rHub = rHub;
wheeldata.rRim = rRim;
wheeldata.wHub = wHub;
wheeldata.ESpokes = ESpokes;
wheeldata.aSpokes = aSpokes;
wheeldata.kSpokes = kSpokes;

%-------------------------------------------------------------------------

u0 = zeros(nSpokes*6,1); % Vector holding deformation before load (0)

%Coordinates of spokes before load are generated
spokes = spokeCoordinates(wheeldata,u0);

figure %Wheel are plotted before deformation
plotWheel(spokes,0)

% Forces are put in a vector, forces are applied in middle of rim
a = nSpokes*6;
f = zeros(a,1); f(a/2-4) = fy; f(a/2-3) = fz; %Definer kr�fter

% Stiffnessmatrix is generated
K = stiffnessmatrix(wheeldata);
% And solved for deflection
u = pinv(K)*f;

% Deformation is added to the spoke coordinates
spokesAfter = spokeCoordinates(wheeldata,u);

% Wheel is plotted after deformation
figure
plotWheel(spokesAfter,1)

% Graphs of deformation are displayed
figure
plotU(wheeldata,K,fy,fz)
