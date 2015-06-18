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
fWeight =-1000;
fTurn = -500;

% Spoke properties
nSpokes = 32; %[-] Number of spokes
ESpokes = 210e9; %[Pa] Youngs module for spokes (steel)
aSpokes = pi*(2e-3)^2; %[mm^2] Area of spoke
t0 = 1e-3; %[N] Pretension in spokes

% Rim porperties for 3mm
Erim = 70e9;                %[N/m^2] Youngs modulus for the rim
EA = Erim * 186.74e-6;       %[N] Area of cross-section (from SW)
EIz = Erim * 4963.02e-12;   %[Nm^2] E-modul for twill (g�t) gange inertimoment fra SW
EIy = Erim * 31331.61e-12;  %[Nm^2] bending stiffness of rim around y axis

% % Rim porperties for 2mm
% Erim = 70e9;                %[N/m^2] Youngs modulus for the rim
% EA = Erim * 123.81e-6;       %[N] Area of cross-section (from SW)
% EIz = Erim * 3285.91e-12;   %[Nm^2] E-modul for twill (g�t) gange inertimoment fra SW
% EIy = Erim * 19958.06e-12;  %[Nm^2] bending stiffness of rim around y axis


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

% Forces are put in a vector, forces are applied in middle of rim
a = nSpokes*6;
f = zeros(a,1); f(a/2-4) = fy; f(a/2-3) = fz; %Definer kr�fter

% Stiffnessmatrix is generated
K = stiffnessmatrix(wheeldata);
% And solved for deflection
u = pinv(K)*f;

% Deformation is added to the spoke coordinates
spokesAfter = spokeCoordinates(wheeldata,u);

figure
subplot(1,2,1)
plotWheel(spokesAfter,1) % Wheel is plotted after deformation
subplot(1,2,2)
plotWheel(spokes,0) %Wheel are plotted before deformation

% Graphs of deformation are displayed
figure
plotU(wheeldata,K,fy,fz)
