clear
clc
close 'all'
%-------------------------------------------------------------------------
%     Initialization
%-------------------------------------------------------------------------
t0 = 1e-3; %[N] Pretension in spokes, assumed constant
nSpokes = 32; %[-] Number of spokes
rHub = 10e-3; %[m] Radius of hub
rRim = 40e-2; %[m] Radius of rim (where spokes are attached)
wHub = 80e-3; %[m] Width of hub (where spokes are attached)
ESpokes = 210e9; %[Pa] Youngs module for spokes (steel)
aSpokes = pi*(2e-3)^2; %[mm^2] Area of spoke
aBracing = asin(wHub/(rRim-rHub)); %Bracing angle of spokes are computed
EIz = 70e6 * 2928.464;%E-modul for twill (gæt) gange inertimoment fra SW
%1560*10^(-3*4)*70*10^9; %[Nm^2] bending stiffness of rim about x axis (as happens from weight
EIy = 70e9 * 17924.40; %[?] bending stiffness of rim about z axis (as happens from tu
EA = 70e9 * 119.15; %For the stiffnessM.
pointFactor = 0; %Faktor, der tilføjer til detaljering af FEM-beregningen.
kSpokes = ESpokes*aSpokes;
%Ting til stivhedsmatricen
nPoints = pointFactor*nSpokes;


%Random
fWeight =-1000;
fTurn = 1000;

fy = fWeight;
fz = fTurn;

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
wheeldata.nPoints = nPoints;
wheeldata.pointFactor = pointFactor;
wheeldata.kSpokes = kSpokes;

%-------------------------------------------------------------------------

%Spoke tension is initialized before external load is applied.
tSpoke = t0.*ones(1,nSpokes); 

u0 = zeros(nSpokes*6,1);

%Coordinates of spokes are generated
spokes = spokeCoordinates(wheeldata,u0);

% Forces are put in a vector, forces applied in middle of rim
a = nSpokes*6;
f = zeros(a,1); f(a/2-4) = fy; f(a/2-3) = fz; %Definer kræfter

% Stiffnessmatrix is generated
K = stiffnessmatrix(wheeldata);
% And solved for deflection
u = pinv(K)*f;


%sMh = basisskiftesMh(K);
%hMs = inv(sMh);
%Ks = sMh*K*hMs;

%Wheel are plotted before deformation
figure
plotWheel(spokes,0)

spokesAfter = spokeCoordinates(wheeldata,u);

figure
plotWheel(spokesAfter,1)

figure
plotU(wheeldata,K,fy,fz)