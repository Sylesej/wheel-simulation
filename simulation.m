clear;
clc;
clf;
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
EIx = 70e9 * 2964e-8;%E-modul for twill (g�t) gange inertimoment fra SW
%1560*10^(-3*4)*70*10^9; %[Nm^2] 
%bending stiffness of rim about x axis (as happens from weight
EIz = 1; %[?] bending stiffness of rim about z axis (as happens from turn)
pointFactor = 0; %Faktor, der tilf�jer til detaljering af FEM-beregningen.

%Ting til stivhedsmatricen
nPoints = pointFactor*nSpokes;

wheeldata.t0 = t0;
wheeldata.EIx = EIx;
wheeldata.EIz = EIz;
wheeldata.nSpokes = nSpokes;
wheeldata.rHub = rHub;
wheeldata.rRim = rRim;
wheeldata.wHub = wHub;
wheeldata.ESpokes = ESpokes;
wheeldata.aSpokes = aSpokes;
wheeldata.nPoints = nPoints;
wheeldata.pointFactor = pointFactor;


%-------------------------------------------------------------------------

%Spoke tension is initialized before external load is applied.
tSpoke = t0.*ones(1,nSpokes); 

%Coordinates of spokes are generated
spokes = spokeCoordinates(wheeldata);
K = stiffnessmatrix(wheeldata);

%Wheel are plotted before deformation
%plotWheel(spokes)

%
f = zeros(nSpokes*6,1); 
f(nSpokes*6/2+2) = 1000; f(nSpokes*6/2+3) = 250; 
%Indfør kræfter i hhv y og z.
u = pinv(K)*f;

subplot(2,2,1)
plot(u(1:6:end))
title('Udbøjninger i x')

subplot(2,2,2)
plot(u(2:6:end))
title('Udbøjninger i y')

subplot(2,2,3)
plot(u(3:6:end))
title('Udbøjninger i z')

fEg = t0 + ESpokes*aSpokes*u(2:6:end);

subplot(2,2,4)
plot(fEg)
title('Egekræfter')