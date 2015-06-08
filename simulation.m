clear
clc
close 'all'
%-------------------------------------------------------------------------
%     Initialization
%-------------------------------------------------------------------------
t0 = 1e-3; %[N] Pretension in spokes, assumed constant
nSpokes = 5; %[-] Number of spokes
rHub = 10e-3; %[m] Radius of hub
rRim = 40e-2; %[m] Radius of rim (where spokes are attached)
wHub = 80e-3; %[m] Width of hub (where spokes are attached)
ESpokes = 210e9; %[Pa] Youngs module for spokes (steel)
aSpokes = pi*(2e-3)^2; %[mm^2] Area of spoke
aBracing = asin(wHub/(rRim-rHub)); %Bracing angle of spokes are computed
EIx = 70e9 * 2964e-8;%E-modul for twill (gæt) gange inertimoment fra SW
%1560*10^(-3*4)*70*10^9; %[Nm^2] bending stiffness of rim about x axis (as happens from weight
EIz = 1; %[?] bending stiffness of rim about z axis (as happens from turn)
pointFactor = 0; %Faktor, der tilføjer til detaljering af FEM-beregningen.

%Ting til stivhedsmatricen
nPoints = pointFactor*nSpokes;


%Random
fWeight = 1000;

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

%Force from a turn is applied
%spokes = deformTurn(spokes,fTurn,wheeldata)