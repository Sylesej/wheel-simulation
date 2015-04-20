clear
clc
close 'all'
%-------------------------------------------------------------------------
%     Initialization
%-------------------------------------------------------------------------
t0 = 100; %[N] Pretension in spokes, assumed constant
nSpokes = 32; %[-] Number of spokes
rHub = 10e-3; %[m] Radius of hub
rRim = 40e-2; %[m] Radius of rim (where spokes are attached)
wHub = 80e-3; %[m] Width of hub (where spokes are attached)
ESpokes = 210e9; %[Pa] Youngs module for spokes (steel)
aSpokes = pi*(2e-3)^2; %[mm^2] Area of spoke
aBracing = asin(wHub/(rRim-rHub)); %Bracing angle of spokes are computed

wheeldata.nSpokes = nSpokes;
wheeldata.rHub = rHub;
wheeldata.rRim = rRim;
wheeldata.wHub = wHub;
%-------------------------------------------------------------------------

%Spoke tension is initialized before external load is applied.
tSpoke = t0.*ones(1,nSpokes); 

%Coordinates of spokes are generated
spokes = spokeCoordinates(wheeldata);

%Wheel are plotted before deformation
plotWheel(spokes)

%Concentrated force is applied
%spokes = deformConc(spokes,
