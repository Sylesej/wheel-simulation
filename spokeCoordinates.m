function spokes=spokeCoordinates(wheeldata)
%This function defines the array spokes from wheeldata.
%Initialize spoke endpoints. x1 for x-coordinate hub end.

wHub = wheeldata.wHub;
rHub = wheeldata.rHub;
rRim = wheeldata.rRim;
nSpokes = wheeldata.nSpokes;

x1Spoke = wHub/2*ones(1,nSpokes);
sign = 1;
for n=1:nSpokes %Change sign for different side of hub
    x1Spoke(n) = x1Spoke(n)*sign;
    sign = sign*(-1);
end
%Does this even work? Yes it does!
x2Spoke = zeros(1,nSpokes); 
y1Spoke = rHub.*sin(linspace(0,2*pi,nSpokes+1));
y1Spoke(nSpokes+1) = [];
y2Spoke = rRim.*sin(linspace(0,2*pi,nSpokes+1));
y2Spoke(nSpokes+1) = [];
z1Spoke = rHub.*cos(linspace(0,2*pi,nSpokes+1));
z1Spoke(nSpokes+1) = [];
z2Spoke = rRim.*cos(linspace(0,2*pi,nSpokes+1));
z2Spoke(nSpokes+1) = [];
spokes = [x1Spoke' x2Spoke' y1Spoke' y2Spoke' z1Spoke' z2Spoke'];