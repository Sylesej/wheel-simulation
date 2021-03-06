function spokes=spokeCoordinates(wheeldata,u)
%This function defines the array spokes from wheeldata.
%Initialize spoke endpoints. x1 for x-coordinate hub end.
a = wheeldata.nSpokes*6;

uX = u(1:6:a);
uX = asin(uX/wheeldata.rRim)';
uX(end+1) = uX(1);
uY = u(2:6:a)';
uY(end+1) = uY(1);
uZ = u(3:6:a)';

wHub = wheeldata.wHub;
rHub = wheeldata.rHub;
rRim = wheeldata.rRim + uY;
nSpokes = wheeldata.nSpokes;

x1Spoke = wHub/2*ones(1,nSpokes);
sign = 1;
for n=1:nSpokes %Change sign for different side of hub
    x1Spoke(n) = x1Spoke(n)*sign;
    sign = sign*(-1);
end
x2Spoke = zeros(1,nSpokes)+uZ; 
y1Spoke = rHub.*sin(uX+linspace(0,2*pi,nSpokes+1));
y1Spoke(nSpokes+1) = [];
y2Spoke = rRim.*sin(uX+linspace(0,2*pi,nSpokes+1));
y2Spoke(nSpokes+1) = [];
z1Spoke = rHub.*cos(uX+linspace(0,2*pi,nSpokes+1));
z1Spoke(nSpokes+1) = [];
z2Spoke = rRim.*cos(uX+linspace(0,2*pi,nSpokes+1));
z2Spoke(nSpokes+1) = [];
spokes = [x1Spoke' x2Spoke' y1Spoke' y2Spoke' z1Spoke' z2Spoke'];