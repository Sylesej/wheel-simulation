function plotU(wheeldata,K,fy,fz)
% This function creates plots of deformations of the rim as well as forces
% in the spokes.
% It takes the arguments wheeldata (struct), the stiffnessmatrix of the
% wheel, K, as well as the forces fy and fz.

rRim = wheeldata.rRim;
rHub = wheeldata.rHub;
wHub = wheeldata.wHub;
k= wheeldata.kSpokes;
theta = atan(wheeldata.wHub/wheeldata.rRim); % vinklen mellem egen og f�lg

[a,~] = size(K);
f = zeros(a,1); f(a/2-4) = fy; f(a/2-3) = fz; %Definer kr�fter

u = pinv(K)*f; % Deformation are found by solving the equation system

% -------------------------------------------------------------------------
% Plots are made!
% -------------------------------------------------------------------------
subplot(2,2,1)
plot(u(1:6:a))
title 'Udb�jning i x-retning'

subplot(2,2,2)
plot(u(2:6:a))
title 'Udb�jning i y-retning'

subplot(2,2,3)
plot(u(3:6:a))
title 'Udbøjning i z-retning'

subplot(2,2,4)
l0 = sqrt((rRim-rHub)^2+wHub^2);
ld = sqrt((rRim+u(2:6:a)-rHub).^2+(wHub+u(3:6:a)).^2); %pythagoras!
egpower =(l0-ld)*k; %pythagoras?
%Forces in spokes are calculated
plot(egpower)
title 'Eg-kr�fter'

% Spoke forces are projected onto main axis to check force equilibrium.
egF = zeros(length(egpower),2);
sigma = 1:wheeldata.nSpokes;
sigma = sigma*2*pi/wheeldata.nSpokes;
%Antager meget lille forskydning af eger

egF(:,1) = sin(sigma).*egpower'; % Spoke forces are projected 
egF(:,2) = cos(sigma).*egpower';

display(['Summen af egekræfter i y: ' num2str(sum(egF(:,2)))])