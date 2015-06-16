function plotU(wheeldata,K,fy,fz)
% This function creates plots of deformations of the rim as well as forces
% in the spokes.
% It takes the arguments wheeldata (struct), the stiffnessmatrix of the
% wheel, K, as well as the forces fy and fz.

k= wheeldata.kSpokes;
theta = atan(wheeldata.wHub/wheeldata.rRim); % vinklen mellem egen og f�lg

[a,~] = size(K);
f = zeros(a,1); f(a/2-4) = fy; f(a/2-3) = fz; %Definer kr�fter

u = pinv(K)*f; % Deformation are found by solving the equation system

hold off
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
title 'Udb�jning i z-retning'

subplot(2,2,4)
u2 = (u(2:6:a)+u(3:6:a)*sin(theta))*k; %Forces in spokes are calculated
plot(u2)
title 'Eg-kr�fter'

% Spoke forces are projected onto main axis to check force equilibrium.
egF = zeros(length(u2),2);
sigma = 1:wheeldata.nSpokes;
sigma = sigma*2*pi/wheeldata.nSpokes;
%Antager meget lille forskydning af eger

egF(:,1) = sin(sigma).*u2'; % Spoke forces are projected 
egF(:,2) = cos(sigma).*u2';

display(['Summen af egekræfter i y: ' num2str(sum(egF(:,2)))])
display(['Summen af egekræfter i z: ' num2str(sum(egF(:,1)))])