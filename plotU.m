function plotU(wheeldata,K,fy,fz)

[a,~] = size(K);

f = zeros(a,1); f(a/2-4) = fy; f(a/2-3) = fz; %Definer kr�fter

u = pinv(K)*f;

k= wheeldata.kSpokes;

theta = atan(wheeldata.wHub/wheeldata.rRim); % vinklen mellem egen og f�lg

u2 = (u(2:6:a)+u(3:6:a)*sin(theta))*k;

hold off

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
plot(u2)
title 'Eg-kr�fter'

egF = zeros(length(u2),2);
sigma = 1:wheeldata.nSpokes;
sigma = sigma*2*pi/wheeldata.nSpokes;
%Antager meget lille forskydning af eger

egF(:,1) = sin(sigma).*u2';
display('Egekr�fter i x og y:')
egF(:,2) = cos(sigma).*u2'
display('Summen af egekr�fter i y')
sum(egF(:,2))