function [thetaOut,newSpokes] = deformWeight(spokes,fWeight,wheeldata)
%This function calculates deformation of a wheel given the force and data
%about the wheel, and outputs changed spoke coordinates.
%Inputs:
%spokes is an array of coordinates of spoke endpoints.
%fWeight is a number specifying the amplitude of force. Assumed
%perpendicular to the rim.
%wheeldata is a struct containing data such as bending stiffness.
nSpokes = wheeldata.nSpokes;
EI = wheeldata.EIx;
forrige = 0;
ea = wheeldata.ESpokes*wheeldata.aSpokes;

l = 2*pi*wheeldata.rRim/nSpokes; %Can change
theta = zeros(1,nSpokes);
u =     zeros(1,nSpokes);
B =     zeros(1,nSpokes);
M =     zeros(1,nSpokes);
Fe =    ones(1,nSpokes).*wheeldata.t0;
B(1) = wheeldata.t0/2;
M(1) = 1;
nRuns = 0;
%nMangeRuns = 0;
while abs(sum(theta) - 2*pi) > 0.1 % Fejlen defineres.
    for n = 1:nSpokes
        if n == 1 %Første gennemkørsel, gravity
            u(1) = M(1)*l^2/(2*EI)+(B(1)-Fe(1)-fWeight)*l^3/(3*EI);
            theta(1) = M(1)*l/(EI)+(B(1)-Fe(1)-fWeight)*l^2/(2*EI);
            B(2) = B(1)-Fe(1)-fWeight;
            M(2) = M(1)+(-1*B(1)+Fe(1)+fWeight)*l;
        else
            %Følgende gennmkørsler
            u(n) = M(n)*l^2/(2*EI)+(B(n)-Fe(n))*l^3/(3*EI);
            theta(n) = M(n)*l/(EI)+(B(n)-Fe(n))*l^2/(2*EI);
            Fe(n) = Fe(n) - ea*u(n);
            B(n+1) = B(n)-Fe(n);
            M(n+1) = M(n)-B(n)*l+Fe(n)*l;
        end
    end
    M(1) = M(end);
    B(1) = B(end);
    nRuns= nRuns + 1;
    if(mod(nRuns,1e5)==0)
        display(['Antal gennemkørsler: ' num2str(nRuns)])
        display(['B er lig: ' num2str(B(1))])
        display(['M er lig: ' num2str(M(1))])
        display(['theta er lig: ' num2str(sum(theta))])
        u
        if sum(theta)>forrige
            display('SATANS!')
        end
        display(' ')
        forrige = sum(theta);
    end
end

nRuns
subplot(2,2,1)
plot(1:32,u(1:32))
title('Udbøjning')
subplot(2,2,2)
plot(1:32,theta)
title('Theta')
subplot(2,2,3)
plot(10:20,B(10:20))
title('B')
subplot(2,2,4)
plot(1:33,M)
title('M')
thetaOut = theta;
newSpokes = spokes;