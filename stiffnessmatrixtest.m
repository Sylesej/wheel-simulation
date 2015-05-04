function [ K ] = stiffnessmatrixtest(pointsFactor, wheeldata)
%STIFFNESSMATRIX Creats a stiffnessmatrix for a wheel with the given
%Bemærk, at det nu antages, at deformationerne er så små, at de ikke har
%indvirkning på vinklen ved det punkt, der ikke påvirkes.
%parameters
nspokes = 4;
npoints = pointsFactor*4;
rRim = wheeldata.rRim;
N = nspokes+npoints;
step = pointsFactor;
K = zeros(N*3, N*3);
Sigma = 2*pi/N;
EI = 1;
EA = 1;
k=2.1e+11;
l = 2*rRim*tan(pi/N); %Regner sidelængden i hjulet
for(n=1:N)
    if(n==1)
        K(n*3-2,n*3-2) = 1;
    elseif(n==N)
        K(n*3-2,n*3-2) = 3;
    elseif(mod(n-1,step+1)==0)
        K(n*3-2,n*3-2) = 2;
    else
       K(n*3-2,n*3-2) = 4;
        
    end
end