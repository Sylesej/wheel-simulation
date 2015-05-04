function [ K ] = stiffnessmatrix(wheeldata)
%STIFFNESSMATRIX Creats a stiffnessmatrix for a wheel with the given
%Bemï¿½rk, at det nu antages, at deformationerne er sï¿½ smï¿½, at de ikke har
%indvirkning pï¿½ vinklen ved det punkt, der ikke pï¿½virkes.
%parameters
nspokes = wheeldata.nSpokes;
npoints = wheeldata.nPoints;
rRim = wheeldata.rRim;
N = nspokes+npoints;
step = pointsFactor;
K = zeros(N*3, N*3);
Sigma = 2*pi/N;
EI = wheeldata.EIx;
EA = wheeldata.ESpokes*wheeldata.aSpokes;
k=2.1e+11;
l = 2*rRim*tan(pi/N); %Regner sidelï¿½ngden i hjulet
for n=1:N
    if(n==1)
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3, n*3) = 8*EI/l;
        K(N*3, n*3) = 2*EI/l;
        K(n*3+3, n*3) = -2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(N*3-2, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(N*3-1, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3) = -sin(Sigma/2)*6*EI/l^2;
        K(n*3+2, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 0;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(N*3-2, n*3-2) = cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(N*3-1, n*3-2) = sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3+1, n*3-2) = -cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(n*3+2, n*3-2) = -sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = 0; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
        %Reaktioner grundet forflytning i y-retning.
        K(N*3-2, n*3-2) = cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(N*3-1, n*3-2) = sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3+1, n*3-2) = -cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(n*3+2, n*3-2) = -sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = k; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
    elseif(n==N)
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3, n*3) = 8*EI/l;
        K(n*3-3, n*3) = 2*EI/l;
        K(3, n*3) = -2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(n*3-5, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3-4, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(3-2, n*3) = -sin(Sigma/2)*6*EI/l^2;
        K(3-1, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 0;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(n*3-5, n*3-2) = cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(n*3-4, n*3-2) = sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(3-2, n*3-2) = -cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(3-1, n*3-2) = -sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = 0; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
        %Reaktioner grundet forflytning i y-retning.
        K(n*3-5, n*3-2) = cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(n*3-4, n*3-2) = sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(3-2, n*3-2) = -cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(3-2, n*3-2) = -sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = k; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
    elseif(mod(n-1,step+1)==0)
<<<<<<< HEAD
        %Reaktionskræfterne i x på grund af vinkelrotationen
=======
        K(n*3-2, n*3-2) = 1;
        %Reaktionskrï¿½fterne i x pï¿½ grund af vinkelrotationen
>>>>>>> ba2fc3c03595302b565cd8487e88f16e74fe8c6d
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3, n*3) = 8*EI/l;
        K(n*3-3, n*3) = 2*EI/l;
        K(n*3+3, n*3) = -2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(n*3-5, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3-4, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3) = -sin(Sigma/2)*6*EI/l^2;
        K(n*3+2, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 0;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(n*3-5, n*3-2) = cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(n*3-4, n*3-2) = sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3+1, n*3-2) = -cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(n*3+2, n*3-2) = -sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = 0; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
        %Reaktioner grundet forflytning i y-retning.
        K(n*3-5, n*3-2) = cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(n*3-4, n*3-2) = sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3+1, n*3-2) = -cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(n*3+2, n*3-2) = -sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = k; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
    else
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3, n*3) = 8*EI/l;
        K(n*3-3, n*3) = 2*EI/l;
        K(n*3+3, n*3) = -2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(n*3-5, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3-4, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3) = -sin(Sigma/2)*6*EI/l^2;
        K(n*3+2, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 0;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(n*3-5, n*3-2) = cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(n*3-4, n*3-2) = sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3+1, n*3-2) = -cos(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        K(n*3+2, n*3-2) = -sin(Sigma/2)*(sin(Sigma/2)*EA/l+cos(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = 0; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
        %Reaktioner grundet forflytning i y-retning.
        K(n*3-5, n*3-2) = cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(n*3-4, n*3-2) = sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3+1, n*3-2) = -cos(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        K(n*3+2, n*3-2) = -sin(Sigma/2)*(cos(Sigma/2)*EA/l+sin(Sigma/2)*12*EI/l^3);
        
        K(n*3-2, n*3-2) = 0;
        K(n*3-1, n*3-2) = k; %Dette er det jeg er mest i tvivl om. At krï¿½fterne bliver 0 i punktet, der bliver flyttet?
        
    end
end
       

