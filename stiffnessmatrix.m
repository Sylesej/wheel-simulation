function [ K ] = stiffnessmatrix(wheeldata)
%STIFFNESSMATRIX Creats a stiffnessmatrix for a wheel with the given
%Bemï¿½rk, at det nu antages, at deformationerne er sï¿½ smï¿½, at de ikke har
%indvirkning pï¿½ vinklen ved det punkt, der ikke pï¿½virkes.
%parameters
nspokes = wheeldata.nSpokes;
npoints = wheeldata.nPoints;
rRim = wheeldata.rRim;
N = nspokes+npoints;
step = wheeldata.pointFactor;
K = zeros(N*3, N*3);
Sigma = 2*pi/N;
EI = wheeldata.EIx;
EA = wheeldata.ESpokes*wheeldata.aSpokes;
k= wheeldata.ESpokes*wheeldata.aSpokes;
l = 2*rRim*sin(pi/N); %Regner sidelï¿½ngden i hjulet
for n=1:N
    if(n==1)
        %Starten er kun i x-y og alle z-kræfter og momenter bliver 0.
        %___________________________________________________________%
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3, n*3) = 8*EI/l;
        K(N*3, n*3) = 2*EI/l;
        K(n*3+3, n*3) = 2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(N*3-2, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(N*3-1, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3+2, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 2*sin(Sigma/2)*6*EI/l^2;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(N*3-2, n*3-2) = -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3;
        K(N*3-1, n*3-2) = sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EI/l^3);
        K(N*3  , n*3-2) = sin(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3-2) = sin(Sigma/2)^2*12*EI/l^3-cos(Sigma/2)^2*EA/l;
        K(n*3+2, n*3-2) = cos(Sigma/2)*sin(Sigma/2)*(-12*EI/l^3-EA/l);
        K(n*3+3, n*3-2) = sin(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3-2) = 2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3);
        K(n*3-1, n*3-2) = 0;
        K(n*3  , n*3-2) = 2*sin(Sigma/2)*6*EI/l^2;
        
        %Reaktioner grundet forflytning i y-retning.
        K(N*3-2, n*3-1) = -(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l);
        K(N*3-1, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
        K(N*3  , n*3-1) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3-1) = cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l;
        K(n*3+2, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
        K(n*3+3, n*3-1) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3-1) = 0;
        K(n*3-1, n*3-1) = ...
            2*(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l)+k;
        K(n*3  , n*3-1) = 0;
              
    elseif(n==N)
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3, n*3) = 8*EI/l;
        K(n*3-3, n*3) = 2*EI/l;
        K(3, n*3) = 2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(n*3-5, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3-4, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(3-2, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(3-1, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 2*sin(Sigma/2)*6*EI/l^2;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(n*3-5, n*3-2) = -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3;
        K(n*3-4, n*3-2) = sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EI/l^3);
        K(n*3-3  , n*3-2) = sin(Sigma/2)*6*EI/l^2;
        
        K(3-2, n*3-2) = sin(Sigma/2)^2*12*EI/l^3-cos(Sigma/2)^2*EA/l;
        K(3-1, n*3-2) = cos(Sigma/2)*sin(Sigma/2)*(-12*EI/l^3-EA/l);
        K(3-0, n*3-2) = sin(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3-2) = 2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3);
        K(n*3-1, n*3-2) = 0;
        K(n*3  , n*3-2) = 2*sin(Sigma/2)*6*EI/l^2;
        
        %Reaktioner grundet forflytning i y-retning.
        K(n*3-5, n*3-1) = -(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l);
        K(n*3-4, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
        K(n*3-3, n*3-1) = -cos(Sigma/2)*6*EI/l^2;
        
        K(3-2, n*3-1) = cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l;
        K(3-1, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
        K(3-0  , n*3-1) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3-1) = 0;
        K(n*3-1, n*3-1) = ...
            2*(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l)+k;
        K(n*3-0, n*3-1) = 0;
        
    else  %if(mod(n-1,step+1)==0)
        
        %Reaktionskræfterne i x på grund af vinkelrotationen
        %Momenterne grundet vinkelrotationen i et punkt.
        K(n*3  , n*3) = 8*EI/l;
        K(n*3-3, n*3) = 2*EI/l;
        K(n*3+3, n*3) = 2*EI/l;
        
        %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
        K(n*3-5, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3-4, n*3) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3) = sin(Sigma/2)*6*EI/l^2;
        K(n*3+2, n*3) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3) = 2*sin(Sigma/2)*6*EI/l^2;
        K(n*3-1, n*3) = 0;
        
        %Reaktioner grundet forflytning i x-retning.
        K(n*3-5, n*3-2) = -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3;
        K(n*3-4, n*3-2) = sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EI/l^3);
        K(n*3-3, n*3-2) = sin(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3-2) = sin(Sigma/2)^2*12*EI/l^3-cos(Sigma/2)^2*EA/l;
        K(n*3+2, n*3-2) = cos(Sigma/2)*sin(Sigma/2)*(-12*EI/l^3-EA/l);
        K(n*3+3, n*3-2) = sin(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3-2) = 2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3);
        K(n*3-1, n*3-2) = 0;
        K(n*3-0, n*3-2) = 2*sin(Sigma/2)*6*EI/l^2;
        
        %Reaktioner grundet forflytning i y-retning.
        K(n*3-5, n*3-1) = -(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l);
        K(n*3-4, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
        K(n*3-3, n*3-1) = -cos(Sigma/2)*6*EI/l^2;
        
        K(n*3+1, n*3-1) = cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l;
        K(n*3+2, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
        K(n*3+3, n*3-1) = cos(Sigma/2)*6*EI/l^2;
        
        K(n*3-2, n*3-1) = 0;
        K(n*3-1, n*3-1) =...
            2*(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l)+k;
        K(n*3-0, n*3-1) = 0;
%     else
%         %Momenterne grundet vinkelrotationen i et punkt.
%         K(n*3, n*3) = 8*EI/l;
%         K(n*3-3, n*3) = 2*EI/l;
%         K(n*3+3, n*3) = 2*EI/l;
%         
%         %Reaktionskrï¿½fterne pï¿½ grund af vinkeldrejningen.
%         K(n*3-5, n*3) = sin(Sigma/2)*6*EI/l^2;
%         K(n*3-4, n*3) = cos(Sigma/2)*6*EI/l^2;
%         
%         K(n*3+1, n*3) = sin(Sigma/2)*6*EI/l^2;
%         K(n*3+2, n*3) = -cos(Sigma/2)*6*EI/l^2;
%         
%         K(n*3-2, n*3) = 2*sin(Sigma/2)*6*EI/l^2;
%         K(n*3-1, n*3) = 0;
%         
%         %Reaktioner grundet forflytning i x-retning.
%         K(n*3-5, n*3-2) = -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3;
%         K(n*3-4, n*3-2) = sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EI/l^3);
%         K(n*3-3, n*3-2) = sin(Sigma/2)*6*EI/l^2;
%         
%         K(n*3+1, n*3-2) = sin(Sigma/2)^2*12*EI/l^3-cos(Sigma/2)^2*EA/l;
%         K(n*3+2, n*3-2) = cos(Sigma/2)*sin(Sigma/2)*(-12*EI/l^3-EA/l);
%         K(n*3+3, n*3-2) = sin(Sigma/2)*6*EI/l^2;
%         
%         K(n*3-2, n*3-2) = 2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EI/l^3);
%         K(n*3-1, n*3-2) = 0;
%         K(n*3-0, n*3-2) = 2*sin(Sigma/2)*6*EI/l^2;
%         
%         %Reaktioner grundet forflytning i y-retning.
%         K(n*3-5, n*3-1) = -(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l);
%         K(n*3-4, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
%         K(n*3-3, n*3-1) = -cos(Sigma/2)*6*EI/l^2;
%         
%         K(n*3+1, n*3-1) = cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l;
%         K(n*3+2, n*3-1) = -sin(Sigma/2)*cos(Sigma/2)*(12*EI/l^3+EA/l);
%         K(n*3+3, n*3-1) = cos(Sigma/2)*6*EI/l^2;
%         
%         K(n*3-2, n*3-1) = 0;
%         K(n*3-1, n*3-1) = ...
%             2*(cos(Sigma/2)^2*12*EI/l^3+sin(Sigma/2)^2*EA/l)+k;
%         K(n*3-0, n*3-1) = 0;
    end
end