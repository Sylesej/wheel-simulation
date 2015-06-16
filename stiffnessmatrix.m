function [ K ] = stiffnessmatrix(wheeldata)
% STIFFNESSMATRIX Creates a stiffnessmatrix for a wheel with the properties
% givin in a wheeldata struct.
% Deformations are assumed small.

% Properties are loaded from wheeldata struct
N = wheeldata.nSpokes; %Number of spokes
rRim = wheeldata.rRim; %Radius of rim
EIz = wheeldata.EIz; %Moments of inertia
EIy = wheeldata.EIy;
EA = wheeldata.EA; %Tensile moment
k=wheeldata.kSpokes; %Stiffness of spokes
l = 2*rRim*sin(pi/N); %Length of a beam in the N-polygon.

K = zeros(N*6, N*6); %Stiffnessmatrix is initialized
Sigma = 2*pi/N; %Angle between two adjunctant beams
theta = atan(wheeldata.wHub/wheeldata.rRim); % Angle between spoke and rim

for n=1:N
    if(n==1)
        %Starten er kun i x-y og alle z-kr�fter og momenter bliver 0.
        %___________________________________________________________%
        %Momenterne om z grundet vinkelrotationen om z-aksen i et punkt.
        K(n*6, n*6) = 8*EIz/l;   %Punkt(n)
        K(N*6, n*6) = 2*EIz/l;   %Punkt(n-1)
        K((n+1)*6, n*6) = 2*EIz/l; %punkt(n+1)
        
        %Reaktionskr�fterne i system p� grund af vinkeldrejningen
        %om z-aksen.
        
        %Punkt(n-1)
        K(N*6-5, n*6) = sin(Sigma/2)*6*EIz/l^2;  %x
        K(N*6-4, n*6) = cos(Sigma/2)*6*EIz/l^2;  %y
        K(N*6-3, n*6) = 0;                      %z
        
        %Punkt(n+1)
        K((n+1)*6-5, n*6) = sin(Sigma/2)*6*EIz/l^2;  %x
        K((n+1)*6-4, n*6) = -cos(Sigma/2)*6*EIz/l^2; %y
        K((n+1)*6-3, n*6) = 0;                      %z
        
        %Punkt(n)
        K(n*6-5, n*6) = 2*sin(Sigma/2)*6*EIz/l^2;    %x
        K(n*6-4, n*6) = 0;                          %y
        K(n*6-3, n*6) = 0;                          %z
        
        %Reaktioner grundet forflytning i x-retning.
        %Punktet inden (n-1)
        K(N*6-5, n*6-5) =...%x
            -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EIz/l^3;
        K(N*6-4, n*6-5) =...%y
            sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EIz/l^3);
        K(N*6-3, n*6-5) = 0; %z
        
        K(N*6  , n*6-5) = sin(Sigma/2)*6*EIz/l^2;  %Moment om Z
        
        %Punktet efter (n+1)
        K((n+1)*6-5, n*6-5) =...%x
            sin(Sigma/2)^2*12*EIz/l^3-cos(Sigma/2)^2*EA/l;
        K((n+1)*6-4, n*6-5) =...%y
            -cos(Sigma/2)*sin(Sigma/2)*(12*EIz/l^3+EA/l);
        
        K((n+1)*6, n*6-5) = sin(Sigma/2)*6*EIz/l^2; %Moment om Z
        
        %Selve punktet (n)
        K(n*6-5, n*6-5) =...%x
            2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EIz/l^3);
        K(n*6-4, n*6-5) = 0;%y
        K(n*6-0, n*6-5) = 2*sin(Sigma/2)*6*EIz/l^2;%Moment om Z
        
        %Reaktioner grundet forflytning i y-retning.
        %Punktet inden (n-1)
        K(N*6-5, n*6-4) =...%x
            -cos(Sigma/2)*sin(Sigma/2)*(EA/l+12*EIz/l^3);
            %-(cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l); %Fejl?
        K(N*6-4, n*6-4) =...%y
            -cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;
            %-sin(Sigma/2)*cos(Sigma/2)*(12*EIz/l^3+EA/l); %Fejl?
        K(N*6-0, n*6-4) = -cos(Sigma/2)*6*EIz/l^2; %moment om z
        
        %Punktet efter (n+1)
        K((n+1)*6-5, n*6-4) =...%x
            cos(Sigma/2)*sin(Sigma/2)*(EA/l+12*EIz/l^3);
            %cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;%Fejl?
        
        K((n+1)*6-4, n*6-4) =...%y
            -cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;
            %-sin(Sigma/2)*cos(Sigma/2)*(12*EIz/l^3+EA/l); %Fejl?
        K((n+1)*6-0, n*6-4) = cos(Sigma/2)*6*EIz/l^2; %Moment om z
        
        %Punktet selv(n)
        K(n*6-5, n*6-4) = 0;%x
        K(n*6-4, n*6-4) = ...%y
            2*(cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l)+k;
        K(n*6-0, n*6-4) = 0;%Moment om z
        
        %Reaktioner grundet forflytning i z-retningen z-asksen regnes
        %positiv "udad".
        
        %Punktet inden (n-1)
        
        K(N*6-3,n*6-3) = -12*EIy/l^3; % z
%         K(N*6-2,n*6-3) = sin(Sigma/2)*6*EIy/l^2; 
% moment om x antages til 0
        K(N*6-1,n*6-3) = -6*EIy/l^2; % moment om y
        K(N*6-0,n*6-3) = 0; % moment om z
        
        %Punktet efter (n+1)
        
        K((n+1)*6-3,n*6-3) = -12*EIy/l^3; % z
%         K((n+1)*6-2,n*6-3) = sin(Sigma/2)*6*EIy/l^2; 
% moment om x antages til 0
        K((n+1)*6-1,n*6-3) = 6*EIy/l^2; % moment om y
        K((n+1)*6-0,n*6-3) = 0; % moment om z
        
        %Punktet selv
        
        K(n*6-3, n*6-3) = 2*12*EIy/l^3 + k*sin(theta); %z
%         K(n*6-2, n*6-3) = 2*sin(Sigma/2)*6*EIy/l^2; 
%moment om x antages til 0
        K(n*6-1, n*6-3) = 0; %moment om y
        K(n*6-0, n*6-3) = 0; %moment om z
        
        %Reaktioner grundet vinkeldrejning om x-aksen
        
%         %Punktet inden
%         K(N*6-3, n*6-2) = sin(Sigma/2)*6*EIy/l^2; %z
%         K(N*6-2, n*6-2) = 
%         K(N*6-1, n*6-2) =
%         K(N*6-0, n*6-2) = 
%                
        %Punktet efter
        
        %Punktet selv
        
        %Reaktioner grundet vinkeldrejning om y-aksen
        
        
        %Punktet inden
        K(N*6-3, n*6-1) = cos(Sigma/2)*6*EIy/l^2; %z
        %K(N*6-2, n*6-1) = ... %Moment om x
        %    2*sin(Sigma/2)*cos(Sigma/2)*2*EIy/l;
        K(N*6-1, n*6-1) = ... %Moment om y
            (cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K(N*6-0, n*6-1) = 0; %Moment om z
        
        %Punktet efter
        K((n+1)*6-3, n*6-1) = -cos(Sigma/2)*6*EIy/l^2; %z
        %K((n+1)*6-2, n*6-1) = ... %moment om x
        %    -2*cos(Sigma/2)*sin(Sigma/2)*2*EIy/l;
        K((n+1)*6-1, n*6-1) = ... %moment om y
            (cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K((n+1)*6-0, n*6-1) = 0;
        
        %Punktet selv
        K(n*6-3, n*6-1) = 0;
        K(n*6-2, n*6-1) = 0;
        K(n*6-1, n*6-1) = 2*4*EIy/l;
        K(n*6-0, n*6-1) = 0;
        
    elseif(n==N)
        %Starten er kun i x-y og alle z-kr�fter og momenter bliver 0.
        %___________________________________________________________%
        %Momenterne om z grundet vinkelrotationen om z-aksen i et punkt.
        K(n*6, n*6) = 8*EIz/l;   %Punkt(n)
        K((n-1)*6, n*6) = 2*EIz/l;   %Punkt(n-1)
        K((1)*6, n*6) = 2*EIz/l; %punkt(n+1)
        
        %Reaktionskr�fterne i system p� grund af vinkeldrejningen
        %om z-aksen.
        
        %Punkt(n-1)
        K((n-1)*6-5, n*6) = sin(Sigma/2)*6*EIz/l^2;  %x
        K((n-1)*6-4, n*6) = cos(Sigma/2)*6*EIz/l^2;  %y
        K((n-1)*6-3, n*6) = 0;                      %z
        
        %Punkt(n+1)
        K((1)*6-5, n*6) = sin(Sigma/2)*6*EIz/l^2;  %x
        K((1)*6-4, n*6) = -cos(Sigma/2)*6*EIz/l^2; %y
        K((1)*6-3, n*6) = 0;                      %z
        
        %Punkt(n)
        K(n*6-5, n*6) = 2*sin(Sigma/2)*6*EIz/l^2;    %x
        K(n*6-4, n*6) = 0;                          %y
        K(n*6-3, n*6) = 0;                          %z
        
        %Reaktioner grundet forflytning i x-retning.
        %Punktet inden (n-1)
        K((n-1)*6-5, n*6-5) =...%x
            -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EIz/l^3;
        K((n-1)*6-4, n*6-5) =...%y
            sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EIz/l^3);
        K((n-1)*6-3, n*6-5) = 0; %z
        
        K((n-1)*6  , n*6-5) = sin(Sigma/2)*6*EIz/l^2;  %Moment om Z
        
        %Punktet efter (n+1)
        K((1)*6-5, n*6-5) =...%x
            sin(Sigma/2)^2*12*EIz/l^3-cos(Sigma/2)^2*EA/l;
        K((1)*6-4, n*6-5) =...%y
            cos(Sigma/2)*sin(Sigma/2)*(-12*EIz/l^3-EA/l);
        
        K((1)*6, n*6-5) = sin(Sigma/2)*6*EIz/l^2; %Moment om Z
        
        %Selve punktet (n)
        K(n*6-5, n*6-5) =...%x
            2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EIz/l^3);
        K(n*6-4, n*6-5) = 0;%y
        K(n*6-0, n*6-5) = 2*sin(Sigma/2)*6*EIz/l^2;%Moment om Z
        
        %Reaktioner grundet forflytning i y-retning.
        %Punktet inden (n-1)
        K((n-1)*6-5, n*6-4) =...%x
            -cos(Sigma/2)*sin(Sigma/2)*(EA/l+12*EIz/l^3);
        K((n-1)*6-4, n*6-4) =...%y
            -cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;
        K((n-1)*6-0, n*6-4) = -cos(Sigma/2)*6*EIz/l^2; %moment om z
        
        %Punktet efter (n+1)
        K((1)*6-5, n*6-4) =...%x
            cos(Sigma/2)*sin(Sigma/2)*(EA/l+12*EIz/l^3);
        K((1)*6-4, n*6-4) =...%y
            -cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;
        K((1)*6-0, n*6-4) = cos(Sigma/2)*6*EIz/l^2; %Moment om z
        %Punktet selv(n)
        K(n*6-5, n*6-4) = 0;%x
        K(n*6-4, n*6-4) = ...%y
            2*(cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l)+k;
        K(n*6-0, n*6-4) = 0;%Moment om z
        
        
        %Reaktioner grundet forflytning i z-retningen z-asksen regnes
        %positiv "udad".
        
        %Punktet inden (n-1)
        
        K((n-1)*6-3,n*6-3) = -12*EIy/l^3; % z
%         K(N*6-2,n*6-3) = sin(Sigma/2)*6*EIy/l^2; 
% moment om x antages til 0
        K((n-1)*6-1,n*6-3) = -6*EIy/l^2; % moment om y
        K((n-1)*6-0,n*6-3) = 0; % moment om z
        
        %Punktet efter (n+1)
        
        K((1)*6-3,n*6-3) = -12*EIy/l^3; % z
%         K((n+1)*6-2,n*6-3) = sin(Sigma/2)*6*EIy/l^2; 
% moment om x antages til 0
        K((1)*6-1,n*6-3) = 6*EIy/l^2; % moment om y
        K((1)*6-0,n*6-3) = 0; % moment om z
        
        %Punktet selv
        
        K(n*6-3, n*6-3) = 2*12*EIy/l^3 + k*sin(theta); %z
%         K(n*6-2, n*6-3) = 2*sin(Sigma/2)*6*EIy/l^2; 
%moment om x antages til 0
        K(n*6-1, n*6-3) = 0; %moment om y
        K(n*6-0, n*6-3) = 0; %moment om z
        
        %Reaktioner grundet vinkeldrejning om x-aksen
        
%         %Punktet inden
%         K(N*6-3, n*6-2) = sin(Sigma/2)*6*EIy/l^2; %z
%         K(N*6-2, n*6-2) = 
%         K(N*6-1, n*6-2) =
%         K(N*6-0, n*6-2) = 
%                
        %Punktet efter
        
        %Punktet selv
        
        %Reaktioner grundet vinkeldrejning om y-aksen
        
        
        %Punktet inden
        K((n-1)*6-3, n*6-1) = cos(Sigma/2)*6*EIy/l^2; %z
        %K((n-1)*6-2, n*6-1) = ... %Moment om x
        %    2*sin(Sigma/2)*cos(Sigma/2)*2*EIy/l;
        K((n-1)*6-1, n*6-1) = ... %Moment om y
            (cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K((n-1)*6-0, n*6-1) = 0; %Moment om z
        
        %Punktet efter
        K((1)*6-3, n*6-1) = -cos(Sigma/2)*6*EIy/l^2; %z
        %K((1)*6-2, n*6-1) = ... %moment om x
        %    -2*cos(Sigma/2)*sin(Sigma/2)*2*EIy/l;
        K((1)*6-1, n*6-1) = ... %moment om y
            (cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K((1)*6-0, n*6-1) = 0;
        
        %Punktet selv
        K(n*6-3, n*6-1) = 0;
        K(n*6-2, n*6-1) = 0;
        K(n*6-1, n*6-1) = 2*4*EIy/l;
        K(n*6-0, n*6-1) = 0;
        
    else  %if(mod(n-1,step+1)==0)
        
        %Starten er kun i x-y og alle z-kr�fter og momenter bliver 0.
        %___________________________________________________________%
        %Momenterne om z grundet vinkelrotationen om z-aksen i et punkt.
        K(n*6, n*6) = 8*EIz/l;   %Punkt(n)
        K((n-1)*6, n*6) = 2*EIz/l;   %Punkt(n-1)
        K((n+1)*6, n*6) = 2*EIz/l; %punkt(n+1)
        
        %Reaktionskr�fterne i system p� grund af vinkeldrejningen
        %om z-aksen.
        
        %Punkt(n-1)
        K((n-1)*6-5, n*6) = sin(Sigma/2)*6*EIz/l^2;  %x
        K((n-1)*6-4, n*6) = cos(Sigma/2)*6*EIz/l^2;  %y
        K((n-1)*6-3, n*6) = 0;                      %z
        
        %Punkt(n+1)
        K((n+1)*6-5, n*6) = sin(Sigma/2)*6*EIz/l^2;  %x
        K((n+1)*6-4, n*6) = -cos(Sigma/2)*6*EIz/l^2; %y
        K((n+1)*6-3, n*6) = 0;                      %z
        
        %Punkt(n)
        K(n*6-5, n*6) = 2*sin(Sigma/2)*6*EIz/l^2;    %x
        K(n*6-4, n*6) = 0;                          %y
        K(n*6-3, n*6) = 0;                          %z
        
        %Reaktioner grundet forflytning i x-retning.
        %Punktet inden (n-1)
        K((n-1)*6-5, n*6-5) =...%x
            -cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EIz/l^3;
        K((n-1)*6-4, n*6-5) =...%y
            sin(Sigma/2)*cos(Sigma/2)*(EA/l+12*EIz/l^3);
        K((n-1)*6-3, n*6-5) = 0; %z
        
        K((n-1)*6  , n*6-5) = sin(Sigma/2)*6*EIz/l^2;  %Moment om Z
        
        %Punktet efter (n+1)
        K((n+1)*6-5, n*6-5) =...%x
            sin(Sigma/2)^2*12*EIz/l^3-cos(Sigma/2)^2*EA/l;
        K((n+1)*6-4, n*6-5) =...%y
            cos(Sigma/2)*sin(Sigma/2)*(-12*EIz/l^3-EA/l);
        
        K((n+1)*6, n*6-5) = sin(Sigma/2)*6*EIz/l^2; %Moment om Z
        
        %Selve punktet (n)
        K(n*6-5, n*6-5) =...%x
            2*(cos(Sigma/2)^2*EA/l+sin(Sigma/2)^2*12*EIz/l^3);
        K(n*6-4, n*6-5) = 0;%y
        K(n*6-0, n*6-5) = 2*sin(Sigma/2)*6*EIz/l^2;%Moment om Z
        
        %Reaktioner grundet forflytning i y-retning.
        %Punktet inden (n-1)
        K((n-1)*6-5, n*6-4) =...%x
            -cos(Sigma/2)*sin(Sigma/2)*(EA/l+12*EIz/l^3);
        K((n-1)*6-4, n*6-4) =...%y
            -cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;
        K((n-1)*6-0, n*6-4) = -cos(Sigma/2)*6*EIz/l^2; %moment om z
        
        %Punktet efter (n+1)
        K((n+1)*6-5, n*6-4) =...%x
            cos(Sigma/2)*sin(Sigma/2)*(EA/l+12*EIz/l^3);
        K((n+1)*6-4, n*6-4) =...%y
            -cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l;           
        K((n+1)*6-0, n*6-4) = cos(Sigma/2)*6*EIz/l^2; %Moment om z
        %Punktet selv(n)
        K(n*6-5, n*6-4) = 0;%x
        K(n*6-4, n*6-4) = ...%y
            2*(cos(Sigma/2)^2*12*EIz/l^3+sin(Sigma/2)^2*EA/l)+k;
        K(n*6-0, n*6-4) = 0;%Moment om z
        
        %Reaktioner grundet forflytning i z-retningen z-asksen regnes
        %positiv "udad".
        
        %Punktet inden (n-1)
        
        K((n-1)*6-3,n*6-3) = -12*EIy/l^3; % z
%         K(N*6-2,n*6-3) = sin(Sigma/2)*6*EIy/l^2; 
% moment om x antages til 0
        K((n-1)*6-1,n*6-3) = -6*EIy/l^2; % moment om y
        K((n-1)*6-0,n*6-3) = 0; % moment om z
        
        %Punktet efter (n+1)
        
        K((n+1)*6-3,n*6-3) = -12*EIy/l^3; % z
%         K((n+1)*6-2,n*6-3) = sin(Sigma/2)*6*EIy/l^2; 
% moment om x antages til 0
        K((n+1)*6-1,n*6-3) = 6*EIy/l^2; % moment om y
        K((n+1)*6-0,n*6-3) = 0; % moment om z
        
        %Punktet selv
        
        K(n*6-3, n*6-3) = 2*12*EIy/l^3 + k*sin(theta); %z
%         K(n*6-2, n*6-3) = 2*sin(Sigma/2)*6*EIy/l^2; 
%moment om x antages til 0
        K(n*6-1, n*6-3) = 0; %moment om y
        K(n*6-0, n*6-3) = 0; %moment om z
        
        %Reaktioner grundet vinkeldrejning om x-aksen
        
%         %Punktet inden
%         K(N*6-3, n*6-2) = sin(Sigma/2)*6*EIy/l^2; %z
%         K(N*6-2, n*6-2) = 
%         K(N*6-1, n*6-2) =
%         K(N*6-0, n*6-2) = 
%                
        %Punktet efter
        
        %Punktet selv
        
        %Reaktioner grundet vinkeldrejning om y-aksen
        
        
        %Punktet inden
        K((n-1)*6-3, n*6-1) = cos(Sigma/2)*6*EIy/l^2; %z
        %K((n-1)*6-2, n*6-1) = ... %Moment om x
        %    2*sin(Sigma/2)*cos(Sigma/2)*2*EIy/l;
        K((n-1)*6-1, n*6-1) = ... %Moment om y
            (cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K((n-1)*6-0, n*6-1) = 0; %Moment om z
        
        %Punktet efter
        K((n+1)*6-3, n*6-1) = -cos(Sigma/2)*6*EIy/l^2; %z
        %K((n+1)*6-2, n*6-1) = ... %moment om x
        %    -2*cos(Sigma/2)*sin(Sigma/2)*2*EIy/l;
        K((n+1)*6-1, n*6-1) = ... %moment om y
            (cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K((n+1)*6-0, n*6-1) = 0;
        
        %Punktet selv
        K(n*6-3, n*6-1) = 0;
        K(n*6-2, n*6-1) = 0;
        K(n*6-1, n*6-1) = 2*(cos(Sigma/2)^2-sin(Sigma/2)^2)*2*EIy/l;
        K(n*6-0, n*6-1) = 0;
    end
end