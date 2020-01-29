function [CD,CL,CN,CA]=sentman(alph,Ti,Tw,accom,epsil,V,m,shape,S)
%calculates the coefficient of drag on a sphere or a plate assuming completely diffuse
%reflections
%    epsil   =   0 means specular fraction
%S is the speed ratio, if S is -1, the speed ratio is computed from other
%parameters
%cylider option is for an un-scaled cylinder


kb          =   1.3806503e-23;

%speed coefficient
% V           =   norm(V);
if S==-1
    S       =   norm(V)/sqrt(2*kb*Ti/m);
end

Tr          =   (m/(3*kb))*V^2*(1-accom)+accom*Tw;

CD          =   0;

%%
%surface (plate with ONE SIDE exposed to flow no specular reflection, Sentmans equations)
if strcmp(shape,'surface')%<<<<<<<<<<<<<Sentman
    %for these equations, alpha is defined as the angle between the normal
    %and the flow
    %Here
    %CA is the force "axial" or normal to the plate
    %CN is the force along the plate
    CA      =   ((cos(alph)^2)+1/(2*S^2))*(1+erf(S*cos(alph)))+cos(alph)*exp(-S^2*cos(alph)^2)/(S*sqrt(pi))+...
                sqrt(Tr/Ti)*((sqrt(pi)/(2*S))*cos(alph)*(1+erf(S*cos(alph)))+exp(-S^2*cos(alph)^2)/(2*S^2));
    CN      =   sin(alph)*cos(alph)*(1+erf(S*cos(alph)))+sin(alph)*exp(-S^2*cos(alph)^2)/(S*sqrt(pi));
    CD      =   CN*sin(alph)+CA*cos(alph);
    CL      =   CN*cos(alph)-CA*sin(alph);
end

%%
%surface (plate with ONE SIDE exposed to flow specular reflection included, Bird 1994 p 165)
if strcmp(shape,'surfacespec')%<<<<<<<<<<<Bird
    %Bird p 165 eq 7.58-7.60 and eq 7.61-7.62
    %transform to Brirds coordinates (angle of attack)
    alph    =   -pi/2-alph;
    %Here:
    %CN denotes the force normal to the plate (note the reversal from Sentman's notation)
    %CA is the force parallel (along the plate)
    CN      =   -(((1+epsil)*S*sin(alph)/sqrt(pi)+0.5*(1-epsil)*sqrt(Tr/Ti))*exp(-S^2*sin(alph)^2)...
                +((1+epsil)*(0.5+S^2*sin(alph)^2)+(0.5*(1-epsil)*sqrt(Tr/Ti)*sqrt(pi)*S*sin(alph)))*(1+erf(S*sin(alph))))/S^2;
    CA      =   -(((1-epsil)*S*cos(alph)/sqrt(pi))*(exp(-S^2*sin(alph)^2)+sqrt(pi)*S*sin(alph)*(1+erf(S*sin(alph)))))/S^2;
    CD      =   abs(CN*sin(alph)+CA*cos(alph));
    CL      =   CN*cos(alph)-CA*sin(alph);
end

%%
%plate with BOTH SIDES exposed to flow, Bird 1994 p 170
if strcmp(shape,'plate')
    CN      =   0;
    CA      =   0;
    %Bird p 170
    %transform to Brirds coordinates (Bird uses angle of attack)
    alph    =   pi/2-alph;
    CD      =   (2*(1-epsil*cos(2*alph))/(sqrt(pi)*S))*exp(-S^2*sin(alph)^2)+(sin(alph)/S^2)*(1+2*S^2+epsil*(1-2*S^2*cos(2*alph)))*erf(S*sin(alph))+((1-epsil)/S)*sqrt(pi)*sin(alph)^2*sqrt(Tr/Ti);
    CL      =   (4*epsil/(sqrt(pi)*S))*sin(alph)*cos(alph)*exp(-S^2*sin(alph)^2)+(cos(alph)/S^2)*(1+epsil*(1+4*S^2*sin(alph)^2))*erf(S*sin(alph))+((1-epsil)/S)*sqrt(pi)*sin(alph)*cos(alph)*sqrt(Tr/Ti);
end

%%
%SPHERE
if strcmp(shape,'sphere')

    CN      =   0;
    CA      =   0;
    
    %Bird p.172 eq. 7.71
    CD      =   ((2*S^2+1)/(sqrt(pi)*S^3))*exp(-S^2)+((4*S^4+4*S^2-1)/(2*S^4))*erf(S)+((2*(1-epsil)*sqrt(pi))/(3*S))*sqrt(Tr/Ti);
    CL      =   0;
end

%%
%CYLINDER, NO END CAPS - from Sentman
if strcmp(shape,'cylinder')
    %answer must be scaled to length times radius
    %alph is the pitch angle away from the axis of revolution
    Spar    =   S^2*sin(alph)^2/2;
    scl     =   0;%bessel function scaling
    %Aref    =   obliqueCylProjection(D,L,alph);%<<<<<<<(pi*D^2/4)*cos(alph)+(D*L)*sin(alph);%<<<<<<obliqueCylProjection(D,L,alph)
    SF      =   1;%(D/2)*L/Aref;
    CN      =   SF*(S*sqrt(pi)*sin(alph)*(2*sin(alph)^2+1/S^2)*exp(-Spar)*(besseli(0,Spar,scl)+besseli(1,Spar,scl))+(2*sqrt(pi)/S)*sin(alph)*exp(-Spar)*besseli(0,Spar,scl)+sqrt(Tr/Ti)*((pi^1.5)/(2*S))*sin(alph));
    CA      =   SF*(2*S*sqrt(pi)*sin(alph)^2*cos(alph)*exp(-Spar)*(besseli(0,Spar,scl)+besseli(1,Spar,scl))+(2*sqrt(pi)/S)*cos(alph)*exp(-Spar)*besseli(0,Spar,scl));
    CD      =   CN*sin(alph)+CA*cos(alph);
    CL      =   CA*sin(alph)+CN*cos(alph);
end

