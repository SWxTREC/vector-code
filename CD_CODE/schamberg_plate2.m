function [CD,CNRM,CTAN] = schamberg_plate2(theta,nu,phi_o,Uinf,Tatm,m_gas,m_surface,Tw,htrhmFlag,set_acqs,accomm)
%returns the drag coefficient for a flat plate according to Schambergs model with
%accommodation comupted using Goodman's formula, one side exposed ot flow
%only
%appropriate for q-s reflection only!
%theta is the angle-of-attack, NOT angle from normal
%htrhmFlag - set to 1 for hyperthermal calculation (quicker but less accurate), 0 otherwise

%Boltzmann constant
kb              =   1.3806503e-23;

%angle of reflection measured from surface plane
theta_out       =   abs(acos(cos(theta)^nu));
if isnan(theta_out)
    theta_out   =   pi/2;
end

%check for reflections back into the surface
if theta_out < phi_o
    phi_o       =   theta_out;
end

%beamwidth function
tmp             =   2*phi_o/pi;
PHI             =    ((1-tmp^2)/(1-4*tmp^2))*(0.5*sin(2*phi_o)-tmp)/(sin(phi_o)-tmp);%phi=15deg ===> PHI=0.99
if isnan(PHI)
    phi_o       =   phi_o + 0.00001;
    tmp       	=   2*phi_o/pi;
    PHI      	=    ((1-tmp^2)/(1-4*tmp^2))*(0.5*sin(2*phi_o)-tmp)/(sin(phi_o)-tmp);
end
%PHI

%mass ratio
mu              =   m_gas/m_surface;

%accommodation coefficient (Goodman)
if set_acqs==0
    accomm          =   abs(3.6*mu*sin(theta)/((1+mu)^2));%abs(3.6*mu*sin(theta)/((1+mu)^2))
end
    
%rms thermal speed
%c               =   sqrt(2*kb*Tatm/m_gas);
c               =   sqrt(3*kb*Tatm/m_gas); %this one is consistent with schamberg's report

%inverse of the speed ratio (Schamberg's c/U)
s               =   c/Uinf;

%kinetic energy expressed in temperature
Tin             =   m_gas*Uinf^2/(3*kb);

%theta prime
thetaP          =   (asin(sin((theta))/sqrt(1+s^2)));%<<<<<ABS?????<<<<<

%delta correction angle
delta           =   atan(s);
[CD1,CN1,CT1]  	=   CDhyperthermal(PHI,theta,Uinf,accomm,nu,Tw,Tin);
if htrhmFlag == 0       %currently set to 0, which means that the nonhyperthermal approximation is being used
    [CD2,CN2,CT2]  	=   CDhyperthermal(PHI,theta+delta,Uinf,accomm,nu,Tw,Tin);
    [CD3,CN3,CT3]  	=   CDhyperthermal(PHI,theta-delta,Uinf,accomm,nu,Tw,Tin);
    [CD4,CN4,CT4]  	=   CDhyperthermal(PHI,thetaP,Uinf,accomm,nu,Tw,Tin);


    %drag coefficient (Joule Gas Approximation)
    CD              =   (1/3)*sqrt(1+s^2)*( sqrt(1+s^2)*CD1 +...
                        0.5*(CD2 +...
                        CD3) +...
                        CD4);
                
    CNRM            =   (1/3)*sqrt(1+s^2)*( sqrt(1+s^2)*CN1 +...
                        0.5*(CN2 +...
                        CN3) +...
                        CN4);    
                
    CTAN            =   (1/3)*sqrt(1+s^2)*( sqrt(1+s^2)*CT1 +...
                        0.5*(CT2 +...
                        CT3) +...
                        CT4);    
end             
if htrhmFlag == 1           
    CD          =   CD1;
    CNRM        =   CN1;
    CTAN        =   CT1;
end


%% INTERNAL FUNCTIONS

%hyperthermal drag coefficient
function [CDinf,CNinf,CTinf]  =   CDhyperthermal(PHI,theta,Vr,accomm,nu,Tw,Tin)
Vout            =   Vr*sqrt(1+accomm*(Tw/Tin-1));
%theta           =   abs(theta);%<<<<<<<<<<<<<<<<

%because this is hyperthermal, negative incoming angles should be neglected
if theta < 0
    CDinf       =   0;
    CNinf       =   0;
    CTinf       =   0;
    return
end

%angle of reflection measured from surface plane
theta_out       =   abs(acos(cos(theta)^nu));
if isnan(theta_out)
    theta_out   =   pi/2;
end

if sin(theta) == 0 || cos(theta) == 0
    theta       =   theta + eps;
end 

CNinf        	=   -2*abs(sin(theta)^2*(1+PHI*(Vout/Vr)*(sin(theta_out)/sin(theta))));%always into surface
CTinf           =   abs(2*sin(theta)*cos(theta)*(1-PHI*(Vout/Vr)*(cos(theta_out)/cos(theta))));%always along flow
CDinf           =   CTinf*cos(theta) - CNinf*sin(theta);

