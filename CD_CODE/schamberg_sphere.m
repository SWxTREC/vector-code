function [CD]=schamberg_sphere(nu,phi_o,Uinf,Tatm,m_gas,m_surface,Tw,htrhmFlag,set_acqs,accomm)
%returns the drag coefficient for a sphere according to Schamberg's model with
%htrhmFlag - set to 1 for hyperthermal calculation (quicker but less accurate), 0 otherwise

%Boltzmann constant
kb              =   1.3806503e-23;

%beamwith function
tmp             =   2*phi_o/pi;
PHI             =    ((1-tmp^2)/(1-4*tmp^2))*(0.5*sin(2*phi_o)-tmp)/(sin(phi_o)-tmp);%phi=15deg ===> PHI=0.99
if isnan(PHI)
    phi_o       =   phi_o + 0.00001;
    tmp       	=   2*phi_o/pi;
    PHI      	=    ((1-tmp^2)/(1-4*tmp^2))*(0.5*sin(2*phi_o)-tmp)/(sin(phi_o)-tmp);
end

%mass ratio
mu              =   m_gas/m_surface;

%accommodation coefficient (Goodman, averaged over a half sphere)
if set_acqs==0
    accomm          =   abs(3.6*mu/((1+mu)^2));%abs(3.6*mu*sin(theta)/((1+mu)^2))
end

%rms thermal speed
%c               =   sqrt(2*kb*Tatm/m_gas);
c               =   sqrt(3*kb*Tatm/m_gas);

%inverse of the speed ratio (Schamberg's c/U)
s               =   c/Uinf;

%kinetic energy expressed in temperature
Tin             =   m_gas*Uinf^2/(3*kb);


CDinf           =   CDhyperthermal(PHI,Uinf,accomm,nu,Tw,Tin);
if htrhmFlag == 1           
    CD          =   CDinf;
end

if htrhmFlag == 0       %currently set to 0, which means that the nonhyperthermal approximation is being used
    CD          =   CDnonhyperthermal(s, CDinf);
end



%% INTERNAL FUNCTIONS

%definite integral
function [I1_eval] = I1(nu) %I1 function can be found in Schamberg 1959 report

if nu == 1
    I1_eval = 1/4;
elseif nu == 500
    I1_eval = 1/3;
end

%hyperthermal drag coefficient
function [CDinf]  =   CDhyperthermal(PHI,Vr,accomm,nu,Tw,Tin)
Vout            =   Vr*sqrt(1+accomm*(Tw/Tin-1));

f_of_nu         =   2*(I1(nu) - (1/(nu+3)));

CDinf           =   2*( 1 + PHI * (Vout/Vr) * f_of_nu );


%nonhyperthermal drag coefficient approximation
function [CDfmf] =  CDnonhyperthermal(s, CDinf)

CDfmf = ((2 + sqrt(1 + s^2))/3) * sqrt(1 + s^2) * CDinf;






