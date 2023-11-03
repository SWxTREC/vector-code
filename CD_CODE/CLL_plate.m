function [CD, CL, CN, CA] = CLL_plate(alph,alpha_n,sigma_t,Uinf,Tatm,m_gas,Tw)
%returns the drag coefficient for a flat plate according to the CLL model
%one side exposed to the flow
%alph is the angle from normal (units of radians)
%htrhmFlag - set to 1 for hyperthermal calculation (quicker but less accurate), 0 otherwise

%Beta is the angle-of-attack, NOT angle from normal (units of radians)
%Beta = pi/2-alph;
Beta = alph - pi/2;
%rad2deg(alph)
%rad2deg(Beta)

attk_angle = -pi/2-alph; %transform to Birds coordinates (angle of attack)

%universal gas constant
R = 8.3145; %kg m^2 s^-2 mol^-1 K^-1

%Avogadro's constant
N_A = 6.0221409e23; %mol^-1

%molecular weight (m_gas in units of kg)
m_molecular = m_gas * N_A; %kg/mol

%average normal velocity
V_w = sqrt(pi*R*Tw/(2*m_molecular)); %m/s

V = Uinf; %m/s

%speed ratio
s = V / sqrt(2*R*Tatm/m_molecular); %unitless

%from Walker et al. 2014 paper, partial relationship between alpha_n and sigma_n
sigma_n = 1 - sqrt(1 - alpha_n);

%parallel to the flow (i, x direction)
CD = (2/s) * (sigma_t*gamma_1(s*sin(Beta)) + ((2-sigma_n)/s)*gamma_2(s*sin(Beta))*sin(Beta) - sigma_t*gamma_1(s*sin(Beta))*(sin(Beta))^2 + sigma_n*(V_w/V)*gamma_1(s*sin(Beta))*sin(Beta));

%perpendicular to the flow (j, y direction)
CL = (2/s) * (((2-sigma_n)/s)*gamma_2(s*sin(Beta))*cos(Beta) - sigma_t*gamma_1(s*sin(Beta))*sin(Beta)*cos(Beta) + sigma_n*(V_w/V)*gamma_1(s*sin(Beta))*cos(Beta));

%CN denotes the force normal to the plate
%CN = 0;
%CN = abs( CD / cos(pi/2-Beta) )*-1
%CN = CL / (sin(pi/2 - Beta))
CN = abs(CD*sin(attk_angle) + CL*cos(attk_angle)) * -1;

%CA is the force parallel (along the plate)
%CA = 0;
CA = abs(CD*cos(attk_angle) - CL*sin(attk_angle));



%% INTERNAL FUNCTIONS

%gamma functions
function [g1] = gamma_1(x)
g1 = 1/(2*sqrt(pi)) * (exp(-x^2) + sqrt(pi)*x*(1+erf(x)));

function [g2] = gamma_2(x)
g2 = 1/(2*sqrt(pi)) * (x*exp(-x^2) + (sqrt(pi)/2)*(1+2*x^2)*(1+erf(x)));



