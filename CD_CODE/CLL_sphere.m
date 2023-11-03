function [CD] = CLL_sphere(alpha_n,sigma_t,Uinf,Tatm,m_gas,Tw)
%Tw is the surface temperature
%m_gas is in kg
%returns the drag coefficient for a sphere according to the CLL model

%universal gas constant
R = 8.3145; %kg m^2 s^-2 mol^-1 K^-1

%Avogadro's constant
N_A = 6.0221409e23; %mol^-1

%molecular weight (m_gas in units of kg)
m_molecular = m_gas * N_A; %kg/mol

%average normal velocity (from Storch 2002 paper)
V_w = sqrt(pi*R*Tw/(2*m_molecular)); %m/s

V = Uinf; %m/s

%speed ratio
s = V / sqrt(2*R*Tatm/m_molecular); %unitless

%from Walker et al. 2014 paper, partial relationship between alpha_n and sigma_n
sigma_n = 1 - sqrt(1 - alpha_n);

%from Storch 2002 paper
CD = (2-sigma_n+sigma_t)/(2*s^3) * (((4*s^4+4*s^2-1)/(2*s))*erf(s) + ((2*s^2+1)/sqrt(pi))*exp(-s^2)) + (4/3)*sigma_n*(V_w/V);


