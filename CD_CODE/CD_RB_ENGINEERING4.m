function [accom,P_o,THETA] = CD_RB_ENGINEERING4(ro,mb,Ta,Vp,ms_amu,accom_mode,Xp,Eb,Kf,Ko)
%ACCOMMODATION ONLY
%coefficient for a tumbling rocket body. 
%this function is written to be called in a loop
%inputs:
%h              - altitude [km]
%ro             - mass density of atomic oxygen
%mb             - mean molecular mass (kg)
%NO_DENS        - number densities [mN2 mO2 mO mHe mH]  
%Vp             - perigee velocity magnitude (winds and corotating atmosphere excluded) [m/s]
%intmethod      - type of interpolation [string]
%accom_mode     - flag for computing accommodation -1:pilinski model, 
%                 0-1:set to constant value,  2: goodman model, 3:Afonso
%                 1985
%outputs:
%CD             - drag coefficient
%alpha          - accommodation coefficient
%P_o            - partial atomic oxygen pressure

%% constants
mO            	=   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
kb              =   1.3806503e-23;  %Boltzmann constant [J/K]
PaInTorr        =   (101325/760);
amu             =   1.6605e-027;
% Eb              =   5.707;%[5.707];%eV
Ta_trans        =   93.31;
%Ko              =   5E6;%2.8416E6<<<fit results, 



%% Compute Effective Pressure Parameter for Impinging Atomic Oxygen
s_o             =   Vp*sqrt(mO/(kb*Ta));
CD_O            =   ((2*s_o^2+1)/(sqrt(pi)*s_o^3))*exp(-s_o^2) + ((4*s_o^4+4*s_o^2-1)/(2*s_o^4))*erf(s_o);
P_o             =   Xp*0.5*ro*CD_O*Vp^2/PaInTorr;


%% Accommodation coefficient: alpha(P_o,Vp)
mu              =   mb/(ms_amu*amu);
accom_s         =   2.434*mu./((1+mu).^2);%3.6*mu./((1+mu).^2)
K               =   langmuirKmodel_v3([Ta_trans, Eb, Ko],Vp,Kf);%Vp

THETA           =   K.*P_o./(1+K.*P_o);
if accom_mode==2%implement the average goodman model
    THETA       =   zeros(size(THETA));
end

accom           =   THETA + (1-THETA).*accom_s;

if accom_mode>=0.0 && accom_mode<=1.0%implement fixed accommodation model
    accom       =   ones(size(accom))*accom_mode;
end


