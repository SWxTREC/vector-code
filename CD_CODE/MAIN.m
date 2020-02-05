function [CD_status, CD, Aout, Fcoef, alpha_out] = MAIN(obj_type,D,L,A,Phi,Theta,Ta,Va,n_O,n_O2,n_N2,n_He,n_H,EA_model,alpha,m_s,POSVEL,fnamesurf)
%Computes CD at an array of points
%
%
%INPUTS:
%obj_type:  Object shape designated by 1-sphere, 2-plate with one side
%           exposed to flow, 3-cylinder, 4-geometry file
%
%D:
%
%L:
%
%EA_model:  Method of energy accommodation coefficient computation
%           -1:SESAM model, 0:set to constant value,  2: goodman model
%
%alpha:     energy accommodation coefficient
%
%m_s:       surface mass in amu
%
%POSVEL:    Nx8 array of time, position, and velocity inputs with the columns 
%           defined as[JD,Lon,Lat,Alt,Speed,F107A,F107,Ap]
%
%fnamesurf: surface element file in .wrl format
%OUTPUTS:   


%% Test Flags
%0-use matlab msis, 1-use matlab executable


%% Constants
set_acqs        =   0;%quasi-specular
ff              =   0;%specular fraction
nu              =   0;%quasi specular "bending" parameter
phi_o           =   0;%quasi specular lobe width
T_w             =   300;%wall temperature
mO           	=   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
mO2         	=   mO*2;
mN2             =   4.6528299e-26;                              %molecular nitrogen mass [kg]
mHe             =   6.6465e-027;
mH              =   1.6737e-027;
kb           	=   1.3806503e-23;  %Boltzmann constant [J/K]
Eb           	=   5.7;%eV
Kf              =   3e4;
Ko              =   5e6;



%% input check
%check size of position and velocity file
[rpv,cpv]       =   size(POSVEL);


%% multi point mode inputs
if rpv > 0 
    Npts        =   rpv;
    
    %retrieve atmospheric properties
    
    
    
end

%% single point mode inputs
if rpv == 0 
    Npts        =   1;
    %retrieve atmospheric properties
    NO_DENS     =   [n_N2 n_O2 n_O n_He n_H];
    V_rel       =   Va;
    T_atm       =   Ta;
    
    
end

%% Energy Accommodation Coefficients
%mean molecular mass calculation
MASS_MAT       	=   [mN2*ones(Npts,1) mO2*ones(Npts,1) mO*ones(Npts,1) mHe*ones(Npts,1) mH*ones(Npts,1)];
m_b             =   sum(MASS_MAT .* NO_DENS,2)./sum(NO_DENS,2);
ro              =   mO*NO_DENS(:,3);

%compute energy accommodation
EA_set          =   EA_model;
if EA_model == 0 
    EA_set      =   alpha;
end
[EA_vec,P_o,THETAsrf]  = CD_RB_ENGINEERING4(ro,m_b,T_atm,V_rel,m_s,EA_set,1.00,Eb,Kf,Ko);
alpha_out       =   EA_vec;

%% CD computation
%% Set Outputs
CD              =   -99*ones(Npts,1);
Aout            =   A*ones(Npts,1);
Fcoef           =   -99*ones(Npts,1);

for k=1:Npts
    %SPHERE
    if obj_type == 1
        
        COEFS       =   CD_sphere2(V_rel(k),NO_DENS(k,:),T_atm(k),T_w,EA_vec(k),0,1,0,m_s,0);
        CD(k)       =   COEFS(2);
        Aout(k)     =   A;
        Fcoef(k)    =   COEFS(2)*A;
        
    end
    
    %PLATE W/ ONE SIDE EXPOSED TO FLOW
    if obj_type == 2
        
        COEFS       =   CD_plate_effective(Phi*pi/180,V_rel(k),NO_DENS(k,:),T_atm(k),T_w,EA_vec(k),0,1,0,m_s,0,A);
        CD(k)       =   COEFS(2);
        Aout(k)     =   A*sin(Phi*pi/180);
        Fcoef(k)    =   COEFS(2)*Aout(k);%fix this later
        
    end
    
    %CYLINDER
    if obj_type == 3
        
        COEFS       =   CD_cyl_effective(V_rel(k),NO_DENS(k,:),T_atm(k),T_w,EA_vec(k),D,L,Phi*pi/180);
        CD(k)       =   COEFS(2);
        Aout(k)     =   obliqueCylProjection(D,L,Phi*pi/180);
        Fcoef(k)    =   COEFS(2)*Aout(k);%fix this later
        
    end
end

CD_status       =   1;
