function [CD_status, CD, Aout, Fcoef, alpha_out] = MAIN(obj_type,D,L,A,Phi,Theta,Ta,Va,n_O,n_O2,n_N2,n_He,n_H,EA_model,alpha,m_s,POSVEL,fnamesurf)
%
%
%Computes the free-molecular flow drag coefficients, force coefficients, 
%and areas of basic shapes as well as more complex geometries as specified 
%by a "plate model" geometry. Gas surface interaction parameters are 
%specified by either the SESAM model, a fixed value defined by the user, or 
%the Goodman model. Values can be computed at a single set of inputs or an
%array of inputs specified by position and velocity.
%
%
%Written by:Marcin D. Pilinski, February 2020
%           Space Weather Technology, Research and Education Center
%           Laboratory for Atmospheric and Space Physics
%           University of Colorado Boulder
%
%=========================================================================
%INPUTS
%obj_type:  Object shape designated by 1-sphere, 2-plate with one side
%           exposed to flow, 3-cylinder, 4-geometry file
%
%D:         Diameter for sphere and cylinder [m]
%
%L:         Length for cylinder [m]
%
%Phi:       Pitch angle [deg]
%
%Theta:     Sideslip anlge [deg]
%
%Ta:        Ambient atmospheric temperature [K]
%
%Va:        Free-stream velocity magnitude, or speed [m/s]
%
%n_O,n_O2,n_N2,n_He,n_H: partial number densities of atmospheric species 
%           atomic oxygen, molecular oxygen, molecular nitrogen, helium,
%           and hydrogen respectively. [m^-3]
%           
%EA_model:  Method of energy accommodation coefficient computation
%           -1:SESAM model, 0:set to constant value,  2: goodman model
%
%alpha:     energy accommodation coefficient
%
%m_s:       surface mass in amu
%
%POSVEL:    Nx8 array of time, position, and velocity inputs with the columns 
%           defined as[JD,Lon,Lat,Alt,Speed,F107A?,F107?,Ap?]
%
%fnamesurf: surface element file in .wrl format
%
%=========================================================================
%OUTPUTS
%CD_status: 1-run successful
%
%CD:        Drag coefficient
%
%Aout:      Cross Sectional Area [m^2]
%
%Fcoef:     Force Coefficient [m^2]
%
%alpha_out: Energy Accomodation Coefficient  


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

Rcm         	=   [-0.01 0.02 0.01]';           %position of center of mass from center of rectangular solid or sphere [m]


%plotting constants
ppscl         	=   1.25;
offsetx     	=   -0.05;
offsety        	=   -0.05;
ppwidth        	=   8.5*ppscl;
ppheight     	=   8.5*ppscl;


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
% Set Outputs
CD              =   -99*ones(Npts,1);
Aout            =   A*ones(Npts,1);
Fcoef           =   -99*ones(Npts,1);

%SURFACE MODEL (PLATE MODEL) GEOMETRY
if obj_type == 4 %geometry file
    
    %open geometry and convert to vertex array
    hf0       	=   figure('visible','off');
    ViewDir     =   [Phi,Theta];%view direction
    geometry_wizard(fnamesurf,'TRIprops.txt',0,1,ViewDir)
    set(gcf,'Color','white')
    set(gca,'FontSize',16)
    set(hf0, 'PaperSize', [ppwidth ppheight], 'PaperPosition', [0+offsetx 0+offsety ppwidth+offsetx ppheight+offsety])
    %print(hf0, '-dpdf', 'geometry.pdf')
    saveas(hf0,'geometry.png')

    %load triangle array
    TRS     	=   load('TRIprops.txt');
    n_triangles	=   length(TRS(:,1)); 
    
end


for k=1:Npts
    %SPHERE
    if obj_type == 1
        
        COEFS       =   CD_sphere2(V_rel(k),NO_DENS(k,:),T_atm(k),T_w,EA_vec(k),0,1,0,m_s,0);
        CD(k)       =   COEFS(2);
        Aout(k)     =   pi*D^2/4;
        Fcoef(k)    =   COEFS(2)*Aout(k);%fix this later
        
    end
    
    %PLATE W/ ONE SIDE EXPOSED TO FLOW
    if obj_type == 2
        
        COEFS       =   CD_plate_effective(Phi*pi/180,V_rel(k),NO_DENS(k,:),T_atm(k),T_w,EA_vec(k),0,1,0,m_s,0,A);
        CD(k)       =   COEFS(2);
        Aout(k)     =   COEFS(4);%abs(A*sin(Phi*pi/180));
        Fcoef(k)    =   COEFS(2)*Aout(k);%fix this later
        
    end
    
    %CYLINDER
    if obj_type == 3
        
        COEFS       =   CD_cyl_effective(V_rel(k),NO_DENS(k,:),T_atm(k),T_w,EA_vec(k),D,L,Phi*pi/180);
        CD(k)       =   COEFS(2);
        Aout(k)     =   COEFS(4);%abs(obliqueCylProjection(D,L,Phi*pi/180));
        Fcoef(k)    =   COEFS(2)*Aout(k);%fix this later
        
    end
    
     %SURFACE MODEL (PLATE MODEL) GEOMETRY
    if obj_type == 4
        
        %velocity vector
     	V_horz     	=   V_rel(k)*cos(Phi(k)*pi/180);
        V_z        	=   V_rel(k)*sin(Phi(k)*pi/180);
        V_x        	=   V_horz*cos(Theta(k)*pi/180);
        V_y       	=   V_horz*sin(Theta(k)*pi/180);
        V_plate_in	=   [V_x,V_y,V_z];
        
        
        EPSILprops	=   zeros(n_triangles,1);%specular fraction
        
        %CD computation
        tic
        COEFS     	=   CD_triFile_effective(TRS,V_plate_in,NO_DENS(k,:),MASS_MAT(k,:),T_atm(k),T_w,EA_vec(k),...
                                             EPSILprops, nu*ones(n_triangles,1), phi_o*ones(n_triangles,1),...
                                             m_s*ones(n_triangles,1),ff,Rcm',set_acqs);
        toc
        CD(k)       =   COEFS(2);
        Aout(k)     =   COEFS(4);%abs(obliqueCylProjection(D,L,Phi*pi/180));
        Fcoef(k)    =   COEFS(2)*Aout(k);%fix this later                                

    end
end

CD_status       =   1;
