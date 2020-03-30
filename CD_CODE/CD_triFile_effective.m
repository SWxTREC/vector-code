function COEFS = CD_triFile_effective(TRS,V_plate_in,NO_DENS,MASS_MAT,T_atm,T_w,accom,EPSILprops,NU,PHI_O,M_SURF,ff,Rcm,set_acqs)
%given an arbitrary triangle file, computes the FMF drag of the entire
%surface, cross sectional area, and total force, weighed over multiple species
%assumes ff is the fraction of quasi specular component
%Rcm        =   cm offset from coordinate system

%constants, 
mO                  =   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
mO2                 =   mO*2;
mN2                 =   4.6528299e-26;                              %molecular nitrogen mass [kg]
mHe                 =   6.6465e-027;
mH                  =   1.6737e-027;

RHO_MAT             =   MASS_MAT.*NO_DENS;
RhoTot              =   dot(MASS_MAT,NO_DENS);
CDpart              =   zeros(1,5);
CXpart              =   zeros(1,5);
CYpart              =   zeros(1,5);
CZpart              =   zeros(1,5);
FXpart              =   zeros(1,5);
FYpart              =   zeros(1,5);
FZpart              =   zeros(1,5);
TQXpart             =   zeros(1,5);
TQYpart             =   zeros(1,5);
TQZpart             =   zeros(1,5);


%contants


%cross sectional area
% if alph == pi/2%prevent zero-crossing singularity
%    alph= pi/2-eps; 
% end

%scan atomic masses
for km=1:5
    if NO_DENS(km) == 0%skip zero number denisties
        continue
    end
    
 	[CDXYZtot,CDtot,Atot,Ftot,TQtot]    =   PLATEaeroCoeffs(TRS, V_plate_in, NO_DENS(km), MASS_MAT(km), T_atm, T_w, accom, EPSILprops, NU, PHI_O,M_SURF,ff,Rcm,set_acqs);

    %[CDqs,~,~]      = schamberg_sphere(nu,phi_o*pi/180,Vt,Tatm,MASS_MAT(km),ms,Tw,0,set_acqs,accom);%<<<future capability
    CDpart(1,km)    =   CDtot;
    CXpart(1,km)    =   CDXYZtot(1);
    CYpart(1,km)    =   CDXYZtot(2);
    CZpart(1,km)    =   CDXYZtot(3);
    FXpart(1,km)   	=   Ftot(1);
    FYpart(1,km)   	=   Ftot(2);
    FZpart(1,km)   	=   Ftot(3);
    TQXpart(1,km)  	=   TQtot(1);
    TQYpart(1,km)  	=   TQtot(2);
    TQZpart(1,km) 	=   TQtot(3);
end

CDL                 =   dot(CDpart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient
CXL                 =   dot(CXpart,RHO_MAT)/RhoTot;%recombine composite CAx ceofficient
CYL                 =   dot(CYpart,RHO_MAT)/RhoTot;%recombine composite CAy ceofficient
CZL                 =   dot(CZpart,RHO_MAT)/RhoTot;%recombine composite CAz ceofficient
FXL                 =   dot(FXpart,RHO_MAT)/RhoTot;%recombine composite CAx ceofficient
FYL                 =   dot(FYpart,RHO_MAT)/RhoTot;%recombine composite CAy ceofficient
FZL                 =   dot(FZpart,RHO_MAT)/RhoTot;%recombine composite CAz ceofficient
TQXL                =   dot(TQXpart,RHO_MAT)/RhoTot;%recombine composite CAx ceofficient
TQYL             	=   dot(TQYpart,RHO_MAT)/RhoTot;%recombine composite CAy ceofficient
TQZL                =   dot(TQZpart,RHO_MAT)/RhoTot;%recombine composite CAz ceofficient

CL1=0;
CL2=0;

COEFS               =   [CL1 CDL CL2 Atot CXL CYL CZL FXL FYL FZL TQXL TQYL TQZL];%placeholders for future