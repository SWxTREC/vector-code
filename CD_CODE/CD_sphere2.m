function COEFS = CD_sphere2(Vt,NO_DENS,Tatm,Tw,accom,ff,nu,phi_o,ms,set_acqs)

%constants, 
mO                  =   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
mO2                 =   mO*2;
mN2                 =   4.6528299e-26;                              %molecular nitrogen mass [kg]
mHe                 =   6.6465e-027;
mH                  =   1.6737e-027;

MASS_MAT            =   [mH,mHe,mO,mN2,mO2];
RHO_MAT             =   MASS_MAT.*NO_DENS;
RhoTot              =   dot(MASS_MAT,NO_DENS);
CDpart              =   zeros(1,5);
dCDpart             =   zeros(1,5);


%scan atomic masses
for km=1:5
    if NO_DENS(km) == 0%skip zero number denisties
        continue
    end
    [CDm,~,~,~]     =   sentman(0,Tatm,Tw,accom,0,Vt,MASS_MAT(km),'sphere',-1,1,1);
    [CDqs,~,~]      =   schamberg_sphere(nu,phi_o*pi/180,Vt,Tatm,MASS_MAT(km),ms,Tw,0,set_acqs,accom);
    CDpart(1,km)    =   (1-ff)*CDm + ff*CDqs;
end

CDL                 =   dot(CDpart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient

COEFS               =   [0 CDL 0];%placeholders for future
