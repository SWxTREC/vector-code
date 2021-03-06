function [COEFS,dCDL] = CD_sphere(Vt,NO_DENS,Tatm,Tw,accom)

%constants, 
mO                  =   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
mO2                 =   mO*2;
mN2                 =   4.6528299e-26;                              %molecular nitrogen mass [kg]
mHe                 =   6.6465e-027;
mH                  =   1.6737e-027;

MASS_MAT            =   [mN2 mO2 mO mHe mH];
RHO_MAT             =   MASS_MAT.*NO_DENS;
RhoTot              =   dot(MASS_MAT,NO_DENS);
CDpart              =   zeros(1,5);
dCDpart             =   zeros(1,5);


%scan atomic masses
for km=1:5
    if NO_DENS(km) == 0%skip zero number denisties
        continue
    end
    [CD,~,~,~]      =   sentman(0,Tatm,Tw,accom,0,Vt,MASS_MAT(km),'sphere',-1);
    CDpart(1,km)    =   CD;
    dCDpart(1,km)   =   derviCDsphere(Tatm,Tw,accom,Vt,MASS_MAT(km));%partial derivative
end

CDL                 =   dot(CDpart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient
dCDL                =   dot(dCDpart,RHO_MAT)/RhoTot;

COEFS               =   [0 CDL 0];%placeholders for future
