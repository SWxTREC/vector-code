function COEFS = CD_sphere2(Vt,NO_DENS,MASS_MAT,Tatm,Tw,accom,ms,set_acqs,htrhmFlag,r,material,GSI_model)

%constants, 
mO                  =   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
mO2                 =   mO*2;
mN2                 =   4.6528299e-26;                              %molecular nitrogen mass [kg]
mHe                 =   6.6465e-027;
mH                  =   1.6737e-027;

RHO_MAT             =   MASS_MAT.*NO_DENS;
RhoTot              =   dot(MASS_MAT,NO_DENS);
CDpart              =   zeros(1,5);
dCDpart             =   zeros(1,5);

%scan atomic masses
for km=1:5
    if NO_DENS(km) == 0%skip zero number denisties
        continue
    end

    if GSI_model == 3 %CLL quasi-specular reflection with alpha_n = 0.75, sigma_t = 0.9
        alpha_n = 0.75;
        sigma_t = 0.9;
        alpha_t = sigma_t*(2-sigma_t);
        [CDm]             =   CLL_sphere(alpha_n,sigma_t,Vt,Tatm,MASS_MAT(km),Tw);
        CDpart(1,km) = CDm;
        alpha_out = (alpha_n+alpha_t)/2;

    elseif GSI_model == 4 %Extrapolated laboratory-derived GSI parameters that have been weighted by surface area
        [cos_frac, accom, alpha_n, sigma_t] = sphere_incident_angle_weighting(r,material);
        [CDqs]             =   CLL_sphere(alpha_n,sigma_t,Vt,Tatm,MASS_MAT(km),Tw);
        [CDm,~,~,~]     =   sentman(0,Tatm,Tw,accom,0,Vt,MASS_MAT(km),'sphere',-1);
        CD_combined = (1-cos_frac)*CDqs + cos_frac*CDm;
        CDpart(1,km) = CD_combined;
        alpha_cos = accom;
        alpha_t = sigma_t*(2-sigma_t);
        alpha_qs = (alpha_n+alpha_t)/2;
        alpha_out = cos_frac*alpha_cos + (1-cos_frac)*alpha_qs;

    else %Sentman diffuse with incomplete and/or variable energy accommodation
        [CDm,~,~,~]     =   sentman(0,Tatm,Tw,accom,0,Vt,MASS_MAT(km),'sphere',-1);
        CDpart(1,km) = CDm;
        alpha_out = accom;
        
    end

end

CDL                 =   dot(CDpart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient

COEFS               =   [0 CDL 0 alpha_out];%placeholders for future







