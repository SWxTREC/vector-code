function COEFS = CD_plate_effective(alph,Vt,NO_DENS,MASS_MAT,Tatm,Tw,accom,ff,nu,phi_o,ms,set_acqs,A,material,GSI_model)

%constants, 
mO                  =   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]
mO2                 =   mO*2;
mN2                 =   4.6528299e-26;                              %molecular nitrogen mass [kg]
mHe                 =   6.6465e-027;
mH                  =   1.6737e-027;

RHO_MAT             =   MASS_MAT.*NO_DENS;
RhoTot              =   dot(MASS_MAT,NO_DENS);
CDpart              =   zeros(1,5);
CNpart              =   zeros(1,5);
CApart              =   zeros(1,5);

%cross sectional area
% if alph == pi/2%prevent zero-crossing singularity
%    alph= pi/2-eps; 
% end
Acrs                =   abs(A*cos(alph));


if material == 'SiO2'
    %range of SiO2 parameters
    incident_angles_data = 180 - [30, 45, 60]; %degrees
    alpha_n_data = [0.98, 0.99, 0.98];
    sigma_t_data = [0.49, 0.81, 0.83];
    cos_frac_data = [0.97, 0.9, 0.79];
    alpha_data = [0.71, 0.61, 0.44];

elseif material == 'aluminum'
    %range of aluminum parameters
    incident_angles_data = 180 - [30, 60]; %degrees
    alpha_n_data = [0.99, 0.98];
    sigma_t_data = [0.59, 0.85];
    cos_frac_data = [0.98, 0.88];
    alpha_data = [0.72, 0.55];

elseif material == 'Teflon'
    %range of teflon parameters
    incident_angles_data = 180 - [30, 45, 60]; %degrees
    alpha_n_data = [0.99, 0.99, 0.98];
    sigma_t_data = [0.69, 0.81, 0.84];
    cos_frac_data = [0.96, 0.87, 0.76];
    alpha_data = [0.62, 0.50, 0.30];

elseif material == 'FR4'
    %range of FR4 parameters
    incident_angles_data = 180 - [30, 45, 60]; %degrees
    alpha_n_data = [0.99, 0.99, 0.98];
    sigma_t_data = [0.59, 0.79, 0.85];
    cos_frac_data = [0.97, 0.93, 0.85];
    alpha_data = [0.68, 0.63, 0.52];

end


%scan atomic masses
for km=1:5
    if NO_DENS(km) == 0%skip zero number denisties
        continue
    end

    if GSI_model == 3 %CLL quasi-specular reflection with alpha_n = 0.75, sigma_t = 0.9
        alpha_n = 0.75;
        sigma_t = 0.9;
        alpha_t = sigma_t*(2-sigma_t);
        [CDm,CLm,CNm,CAm] =      CLL_plate(pi - alph,alpha_n,sigma_t,Vt,Tatm,MASS_MAT(km),Tw);
        CNpart(1,km) = CNm;
        CApart(1,km) = CAm;
        CDpart(1,km) = CDm;
        alpha_out = (alpha_n+alpha_t)/2;

    elseif GSI_model == 4 %Extrapolated laboratory-derived GSI parameters that have been weighted by surface area
        cos_frac = interp1(incident_angles_data, cos_frac_data, rad2deg(pi-alph), 'linear', 'extrap');
	    accom = interp1(incident_angles_data, alpha_data, rad2deg(pi-alph), 'linear', 'extrap');
	    alpha_n = interp1(incident_angles_data, alpha_n_data, rad2deg(pi-alph), 'linear', 'extrap');
	    sigma_t = interp1(incident_angles_data, sigma_t_data, rad2deg(pi-alph), 'linear', 'extrap');

        if cos_frac > 1
            cos_frac = 1;
        elseif cos_frac < 0
            cos_frac = 0;
        end
        if accom > 1
            accom = 1;
        elseif accom < 0
            accom = 0;
        end
        if alpha_n > 1
            alpha_n = 1;
        elseif alpha_n < 0
            alpha_n = 0;
        end
        if sigma_t > 1
            sigma_t = 1;
        elseif sigma_t < 0
            sigma_t = 0;
        end


        [CDm,~,CNm,CAm]     =   sentman(pi - alph,Tatm,Tw,accom,0,Vt,MASS_MAT(km),'surfacespec',-1);
        [CDqs,CLqs,CNqs,CAqs] =      CLL_plate(pi - alph,alpha_n,sigma_t,Vt,Tatm,MASS_MAT(km),Tw);
    
        %combined coefficients
        CNpart(1,km)                  =   (1-cos_frac)*CNqs + cos_frac*CNm;
        CApart(1,km)                  =   (1-cos_frac)*CAqs + cos_frac*CAm;
        CDpart(1,km)                  =   (1-cos_frac)*CDqs + cos_frac*CDm;

        alpha_cos = accom;
        alpha_t = sigma_t*(2-sigma_t);
        alpha_qs = (alpha_n+alpha_t)/2;
        alpha_out = cos_frac*alpha_cos + (1-cos_frac)*alpha_qs;

    else %Sentman diffuse with incomplete and/or variable energy accommodation
        [CDm,~,CNm,CAm]     =   sentman(pi - alph,Tatm,Tw,accom,0,Vt,MASS_MAT(km),'surfacespec',-1);
        alpha_out = accom;

        CNpart(1,km) = CNm;
        CApart(1,km) = CAm;
        CDpart(1,km) = CDm;

    end

end

CDL                 =   dot(CDpart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient
CAL                 =   dot(CApart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient
CNL                 =   dot(CNpart,RHO_MAT)/RhoTot;%recombine composite drag ceofficient

CDtot               =   CDL*A/Acrs;%scaled coefficients

COEFS               =   [CAL CDtot CNL Acrs alpha_out];%placeholders for future
