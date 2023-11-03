function [CDXYZtot,CDtot,Atot,Ftot,TQtot,alpha_out]=PLATEaeroCoeffs(TRI,V,n,m,Ti,Tw,accom,EPSIL,nu,phi_o,m_surface,ff,Rcm,set_acqs,material,GSI_model)
%given an arbitrary triangle file, computes the FMF drag of the entire
%surface, cross sectional area, and total force, for a single species assumes
%ff is the fraction of quasi specular component
%Rcm        =   cm offset from coordinate system
%
%example:
%%compute plate-model CD --------------TRSPROPS(k,:)=>[diffuse_fraction qs_fraction specular_fraction nu phi_o m_amu]
%[CDXYZtotP,CDtotP,ArefFMF,FfmfTot,TQ] = PLATEaeroCoeffs(TRS,Vblk,n,m,Tatm,Tsurf,accom,TRSPROPS(:,3),1,m_gls*ones(n_triangles,1),0);

%constants
amu                     =   1.6605e-27;
Vmag                    =   norm(V);
[ntri,~]                =   size(TRI);
kb                      =   1.3806503e-23;
theta_therm             =   atan(sqrt(2*kb*Ti/m)/Vmag);

%allocate storage vectors
ACROSS                  =   zeros(ntri,1);
AFLAG                   =   zeros(ntri,1);%into 1, or out of 0 flow
ASCALE                  =   zeros(ntri,1);
NVEC                    =   zeros(ntri,3);%normal vecor
AVEC                    =   zeros(ntri,3);%axial or along plate vector
CNVEC                   =   zeros(ntri,1);
CAVEC                   =   zeros(ntri,1);
CDXYZ                   =   zeros(ntri,3);
CDVEC                   =   zeros(ntri,1);
FXYZ                    =   zeros(ntri,3);
CENTROIDS               =   zeros(ntri,3);


cos_frac_vals = zeros(ntri,1);
accom_vals = zeros(ntri,1);
alpha_n_vals = zeros(ntri,1);
sigma_t_vals = zeros(ntri,1);


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


%compute the cross sectional (Across) and planform (Ascale) areas and the
%angles of attack for each element
%loop over triangles
for k=1:ntri
    [N,Ct,A]            =   getTRIprops(TRI(k,1:9));%triangle normal and center, and planform area
    %[PAMb,J_hat,K_hat]  =   SURFACE_GEOMETRY2(TRI(k,1:9),N,V/Vmag,0,0);
    [PAM,~,~]           =   surfacegeometry_mx(TRI(k,1:9),N,V/Vmag,0);%compute the projected area map of a single triangle
    CENTROIDS(k,1:3)    =   Ct+Rcm;
    dir_flag            =   sign(dot(V,N));%+1 facing away from flow, -1 facing into flow
    XP                  =   [PAM(1,1),PAM(1,3),PAM(1,5)];
    YP                  =   [PAM(1,2),PAM(1,4),PAM(1,6)];
    ACROSS(k,1)         =   polyarea(XP,YP);%projected area of triangle
    %ACROSS(k,1)         =   1.149;
    AFLAG(k,1)          =   0;
    if dir_flag==-1
        AFLAG(k,1)      =   1;
    end
    ASCALE(k,1)         =   A;
    
    %orientation
    cross_vec           =   cross(V,N);
    if norm(cross_vec)>eps
        cross_vec       =   cross_vec/norm(cross_vec);
    end
    if norm(cross_vec)<=eps
        V               =   V+[eps -eps eps];
        cross_vec       =   cross(V,N);%added /norm(cross_vec)
    end
    a_vec               =   cross(N,cross_vec);
    a_vec               =   a_vec/norm(a_vec);
    NVEC(k,1:3)         =   N(1:3)/norm(N);
    AVEC(k,1:3)         =   a_vec(1:3)/norm(a_vec);%vector along the surface 

    %compute the in plane and out of plane coefficients
    alph                =   acos(dot(V/Vmag,N/norm(N)));%?
    attk                =   -1*dir_flag*acos(dot(V/Vmag,a_vec));%angle of attack

    dotted_VN(k,1) = acos(dot(V/Vmag,N/norm(N)));


    if GSI_model == 3 %CLL quasi-specular reflection with alpha_n = 0.75, sigma_t = 0.9
        alpha_n = 0.75;
        sigma_t = 0.9;
        alpha_t = sigma_t*(2-sigma_t);
        [CDm,CLm,CNm,CAm] =      CLL_plate(alph,alpha_n,sigma_t,Vmag,Ti,m,Tw);
        CN = CNm;
        CA = CAm;
        CD = CDm;
        alpha_out = (alpha_n+alpha_t)/2;

    elseif GSI_model == 4 %Extrapolated laboratory-derived GSI parameters that have been weighted by surface area
        cos_frac = interp1(incident_angles_data, cos_frac_data, rad2deg(dotted_VN(k,1)), 'linear', 'extrap');
	    accom = interp1(incident_angles_data, alpha_data, rad2deg(dotted_VN(k,1)), 'linear', 'extrap');
	    alpha_n = interp1(incident_angles_data, alpha_n_data, rad2deg(dotted_VN(k,1)), 'linear', 'extrap');
	    sigma_t = interp1(incident_angles_data, sigma_t_data, rad2deg(dotted_VN(k,1)), 'linear', 'extrap');

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

        cos_frac_vals(k,1) = cos_frac;
        accom_vals(k,1) = accom;
        alpha_n_vals(k,1) = alpha_n;
        sigma_t_vals(k,1) = sigma_t;

        epsil               =   EPSIL(k);
        [CDm,~,CNm,CAm]     =   sentman(alph,Ti,Tw,accom,epsil,Vmag,m,'surfacespec',-1);
        [CDqs,CLqs,CNqs,CAqs] =      CLL_plate(alph,alpha_n,sigma_t,Vmag,Ti,m,Tw);
    
        %combined coefficients
        CN                  =   (1-cos_frac)*CNqs + cos_frac*CNm;
        CA                  =   (1-cos_frac)*CAqs + cos_frac*CAm;
        CD                  =   (1-cos_frac)*CDqs + cos_frac*CDm;

        alpha_cos = accom;
        alpha_t = sigma_t*(2-sigma_t);
        alpha_qs = (alpha_n+alpha_t)/2;
        alpha_out = cos_frac*alpha_cos + (1-cos_frac)*alpha_qs;

    else %Sentman diffuse with incomplete and/or variable energy accommodation
        epsil               =   EPSIL(k);
        [CDm,~,CNm,CAm]     =   sentman(alph,Ti,Tw,accom,epsil,Vmag,m,'surfacespec',-1);
        alpha_out = accom;

        CN = CNm;
        CA = CAm;
        CD = CDm;

    end
    
    %store coefficients
    CNVEC(k,1)          =   CN;
    CAVEC(k,1)          =   CA;

    %compute force coefficeints along the coordinate axes and along the
    %velocity vector
    CDXYZ(k,1)          =   dot(NVEC(k,1:3),[1 0 0])*CN+dot(AVEC(k,1:3),[1 0 0])*CA;
    CDXYZ(k,2)          =   dot(NVEC(k,1:3),[0 1 0])*CN+dot(AVEC(k,1:3),[0 1 0])*CA;
    CDXYZ(k,3)          =   dot(NVEC(k,1:3),[0 0 1])*CN+dot(AVEC(k,1:3),[0 0 1])*CA;
    CDVEC(k,1)          =   CD;
    
    
end


%compute the total cross sectional area, Atot
Atot                    =   sum(ACROSS.*AFLAG);

%compute the total coefficients
CDtot                   =   sum(CDVEC.*ASCALE)/Atot;
CDXYZtot(1,1)           =   sum(CDXYZ(1:end,1).*ASCALE)/Atot;
CDXYZtot(2,1)           =   sum(CDXYZ(1:end,2).*ASCALE)/Atot;
CDXYZtot(3,1)           =   sum(CDXYZ(1:end,3).*ASCALE)/Atot;

%compute element forces
FXYZ(:,1)               =   -0.5*CDXYZ(1:end,1).*ASCALE*Vmag^2*(n*m);
FXYZ(:,2)               =   -0.5*CDXYZ(1:end,2).*ASCALE*Vmag^2*(n*m);
FXYZ(:,3)               =   -0.5*CDXYZ(1:end,3).*ASCALE*Vmag^2*(n*m);
TORQUES                 =   cross(CENTROIDS,FXYZ,2);
TQtot                   =   sum(TORQUES,1);

%compute the total force
Ftot                    =   -0.5*Atot*CDXYZtot*Vmag^2*(n*m);





%%
%INTERNAL FUNCTION: COMPUTE TRIANGLE NORMAL, CENTER, AND PLANFORM AREA
function [N,Ct,AA] = getTRIprops(TRI)
%given a row of the TRI matrix computes its normal (and center pos)

%get veritces
A               =   [TRI(1:end,1),TRI(1:end,2),TRI(1:end,3)];
B               =   [TRI(1:end,4),TRI(1:end,5),TRI(1:end,6)];
C               =   [TRI(1:end,7),TRI(1:end,8),TRI(1:end,9)];

l1              =   B-A;
l2              =   C-A;

%right handed normal?
N               =   cross(l1(1:end,1:3),l2(1:end,1:3));

%normalize each entry
NRMSUM          =   N(1:end,1).^2+N(1:end,2).^2+N(1:end,3).^2;%normalize
NRMSUM          =   NRMSUM(1:end,1).^0.5;
N(1:end,1)      =   N(1:end,1)./NRMSUM;
N(1:end,2)      =   N(1:end,2)./NRMSUM;
N(1:end,3)      =   N(1:end,3)./NRMSUM;

%triangle centroid
Ct              =   (A(1:end,1:3)+B(1:end,1:3)+C(1:end,1:3))./3;

%triangle area
base            =   l1;
hdir            =   cross(N,base);
h               =   dot(l2,hdir)/norm(hdir);
AA              =   0.5*norm(base)*h;
