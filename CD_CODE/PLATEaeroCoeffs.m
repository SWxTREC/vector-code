function [CDXYZtot,CDtot,Atot,Ftot,TQtot]=PLATEaeroCoeffs(TRI,V,n,m,Ti,Tw,accom,EPSIL,nu,phi_o,m_surface,ff,Rcm,set_acqs)
%given an arbitrary triangle file, computes the FMF drag of the entire
%surface, cross sectional area, and total force, assumes
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

%compute the cross sectional (Across) and planform (Ascale) areas and the
%anles of attack for each element
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
    
    %the Maxwellian component
    epsil               =   EPSIL(k);
    [CDm,~,CNm,CAm]     =   sentman(alph,Ti,Tw,accom,epsil,Vmag,m,'surfacespec',-1,1,1);
    
    %the quasi-specular component (Schamberg)
    CDqs                =   0;    
    CNqs                =   0;
    CAqs                =   0;
    if ff > 0%(pi/2-alph)>-theta_therm && 
        %CDqs         	=   schamberg_plate(pi/2-alph,nu,Vmag,Ti,m,m_surface(k));%m_surfase is an indexed vector
        [CDqs,CNqs,CAqs]   	=   schamberg_plate2(attk, nu(k),phi_o(k)*pi/180,Vmag,Ti,m,m_surface(k)*amu,Tw,0,set_acqs,accom);%m_surfase is an indexed vector
                               %schamberg_plate2(theta,nu,   phi_o,          Uinf,Tatm,m_gas,m_surface,        Tw,htrhmFlag,set_acqs,accomm) 
    end
    %CNqs                =   CDqs*cos(alph);
    %CAqs                =   CDqs*sin(alph);
    
    %combined coefficients    
    CN                  =   (1-ff)*CNm + ff*CNqs;
    CA                  =   (1-ff)*CAm + ff*CAqs;
    CD                  =   (1-ff)*CDm + ff*CDqs;
    
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