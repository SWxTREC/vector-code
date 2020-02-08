function [CD_tot,CL_tot,ASC_cyl] = CD_cyl(D,L,accomm,pitch,m,V_inf,T_inf,T_surf)
%computes the drag coefficient of a capped cylinder (pitch in radians)

ASC_cyl         =   abs(obliqueCylProjection(D,L,pitch));%cross section of cylinder
Aref_surf       =   pi*D^2/4;
[CD,CL,CN,CA]   =   sentman(pitch,T_inf,T_surf,accomm,0,V_inf,m,'cylinder',-1);%,D,L);
CD_sides        =   CD;
CL_sides        =   CL;

[CD,CL,CN,CA]   =   sentman(pitch,T_inf,T_surf,accomm,0,V_inf,m,'surfacespec',-1);%,D,L);
CD_capF         =   CD;
CL_capF         =   CL;

[CD,CL,CN,CA]   =   sentman(pitch+pi,T_inf,T_surf,accomm,0,V_inf,m,'surfacespec',-1);%,D,L);
CD_capB         =   CD;
CL_capB         =   CL;

CD_tot          =   (CD_capF*Aref_surf + CD_capB*Aref_surf + CD_sides*L*D/2)/ASC_cyl;
CL_tot          =   (CL_capF*Aref_surf + CL_capB*Aref_surf + CL_sides*L*D/2)/ASC_cyl;

