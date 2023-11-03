%TEST SCRIPT

%obj_type:  Object shape designated by 1-sphere, 2-plate with one side
%           exposed to flow, 3-cylinder, 4-geometry file
%EA_model:  Method of energy accommodation coefficient computation
%           -1:SESAM model, 0:set to constant value,  2: goodman model

%% Pitch test
obj_type=2; D=1; L=0; A=pi*D^2/4; Phi=30; Theta=0; Ta=10; Va=15000; n_O=1E14; n_O2=0; n_N2=0; n_He=0; n_H=0; EA_model=2; alpha=0.93; m_s=65; POSVEL=[]; fnamesurf="tmp.txt";
Pitch_ARR   = 0:2.5:360;
CDplt_ARR  	= zeros(length(Pitch_ARR),1);
FCplt_ARR  	= zeros(length(Pitch_ARR),1);

%plate
for k=1:length(Pitch_ARR)
    Phi         =   Pitch_ARR(k);
    
    [CD_status, CD, Aout, Fcoef, alpha_out] = MAIN(obj_type,D,L,A,Phi,Theta,Ta,Va,n_O,n_O2,n_N2,n_He,n_H,EA_model,alpha,m_s,POSVEL,fnamesurf);
    
    CDplt_ARR(k)=   CD;
    FCplt_ARR(k)=   Fcoef;
end

%0 width cylinder
obj_type=3;
CDcyl_ARR	= zeros(length(Pitch_ARR),1);
FCcyl_ARR  	= zeros(length(Pitch_ARR),1);

for k=1:length(Pitch_ARR)
    Phi         =   Pitch_ARR(k);
    
    [CD_status, CD, Aout, Fcoef, alpha_out] = MAIN(obj_type,D,L,A,Phi,Theta,Ta,Va,n_O,n_O2,n_N2,n_He,n_H,EA_model,alpha,m_s,POSVEL,fnamesurf);
    
    CDcyl_ARR(k)=   CD;
    FCcyl_ARR(k)=   Fcoef;
end

figure(); 
subplot(2,1,1)
plot(Pitch_ARR,CDplt_ARR,'ok','linewidth',2)
hold on
plot(Pitch_ARR,CDcyl_ARR,'xr','linewidth',2)
axis([0,max(Pitch_ARR)*1.01,0,4.0])
ylabel('C_{D}','FontSize',16)
set(gca,'FontSize',16)
set(gcf,'Color','white')
subplot(2,1,2)
plot(Pitch_ARR,FCplt_ARR,'ok','linewidth',2)
hold on
plot(Pitch_ARR,FCcyl_ARR,'xr','linewidth',2)
ylabel('C_{F}','FontSize',16)
xlabel('Pitch [deg]','FontSize',16)
legend('Plate w/ one Side Exposed','Zero-Width Cylinder (Plate w/ Both Sides Exposed)')
set(gca,'FontSize',16)
set(gcf,'Color','white')
axis([0,max(Pitch_ARR)*1.01,0,4.0])