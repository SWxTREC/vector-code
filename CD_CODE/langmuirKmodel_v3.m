function K_model = langmuirKmodel_v3(params,x,Kf_1)
%params:    [Tatm, Eb_1, Ko_1]
J2eV                    =   6.24150974e18;%eV per Joule
%x is the input velocity
Eb_1                    =   params(2)/J2eV;
Tatm_1                  =   params(1);
Ko_1                    =   params(3);%3.3e6, 
% Kf_1                    =   params(4);%7.0e3, ;
% Kf_1                    =   3.0e4;%[7.0e3], 2.0e4, 6.5e4, 7e4;

if Tatm_1<0
    Tatm_1              =   eps;
end
if Eb_1<0
    Eb_1                =   eps;
end
if Kf_1<5.0e3
    Kf_1                =   5.0e3;
end
kb                      =   1.3806503e-23;
mO                      =   2.6560178e-26;                              %atomic oxygen mass (~16 amu) [kg]

xe                      =   0.5*mO*x.^2;%J energy input


test_section            =   exp(2*sqrt(Eb_1*xe)/(kb*Tatm_1));
k_range                 =   find(test_section==Inf,1);
if ~isempty(k_range)
    test_section        =   1.00e+306;
end
so_3                    =   (sqrt(pi*kb*Tatm_1*xe).*(erf((sqrt(Eb_1)-sqrt(xe))/sqrt(kb*Tatm_1))+erf(sqrt(xe/(kb*Tatm_1))))...
                            +kb*Tatm_1*exp(-(Eb_1+xe)/(kb*Tatm_1)).*(exp(Eb_1/(kb*Tatm_1))-test_section))...
                            ./(sqrt(pi*kb*Tatm_1.*xe).*(erf(sqrt(xe/(kb*Tatm_1)))+1) + kb*Tatm_1*exp(-xe/(kb*Tatm_1)));
K_model                 =   so_3*Ko_1+Kf_1;

knaninf                 =   find(isnan(K_model) | K_model==Inf);
if ~isempty(knaninf)
    K_model(knaninf)    =   Kf_1;
end