function Atotal = obliqueCylProjection(D,L,pitch)
%fids the oblique projection of a cylinder with angle pitch in radians
if pitch==0
    pitch   =   eps;
end
if pitch<0
    pitch   =   abs(pitch);
end

%ellipse parameters
a           =   D/2;
b           =   D*cos(pitch)/2;
Aellipse    =   a*b*pi;%single ellipse area


%sector angle
Phi         =   asin(L*sin(pitch)/(2*b));

%sector area
Asec        =   (a*b)*(atan((b/a)*tan(pi/2)) - atan((b/a)*tan(Phi)));

%triangle area
y           =   L*sin(pitch)/2;
base        =   2*y/tan(Phi);
Atri        =   0.5*base*y;

%segment area
Aseg        =   Asec-Atri;

%front cap oblique projection
Afront      =   Aellipse;

%back cap oblique projection
Aback       =   Aellipse;
if y<b
    Aback   =   Aback - 2*Aseg;
end

%find the projected midsection area - overlap
Amid        =   L*sin(pitch)*D;
if y>=b
    Amid    =   L*sin(pitch)*D - Aellipse;
end
if y<b
    Amid    =   L*sin(pitch)*D - (Aellipse - 2*Aseg);
end

Atotal      =   Afront + Amid + Aback;
