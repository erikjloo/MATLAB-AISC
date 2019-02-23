%% Cantilever Beam
bf = 135; tf = 10.2; r = 15;
h = 270; tw = 6.6;
d = h-2*tf-2*r;

E = 210000;
Fy = 355;
Zx = 484000;

Iy = 419.9e4;
J = 1.15*(2*bf*tf^3/3+(h-tf)*tw^3/3)
Cw = Iy*(h-tf)^2/4

Lb = 3390;
Lr = 10000;
Lg = 3390;
Lst = Lg;
Lkip = Lg;
C = 1;
C1 = 1.45;
C2 = -0.56;
kred = 1;

[phiMnx] = GiesenLoo_LTB(E,Fy,Zx,Iy,J,Cw,Lg,Lkip,C1,C2,kred)


[phiMnx,phiMny] = GiesenLoo_FlexuralStrength(E,Fy,Zx,Sx,Zy,Sy,Iy,J,Cw,Lb,Lp,Lr,Cb,...
    type,Display,AISC16,Tension,lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,d,Sxc)
%% Lateral Stability of a Continuous Beam
clc
bf = 180; tf = 13.5; r = 18;
h = 400; tw = 8.6;
d = h-2*tf-2*r;

E = 210000;
Fy = 275;
Zx = 1308000;

Ix = 2313e5;
Iy = 1318e4;
J = 1.15*(2*bf*tf^3/3+(h-tf)*tw^3/3)
Cw = Iy*(h-tf)^2/4

Lb = 10000;
Lr = 10000;
Lg = 10000;
Lst = Lg;
Lkip = Lg;
C = 1.7;
C1 = 1.7;
C2 = -0.9;
kred = 1;

[phiMnx] = GiesenLoo_LTB(E,Fy,Zx,Iy,J,Cw,Lg,Lkip,C1,C2,kred)
[phiMnx,phiMny] = GiesenLoo_FlexuralStrength(E,Fy,Zx,Sx,Zy,Sy,Iy,J,Cw,Lb,Lp,Lr,Cb,...
    type,Display,AISC16,Tension,lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,d,Sxc)