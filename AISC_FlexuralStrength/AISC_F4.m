function [phiMnx,Mnx] = AISC_F4(E,Fy,Zx,Sx,J,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,ho,rt)
%Author: Erik Gisesen Loo
%Date: 2/13/17
%Description: Computes Flexural Strength per AISC Ch. F4

%% For doubly symmetric I-beams only:
Sxc = Sx;  hc = ho; 
Mp = min(Fy*Zx,1.6*Fy*Sx);     Myc = Fy*Sxc;

%% Web plastification factor
if lambda_w > lambda_pw %& lambda_w < lambda_rw
    Rpc = Mp/Myc - (Mp/Myc-1)*(lambda_w-lambda_pw)/(lambda_rw-lambda_pw);
else %lambda_w < lambda_pw
    Rpc = Mp/Myc;
end

%% Plastic Moment
x = 1;  Mnx(x) = Rpc*Myc;  
fprintf('Mnx(1) = %1.2f k-in.\n',Mnx(x))
 
%% Lateral Torsional Buckling
x = x+1;
if Lb > Lr;
    Fcr = (Cb*pi^2*E)/((Lb/rt)^2)*sqrt(1+0.078*J/(Sxc*ho)*(Lb/rt)^2);
	Mnx(x)= min(Fcr*Sxc,Rpc*Myc);
    fprintf('The critical elastic LTB moment Cb*Mcr is %1.2f \n',Mnx(x))
elseif Lb > Lp
    Mnx(x)= min(Cb*(Rpc*Myc-(Rpc*Myc-0.7*Fy*Sxc)*((Lb-Lp)/(Lr-Lp))),Rpc*Myc);
    fprintf('The critical inelastic LTB moment Cb*Mcr is %1.2f \n',Mnx(x))
else
    Mnx(x) = Mnx(1);
end

%% Flange local buckling
if lambda_f > lambda_rf
    x = x+1;
	kc=4/sqrt(ho/tw); kc = min(max(kc,0.35),0.76);
    fprintf('kc = 4/sqrt(h_o/t_w) = %1.2f',kc); 
    Mnx(x)=0.9*E*kc*Sx/(lambda_f^2);
    fprintf('The slender flange buckling moment is %1.2f \n',Mnx(x))
elseif lambda_f > lambda_pf
    x = x+1;
    Mnx(x) = Rpc*Myc-(Rpc*Myc-0.7*Fy*Sxc)*((lambda_f-lambda_pf)/(lambda_rf-lambda_pf));
    fprintf('The non-compact flange buckling moment is %1.2f \n',Mnx(x))
end

%%
Mnx = min(Mnx);
fprintf('The nominal moment is %1.2f \n',Mnx)
phiMnx=0.9*Mnx;   fprintf('The design moment is %1.2f \n',phiMnx)