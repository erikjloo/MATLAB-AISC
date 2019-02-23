function [phiMnx,Mnx] = AISC_F2F3(E,G,Fy,Zx,Sx,Iy,J,Cw,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,ho)
%Author: Erik Gisesen Loo
%Date: 2/21/15
%Description: Computes Flexural Strength per AISC Ch. F2 and F3

%% Plastic Moment (compact/non-compact/slender flanges)
    x = 1;
    Mnx(x) = Fy*Zx; 
    fprintf('Mnx(1) = %1.2f \n',Mnx(x))
    
% Flange local buckling
if lambda_f > lambda_rf
	x = x+1;
    kc=4/sqrt(ho/tw); kc = min(max(kc,0.35),0.76);
    fprintf('kc = 4/sqrt(ho/t_w) = %1.2f',kc);
    Mnx(x)=0.9*E*kc*Sx/(lambda_f^2);
    fprintf('The slender flange buckling moment is %1.2f \n',Mnx(x))
elseif lambda_f > lambda_pf
	x = x+1;
    Mnx(x) =Mp-(Mp-0.7*Fy*Sx)*(lambda_f-lambda_pf)/(lambda_rf-lambda_pf);
    fprintf('The non-compact flange buckling moment is %1.2f \n',Mnx(x))
end
    
% Lateral Torsional Buckling
x = x+1;
if Lb > Lr
    Mnx(x)= min(Cb*pi/Lb*sqrt(E*Iy*G*J+(pi*E/Lb)^2*Iy*Cw),Mnx(1));
    fprintf('The critical elastic LTB moment Cb*Mcr is %1.2f \n',Mnx(x))
elseif Lb > Lp
    Mnx(x)= min(Cb*(Mp-(Mp-0.7*Fy*Sx)*((Lb-Lp)/(Lr-Lp))),Mnx(1));
    fprintf('The critical inelastic LTB moment Cb*Mcr is %1.2f \n',Mnx(x))
else
    Mnx(x) =Mp;
end

Mnx = min(Mnx);
fprintf('The nominal moment is %1.2f \n',Mnx)
phiMnx=0.9*Mnx;   fprintf('The design moment is %1.2f \n',phiMnx)