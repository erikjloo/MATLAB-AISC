function [phiMnx,phiMny] = AISC_F5(E,Fy,Sx,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,aw,ho,rt)
%Author: Erik Gisesen Loo
%Date: 2/21/15
%Description: Computes Flexural Strength per AISC Ch. F5

% For doubly symmetric I-beams only:
Sxc = Sx; hc = ho;
Myc = Fy*Sxc;

% Bending Strength Reduction Factor
Rpg = min(1- (aw/(1200+300*aw))*(hc/tw-5.7*sqrt(E/Fy)),1);
fprintf('The bending strength reduction factor Rpg is %1.2f \n',Rpg);
% Plastic Moment
x = 1;  Mnx(x) = Rpg*Myc;  
fprintf('Mnx(1) = %1.2f k-in.\n',Mnx(x))

% Lateral Torsional Buckling
x = x+1;
if Lb > Lr
    Mnx(x) = Rpg*Sxc*min(Cb*pi^2*E/((Lb/rt)^2),Fy);
    fprintf('The critical elastic LTB moment Cb*Mcr is %1.2f \n',Mnx(x))
elseif Lb > Lp
    Mnx(x) = min(Cb*(Fy-0.3*Fy((Lb-Lp)/(Lr-Lp))),Fy);
    fprintf('The critical inelastic LTB moment Cb*Mcr is %1.2f \n',Mnx(x))
else
    Mnx(x) = Mnx(1);
end

% Flange local buckling
if lambda_f > lambda_rf
    x = x+1;
	kc=4/sqrt(ho/tw); kc = min(max(kc,0.35),0.76);
    fprintf('kc = 4/sqrt(h_o/t_w) = %1.2f',kc); 
    Mnx(x)=0.9*E*kc*Sx/(lambda_f^2);
    fprintf('The slender flange buckling moment is %1.2f \n',Mnx(x))
elseif lambda_f > lambda_pf
    x = x+1;
    Mnx(x) = Rpg*Sxc*((Fy-0.3*Fy)*(lambda_f-lambda_pf)/(lambda_rf-lambda_pf));
    fprintf('The non-compact flange buckling moment is %1.2f \n',Mnx(x))
end

Mnx = min(Mnx);
fprintf('The nominal moment is %1.2f \n',Mnx)
phiMnx=0.9*Mnx;   fprintf('The design moment is %1.2f \n',phiMnx)