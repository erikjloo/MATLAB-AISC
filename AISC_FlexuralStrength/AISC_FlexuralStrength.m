function [phiMnx,phiMny,Mnx,Mny] = AISC_FlexuralStrength(E,Fy,Zx,Sx,Zy,Sy,Iy,J,Cw,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,ho,d)
%Author: Erik Gisesen Loo
%Date: 2/21/15
%Description: Computes Flexural Strength per AISC Ch. F
if E == 210000; G = 81000; else G = E/2.6; end


%% Major Axis Flexural Capacity
if lambda_w < lambda_pw
    [phiMnx,Mnx] = AISC_F2F3(E,G,Fy,Zx,Sx,Iy,J,Cw,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,ho);
elseif lambda_w < lambda_rw
    [phiMnx,Mnx] = AISC_F4(E,Fy,Zx,Sx,J,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,ho,rt);
else % lambda_w > lambda_rw
    AISC_F5
end

%% Minor Axis Flexural Capacity
Mny(1) = min(Fy*Zy,1.6*Fy*Sy);
fprintf('Mny(1) = %1.2f k-in.\n',Mny(1))
    
 % Check flange local buckling
 if lambda_f > lambda_rf
    fprintf('%1.2f > %1.2f --> Slender flange. \n',lambda_f,lambda_rf)
    Fcr = 0.69*E/((2*lambda_f)^2);
    Mny(2) = Fcr*Sy; fprintf('Mny = 0.69*E/((2*lambda_f)^2)*Sy = %1.2 ksi.\n',Mny(2))
elseif lambda_f > lambda_pf
	fprintf('%1.2f > %1.2f --> Non-compact flange. \n',lambda_f,lambda_pf)
	Mny(2)=Mpy-(Mpy-0.7*Fy*Sy)*(lambda_f-lambda_pf)/(lambda_rf-lambda_pf);
    fprintf('Mny = Mpy-(Mpy-0.7*Fy*Sy)*(lam_f-lam_pf)/(lam_rf-lam_pf) = %1.2f k-ft.\n',Mny(2))
end
phiMny = 0.9*min(Mny);   fprintf('Mny = %1.2f k-in. -->  phiMny = %1.2f k-in.\n\n',Mny,phiMny)
end