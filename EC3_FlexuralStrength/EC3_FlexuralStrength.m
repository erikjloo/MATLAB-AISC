function [phiMnx] = EC3_FlexuralStrength(E,Fy,Zx,Iy,J,Cw,Lg,Lkip,C1,C2,kred)
%[phiMnx] = GiesenLoo_LTB(E,Fy,Zx,Iy,J,Cw,Lg,Lkip,C1,C2,kred)
%
% Author: Erik Gisesen Loo;   Date: 5/25/15
% Description: Computes & compares Flexural Strength

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if E == 210000; G = 81000; 
else G = E/2.6; 
end
S = sqrt(E*Cw/G/J); 

% Compute C
A = (pi*C1*Lg/Lkip); B = sqrt(1+((pi^2*S^2)/(Lkip^2)*(C2^2+1)));
D = (pi*C2*S/Lkip); C = A*(B+D);

% Find Mcr
Mcr = kred*C/Lg*sqrt(E*Iy*G*J); 
lambda = sqrt(Zx*Fy/Mcr);
a = 0.34;  b = 0.75; lambda_0 = 0.4;
phi = 0.5*(1+a*(lambda - lambda_0)+b*lambda^2);
Chi = min(1/(phi+sqrt(phi^2-lambda^2)),min(1/lambda^2,1));
phiMnx = Chi*Zx*Fy;

if Display
fprintf(' S = sqrt((E*I_w)/(G*I_t)) = %1.2f \n\n',S);
fprintf(' C = pi*C_1*L_g/L_kip*[sqrt(1+((pi^2*S^2)/(L_kip^2)*(C_2^2+1))) + (pi*C_2*S/L_kip)] = %1.3f \n\n',C); end
fprintf('M_cr = k_red*C/L_g*sqrt(E*I_y*G*J) = %1.2f \n\n', Mcr);
fprintf('lambda = sqrt(Z_x*F_y/M_cr) = %1.2f \n\n', lambda);
fprintf('phi = 0.5*(1+a*(lambda - lambda_0)+b*lambda^2) = %1.2f', phi);
fprintf('Chi = min(1/(phi+sqrt(phi^2-lambda^2)),min(1/lambda^2,1)) = %1.2f',Chi);
end