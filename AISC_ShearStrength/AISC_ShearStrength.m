    %% Check for Shear Strength
    kv = 5;
    if lambda_w >= 1.37*sqrt(kv*E/Fy);      Cv = 1.51*kv*E/((lambda_w^2)*Fy);
    elseif lambda_w >= 1.10*sqrt(kv*E/Fy);  Cv = 1.10*sqrt(kv*E/Fy)/lambda_w;
    else    Cv = 1; 
    end

    fprintf('Cv = %1.2f.\n',Cv)
    Vn = 0.6*Fy*A_w*Cv; fprintf('The nominal shear strength is %1.2f k.\n',Vn)
    phiVn=0.9*Vn;       fprintf('The phiVn is %1.2f k.\n',phiVn)
    close(f);