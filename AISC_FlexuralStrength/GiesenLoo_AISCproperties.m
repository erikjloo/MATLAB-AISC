% GiesenLoo_AISCproperties
if ~exist('slabel');  slabel = input('Section name, e.g. W12X26 (use apostrophes around name):   '); end
if ~exist('E'); E = input('E (ksi) = '); end
if ~exist('Fy'); Fy = input('Fy (ksi) = '); end
if ~exist('AISClabel'); load('AISCShapes'); end
if ~exist('Display'); Display = 1; end

if all(slabel([1 2]) == 'WT') % strcmp(slabel,'WT')
    % Find the number of the property
	snum = find(strcmp(slabel,AISClabel),1);    
    % Area
    Ag = AISCshapes(snum,2);
    % Flange properties
    bf = AISCshapes(snum,5);    tf = AISCshapes(snum,7);
    % Stem properties
    d = AISCshapes(snum,3);     tw = AISCshapes(snum,6);
    % axis X-X properties
    Ix = AISCshapes(snum,15);   Zx = AISCshapes(snum,16);    
    Sx = AISCshapes(snum,17);   rx = AISCshapes(snum,18);    Sxc = Ix/AISCshapes(snum,9);
    % axis Y-Y properties
    Iy = AISCshapes(snum,19);   Zy = AISCshapes(snum,20);    
    Sy = AISCshapes(snum,21);   ry = AISCshapes(snum,22);
    % Torsional Properties
    J = AISCshapes(snum,26);    Cw = AISCshapes(snum,27);
    ro2 = (AISCshapes(snum,28))^2;   H = AISCshapes(snum,29);  
    Qs = AISCshapes(snum,30);
	% Compact section criteria (L-shape):
    lambda_f = AISCshapes(snum,13);  lambda_pf = 0.38*sqrt(E/Fy); lambda_rf = 1.0*sqrt(E/Fy);
    lambda_w = AISCshapes(snum,14);  lambda_pw = 0.84*sqrt(E/Fy); lambda_rw = 1.03*sqrt(E/Fy);
    if AISC16; lambfa_rw = 1.52*sqrt(E/Fy); end
    
	% Remaining flexural torsional properties
    Lp=1.76*ry*sqrt(E/Fy);
    Lr=1.95*(E/Fy)*sqrt(Iy*J)/Sx*sqrt(2.36*(Fy/E)*d*Sx/J+1);
    
    % Dummy variables
    rts = 0; ho = 0; rz = 0; a = 0;
    
    if Display % Display properties in two tables
    variablenames1={'Ag','d','tw','bf','tf','bf_2tf','H_tw','Sxc','ro2','H','Qs'};
    variablenames2={'Ix','Sx','rx','Zx','Iy','Sy','ry','Zy','J','Cw','Lp_ft','Lr_ft'};
    T1 = table(Ag,d,tw,bf,tf,lambda_f,lambda_w,Sxc,ro2,H,Qs,'VariableNames',variablenames1);
    T2 = table(Ix,Sx,rx,Zx,Iy,Sy,ry,Zy,J,Cw,Lp/12,Lr/12,'VariableNames',variablenames2);
    fprintf('Dimensions & Properties:\n');  disp(T1);   disp(T2)
    end

elseif all(slabel([1 2 3 4]) == 'WHEB') || all(slabel([1 2 3 4]) == 'WIPE')
    Ag = sect_info(elem_info(elem(1),3),1);
    Ix = sect_info(elem_info(elem(1),3),2); Zx = sect_info(elem_info(elem(1),3),6);  rx = sqrt(Ix/Ag);
    Iy = sect_info(elem_info(elem(1),3),3); Zy = sect_info(elem_info(elem(1),3),7);  ry = sqrt(Iy/Ag);
    J = sect_info(elem_info(elem(1),3),4); Cw = sect_info(elem_info(elem(1),3),5);
    bf = 1;tf = 1;ho = 1;tw = 1;lambda_f = 1;lambda_w = 1; %dummy variables
    
    lambda_pf = 0.38*sqrt(E/Fy); lambda_rf = 1.0*sqrt(E/Fy);
    lambda_pw = 3.76*sqrt(E/Fy); lambfa_rw = 5.70*sqrt(E/Fy);

    c = 1;
    Lp=1.76*ry*sqrt(E/Fy); % fprintf('Lp = 1.76*ry*sqrt(E/Fy)/12 = %1.2f ft.\n',Lp);
    Lr = 2*Lp; % assume Lr = 2*Lp
	variablenames={'Ag','Ix','rx','Zx','Iy','ry','Zy','J','Cw'};
    T = table(Ag,Ix,rx,Zx,Iy,ry,Zy,J,Cw,'VariableNames',variablenames);
    disp(T)

elseif slabel(1) == 'W'
    % Find the number of the property
	snum = find(strcmp(slabel,AISClabel),1);
    % Area
    Ag = AISCshapes(snum,2);
    % Flange properties
    bf = AISCshapes(snum,5); tf = AISCshapes(snum,7);
    % Web properties
    d = AISCshapes(snum,3);  tw = AISCshapes(snum,6); 
    ho = AISCshapes(snum,4);    xo = 0;    yo = 0;
    % axis X-X properties
    Ix = AISCshapes(snum,15);   Zx = AISCshapes(snum,16);    
    Sx = AISCshapes(snum,17);   rx = AISCshapes(snum,18);   Sxc = Sx;
    % axis Y-Y properties
    Iy = AISCshapes(snum,19);   Zy = AISCshapes(snum,20);    
    Sy = AISCshapes(snum,21);   ry = AISCshapes(snum,22);
    % Torsional Properties
    J = AISCshapes(snum,26);    Cw = AISCshapes(snum,27);
    % Compact section criteria (W-shape):
    lambda_f = AISCshapes(snum,13); lambda_pf = 0.38*sqrt(E/Fy); lambda_rf = 1.0*sqrt(E/Fy);
    lambda_w = AISCshapes(snum,14); lambda_pw = 3.76*sqrt(E/Fy); lambda_rw = 5.70*sqrt(E/Fy);
    
    % Compute remaining geometric properties:
    rts = sqrt(sqrt(Iy*Cw)/Sx); % fprintf('rts = sqrt(sqrt(I_y*C_w)/S_x) = %1.2f.\n',rts);
    Aw=Ag-(bf*tf); % web depth & shear area
    
    % Compute flexural strength geometric properties:
    c = ho/2*sqrt(Iy/Cw); 
    % fprintf('c = h_o/2*sqrt(I_y/C_w) = %1.2f.\n',c);
    Lp=1.76*ry*sqrt(E/Fy); 
    % fprintf('Lp = 1.76*ry*sqrt(E/Fy)/12 = %1.2f ft.\n',Lp);
    Lr=1.95*rts*(E/(0.7*Fy))*sqrt(J*c/(Sx*ho)+sqrt((J*c/(Sx*ho))^2+6.76*(0.7*Fy/E)^2));
    % fprintf('Lr = 1.95*rts*(E/(0.7*Fy))*sqrt(J*c/(Sx*ho)+sqrt((J*c/(Sx*ho))^2+6.76*(0.7*Fy/E)^2))/12 = %1.2f ft.\n\n',Lr);

	% Dummy variables
    ro2 = 0; H = 0; rz = 0; a = 0; Qs = 1;
    
    if Display % Display properties in two tables
    variablenames1={'Ag','d','tw','bf','tf','bf_2tf','H_tw','rts','ho'};
    variablenames2={'Ix','Sx','rx','Zx','Iy','Sy','ry','Zy','J','Cw','Lp_ft','Lr_ft'};
    T1 = table(Ag,d,tw,bf,tf,lambda_f,lambda_w,rts,ho,'VariableNames',variablenames1);
    T2 = table(Ix,Sx,rx,Zx,Iy,Sy,ry,Zy,J,Cw,Lp/12,Lr/12,'VariableNames',variablenames2);
    fprintf('Dimensions & Properties:\n');  disp(T1);   disp(T2)
    end
    
elseif all(slabel([1 2]) == '2L')
    % Find the number of the property
	snum = find(strcmp(slabel,AISClabel),1);
    % Find the number of the corresponding single angle
    sngl = find(slabel == 'X');
    if size(sngl,2)<3 && slabel(end) == 'B'; sngl = find(strcmp(slabel(2:end-4),AISClabel),1);
	elseif size(sngl,2)<3;  sngl = find(strcmp(slabel(2:end),AISClabel),1);
    else sngl = sngl(end); sngl = find(strcmp(slabel(2:sngl-1),AISClabel),1);
    end
    
    % Area
    Ag = AISCshapes(snum,2);
    % Long leg properties
    bf = AISCshapes(snum,5); tf = AISCshapes(snum,7);
    % Short leg properties
    d = AISCshapes(snum,3);  tw = tf;    ho = 0;    xo = 0;    yo = 0;
    % axis X-X properties
    Ix = AISCshapes(snum,15);   Zx = AISCshapes(snum,16);    
    Sx = AISCshapes(snum,17);   rx = AISCshapes(snum,18);    Sxc = Ix/AISCshapes(snum,9);
    % axis Y-Y properties
    Iy = AISCshapes(snum,19);   Zy = AISCshapes(snum,20);    
    Sy = AISCshapes(snum,21);   ry = AISCshapes(snum,22);
    % Torsional Properties
    J = 2*AISCshapes(sngl,26);    Cw = 2*AISCshapes(sngl,27);
    % Flexural-Torsional Properties
    rz = AISCshapes(sngl,24);     ro2 = (AISCshapes(snum,28))^2;   
    H = AISCshapes(snum,29);      Qs = AISCshapes(snum,30);
    b_t = AISCshapes(snum,13);
    if all(slabel(end-3:end)=='LLBB');  lambda_f = d/tf; lambda_w = b_t;
    else    lambda_f = b_t; lambda_w = d/tf;
    end
    % Compact section criteria (L-shape) for the long leg:
    lambda_pf = 0.54*sqrt(E/Fy); lambda_rf = 0.91*sqrt(E/Fy);
    % Compact section criteria (L-shape) for the short leg:
    lambda_pw = 0.54*sqrt(E/Fy); lambda_rw = 0.91*sqrt(E/Fy);

    % Remaining flexural torsional properties
    a = min(L/2,48);
    Lp=1.76*ry*sqrt(E/Fy);
    Lr=1.95*(E/(Fy))*sqrt(Iy*J)/Sx*sqrt(2.36*(Fy/E)*bf*Sx/J+1);
    
    if Display % Display properties in two tables
    variablenames1={'Ag','t','b','b_t','Sxc','ro2','H','rz','a','Qs'};
    variablenames2={'Ix','Sx','rx','Zx','Iy','Sy','ry','Zy','J','Cw','Lp_ft','Lr_ft'};
    T1 = table(Ag,tw,bf,b_t,Sxc,ro2,H,rz,a,Qs,'VariableNames',variablenames1);
    T2 = table(Ix,Sx,rx,Zx,Iy,Sy,ry,Zy,J,Cw,Lp/12,Lr/12,'VariableNames',variablenames2);
    fprintf('Dimensions & Properties:\n');  disp(T1);   disp(T2)
    end
    
elseif all(slabel([1 2]) == 'HS')
    % Reserved for future code

else
    fprintf(2,'Shape dimensions & properties not available in GL_AISCproperties.\n');
    % sect_info = [Ag, Izz, Iyy, J, Cw, Zzz, Zyy, Shear area Ayy, Shear area Azz, ...]
    Ag = sect_info(elem_info(elem(1),3),1);
    Ix = sect_info(elem_info(elem(1),3),2); Zx = sect_info(elem_info(elem(1),3),6);  rx = sqrt(Ix/Ag);
    Iy = sect_info(elem_info(elem(1),3),3); Zy = sect_info(elem_info(elem(1),3),7);  ry = sqrt(Iy/Ag);
    J = sect_info(elem_info(elem(1),3),4); Cw = sect_info(elem_info(elem(1),3),5);
    variablenames={'Ag','Ix','rx','Zx','Iy','ry','Zy','J','Cw'};
    T = table(Ag,Ix,rx,Zx,Iy,ry,Zy,J,Cw,'VariableNames',variablenames);
    
    bf = 1;tf = 1;ho = 1;tw = 1;lambda_f = 1;lambda_w = 1; %dummy variables\
    
    lambda_pf = 0.38*sqrt(E/Fy); lambda_rf = 1.0*sqrt(E/Fy);
    lambda_pw = 3.76*sqrt(E/Fy); lambfa_rw = 5.70*sqrt(E/Fy);

    c = 1;
    Lp=1.76*ry*sqrt(E/Fy); % fprintf('Lp = 1.76*ry*sqrt(E/Fy)/12 = %1.2f ft.\n',Lp);
    Lr = 2*Lp; % assume Lr = 2*Lp
    disp(T)
end