function []= BeamDesign()
sc = get(0,'ScreenSize');
BeamDesignGui = figure('Visible','off','MenuBar','none','Toolbar','none',...
    'Position',[sc(3)/2-110 sc(4)/2-90, 400, 130],'name','Create Node');

%  Construct the components.
beam.Button1 = uicontrol('Style','pushbutton','String','Input Properties',...
           'Position',[40, 40 ,120,30],'Callback',{@hpropbutton_Callback});
beam.Button2 = uicontrol('Style','pushbutton','String','Input Section Name',...
           'Position',[240 40 ,120,30],'Callback',{@hsectbutton_Callback});
beam.text1 = uicontrol('Style','text','String','or',...
           'Position',[185,45,30,15]);
beam.text2 = uicontrol('Style','text','String','Steel Flexural Strength',...
           'Position',[50,90,300,15]);
%Make the UI visible.
set(BeamDesignGui,'Visible','on')
set(beam.Button1,'callback',{@Button_call1,beam}) % Set the callback, pass hands.
set(beam.Button2,'callback',{@Button_call2,beam}) % Set the callback, pass hands.
end

function Button_call1(varargin)% Callback for pushbutton.
    load('GL_AISCproperties.mat'); 
    load('GiesenLoo_AISCname.mat');
    dlg_title = 'Input';
    num_lines = 1;
    prompt = {'Modulus of Elasticity (ksi):','Steel grade Fy (ksi):','X-sectional Area (in^2):',...
                'Design depth of the beam (in):','Web thickness (in):','Flange width (in):',...
                'Flange thickness (in):','b_f/(2*t_f):   ','h_o/t_w :   '};
    def = {'29000','50','16.2','23.6','0.395','7.01','0.505','6.94','54.6'}; % default values of each input
    answer = inputdlg(prompt, dlg_title, num_lines, def); % open dialog box
    if ~isempty(answer) % only if OK button was clicked
        E = str2double(answer{1});  G = E/(2*(1+0.3));
        Fy = str2double(answer{2});
        Ag = str2double(answer{3});
        d = str2double(answer{4});
        tw = str2double(answer{5});
        bf = str2double(answer{6});
        tf = str2double(answer{7});       
        ho=d-tf;  %distance between flange centroids
        Aw=Ag-2*(bf*tf);
        lambda_f = str2double(answer{8}); lambda_pf = 0.38*sqrt(E/Fy); lambda_rf = 1.0*sqrt(E/Fy);
        lambda_w = str2double(answer{9}); lambda_pw = 3.76*sqrt(E/Fy);  lambda_rw = 5.70*sqrt(E/Fy);
    end
    
    prompt = {'Elastic x-axis Section modulus Sx :',...
	'Plastic x-axis Section Modulus Zx :',...
    'Moment of inertia Iy :','Elastic y-axis Section modulus Sy  :',...
    'Plastic y-axis Section Modulus Zy :','Polar moment of inertia J :',...
    'Warping torsional property Cw :','unbraced length :','Enter C_b value:   '};
    def = {'114','134','29.1','8.3','13.3','1.18','3870','360','1'}; % default values of each input
    answer = inputdlg(prompt, dlg_title, num_lines, def); % open dialog box
    if ~isempty(answer) % only if OK button was clicked
        Sx = str2double(answer{1});
        Zx = str2double(answer{2});
        Iy = str2double(answer{3});    ry = sqrt(Iy/Ag);
        Sy = str2double(answer{4});
        Zy = str2double(answer{5});
        J = str2double(answer{6});
        Cw = str2double(answer{7});
        Lb = str2double(answer{8});
        Cb = str2double(answer{9});
    end
    
    rts = sqrt(sqrt(Iy*Cw)/Sx); fprintf('rts = sqrt(sqrt(I_y*C_w)/S_x) = %1.2f.\n',rts);
    c = ho/2*sqrt(Iy/Cw); fprintf('c = h_o/2*sqrt(I_y/C_w) = %1.2f.\n',c);
    Lp=1.76*ry*sqrt(E/Fy)/12; fprintf('L_p = 1.76*r_y*sqrt(E/Fy)/12 = %1.2f ft.\n',Lp);
    Lr=1.95*rts*(E/(0.7*Fy))*sqrt(J*c/(Sx*ho)+sqrt((J*c/(Sx*ho))^2+6.76*(0.7*Fy/E)^2))/12;
    fprintf('L_r = 1.95*rts*(E/(0.7*Fy))*sqrt(J*c/(S_x*h_o)+sqrt((J*c/(S_x*h_o))^2+6.76*(0.7*Fy/E)^2))/12 = %1.2f ft.\n\n',Lr);

    %% Display properties in two tables nicely
    variablenames1={'A_in2','d_in','t_w_in','b_f_in','t_f_in','h_o_in','lambda_f','lambda_p'};
    variablenames2={'Sx_in3','Zx_in3','Iy_in4','ry_in','J_in4','Cw_in6','L_p_ft','L_r_ft'};
    T1 = table(Ag,d,tw,bf,tf,ho,lambda_f,lambda_w,'VariableNames',variablenames1);
    T2 = table(Sx,Zx,Iy,ry,J,Cw,Lp,Lr,'VariableNames',variablenames2);
    fprintf('Geometric Properties:\n');  disp(T1);   disp(T2)
    
    %% FlexuralStrength
    [phiMnx,phiMny] = AISC_FlexuralStrength(E,Fy,Zx,Sx,Zy,Sy,Iy,J,Cw,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,d)
end

function Button_call2(varargin)% Callback for pushbutton.
    dlg_title = 'Input';
    num_lines = 1;
    prompt = {'Section Name, e.g. "W8X35"','Modulus of Elasticity','Yield Stress','Buckling length','Cb'};
    def = {'W8X35','29000','50','360','1'};
    answer = inputdlg(prompt, dlg_title, num_lines, def);
    if ~isempty(answer) % only if OK button was clicked
        slabel = answer{1};
        E = str2double(answer{2});
        Fy = str2double(answer{3});
        Lb = str2double(answer{4});
        Cb = str2double(answer{5});
    end
	GiesenLoo_AISCproperties
    [phiMnx,phiMny] = AISC_FlexuralStrength(E,Fy,Zx,Sx,Zy,Sy,Iy,J,Cw,Lb,Lp,Lr,Cb,...
     lambda_f,lambda_pf,lambda_rf,lambda_w,lambda_pw,lambda_rw,d)
end
