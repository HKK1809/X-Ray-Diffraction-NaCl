% NaCl_100, hkl = 002
Theta100 = NaCl100(:,1);
Theta100 = table2array(Theta100);

Imp_s_100 = NaCl100(:,2);
Imp_s_100 = table2array(Imp_s_100);

lamda_b = (12398/8903)*(10^(-10));
lamda_a = (12398/8036.5)*(10^(-10));

figure(1)
plot(Theta100, Imp_s_100)
xlabel('Angle (Degrees)','FontSize',18);
ylabel('Pulse Count','FontSize',18);
title ('Number of photons Observed with Respect to Angle for NaCl(100) ','FontSize',18)

% For then gauss fit
    %n = 1 
    figure(2)
    Theta_100_1 = Theta100(90:150,1);
    Imp_s_100_1 = Imp_s_100(90:150,1);
    [F_100_1, gof_100_1, fit_output_100_1] = fit( Theta_100_1, Imp_s_100_1, 'gauss2');
    plot( F_100_1,Theta_100_1, Imp_s_100_1);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(100),n = 1','FontSize',18)


    %n = 2 
    figure(3)
    Theta_100_2 = Theta100(240:330,1);
    Imp_s_100_2 = Imp_s_100(240:330,1);
    [F_100_2,gof_100_2,fit_output_100_2] = fit( Theta_100_2, Imp_s_100_2, 'gauss2');
    plot( F_100_2,Theta_100_2, Imp_s_100_2);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(100),n = 2 ','FontSize',18)

% Calculation spacing between two faces 'd'
% First set of peaks, n=1
% Beta
thetamax100b1 = (F_100_1.b2*pi)/180;  %In radians
d_b100_1 = lamda_b/(2*sin(thetamax100b1)); % d = ((n*lamda)/(2*sin(theta)))
b_100_1 = d_b100_1*2; % a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_100b1 = (F_100_1.c2*pi)/180;
    err_d_b100_1 = err_Theta_rad_100b1*((lamda_b*cos(thetamax100b1))/(2*sin(thetamax100b1)*sin(thetamax100b1)));
    err_b_100_1 = 2*err_d_b100_1;

%Alpha
thetamax100a1 = (F_100_1.b1*pi)/180;  %In radians
d_a100_1 = (lamda_a)/(2*sin(thetamax100a1));
a_100_1 = d_a100_1*2; %a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_100a1 = (F_100_1.c1*pi)/180;
    err_d_a100_1 = err_Theta_rad_100a1*((lamda_a*cos(thetamax100a1))/(2*sin(thetamax100a1)*sin(thetamax100a1)));
    err_a_100_1 = 2* err_d_a100_1;

% Mean error
    err_d_100_1 = (sqrt((err_d_a100_1*err_d_a100_1)+(err_d_b100_1*err_d_b100_1)))/2;
    err_ab_100_1 = (sqrt((err_a_100_1*err_a_100_1)+(err_b_100_1*err_b_100_1)))/2;

% Second set of peaks, n=2
% Beta
thetamax100b_2 = (F_100_2.b2*pi)/180;  %In radians
d_b100_2 = (2*lamda_b)/(2*sin(thetamax100b_2)); % d = ((n*lamda)/(2*sin(theta)))
b_100_2 = d_b100_2*2; %a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_100b2 = (F_100_2.c2*pi)/180;
    err_d_b100_2 = err_Theta_rad_100b2*((2*lamda_b*cos(thetamax100b_2))/(2*sin(thetamax100b_2)*sin(thetamax100b_2)));
    err_b_100_2 = 2* err_d_b100_2;

%Alpha
% Second set of peaks
thetamax100a_2 = (F_100_2.b1*pi)/180;  %In radians
d_a100_2 = (2*lamda_a)/(2*sin(thetamax100a_2));
a_100_2 = d_a100_2*2; %a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_100a2 = (F_100_2.c1*pi)/180;
    err_d_a100_2 = err_Theta_rad_100a2*((2*lamda_a*cos(thetamax100a_2))/(2*sin(thetamax100a_2)*sin(thetamax100a_2)));
    err_a_100_2 = 2* err_d_a100_2;

% Mean error
err_d_100_2 = (sqrt((err_d_a100_2*err_d_a100_2)+(err_d_b100_2*err_d_b100_2)))/2;
err_ab_100_2 = (sqrt((err_a_100_2*err_a_100_2)+(err_b_100_2*err_b_100_2)))/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NaCl_110, hkl = 022
Theta110 = NaCl110(:,1);
Theta110 = table2array(Theta110);

Imp_s_110 = NaCl110(:,2);
Imp_s_110 = table2array(Imp_s_110);

lamda_b = (12398/8903)*(10^(-10));
lamda_a = (12398/8036.5)*(10^(-10));

figure(4)
plot(Theta110, Imp_s_110)
xlabel('Angle (Degrees)','FontSize',18);
ylabel('Pulse Count','FontSize',18);    
title ('Number of photons Observed with Respect to Angle for NaCl(110) ','FontSize',18)


% For then gauss fit
    %n = 1 
    figure(5)
    Theta_110_1 = Theta110(150:220,1);
    Imp_s_110_1 = Imp_s_110(150:220,1);
    [F_110_1, gof_110_1, fit_output_110_1] = fit( Theta_110_1, Imp_s_110_1, 'gauss2');
    plot( F_110_1,Theta_110_1, Imp_s_110_1);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(110), n = 1 ','FontSize',18)

    %n = 2 
    figure(6)
    Theta_110_2 = Theta110(400:500,1);
    Imp_s_110_2 = Imp_s_110(400:500,1);
    [F_110_2,gof_110_2,fit_output_110_2] = fit( Theta_110_2, Imp_s_110_2, 'gauss2');
    plot( F_110_2,Theta_110_2, Imp_s_110_2);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(110), n = 2 ','FontSize',18)


% Calculation spacing between two faces 'd'
% First set of peaks, n=1
% Beta
root8 = sqrt(8);
thetamax110b1 = (20.3*pi)/180;  %In radians
d_b110_1 = lamda_b/(2*sin((thetamax110b1))); % d = ((n*lamda)/(2*sin(theta)))
b_110_1 = d_b110_1*root8; % a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_110b1 = (F_110_1.c2*pi)/180;
    err_d_b110_1 = err_Theta_rad_110b1*((lamda_b*cos(thetamax110b1))/(2*sin(thetamax110b1)*sin(thetamax110b1)));
    err_b_110_1 = root8*err_d_b110_1;

%Alpha
thetamax110a1 = (22.6*pi)/180;  %In radians
d_a110_1 = (lamda_a)/(2*sin((thetamax110a1)));
a_110_1 = d_b110_1*root8; %a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_110a1 = (F_110_1.c1*pi)/180;
    err_d_a110_1 = err_Theta_rad_110a1*((lamda_a*cos(thetamax110a1))/(2*sin(thetamax110a1)*sin(thetamax110a1)));
    err_a_110_1 = root8*err_d_a110_1;

% Mean error
err_d_110_1 = (sqrt((err_d_a110_1*err_d_a110_1)+(err_d_b110_1*err_d_b110_1)))/2;
err_ab_110_1 = (sqrt((err_a_110_1*err_a_110_1)+(err_b_110_1*err_b_110_1)))/2;

% Second set of peaks, n=2
% Beta
root2 = sqrt(2);
thetamax110b_2 = (F_110_2.b2*pi)/180;  %In radians
d_b110_2 = (2*lamda_b)/(2*sin(thetamax110b_2)); % d = ((n*lamda)/(2*sin(theta))) 
b_110_2 = d_b110_2*root8; % a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_110b2 = (F_110_2.c2*pi)/180;
    err_d_b110_2 = err_Theta_rad_110b2*((2*lamda_b*cos(thetamax110b_2))/(2*sin(thetamax110b_2)*sin(thetamax110b_2)));
    err_b_110_2 = root8*err_d_b110_2;

%Alpha
% Second set of peaks
thetamax110a_2 = (F_110_2.b1*pi)/180;  %In radians
d_a110_2 = (2*lamda_a)/(2*sin(thetamax110a_2));
a_110_2 = d_a110_2*root8; %a = d*sqrt((h*h)+(k*k)+(l*l))
    %Error Analysis
    err_Theta_rad_110a2 = (F_110_2.c1*pi)/180;
    err_d_a110_2 = err_Theta_rad_110a2*((2*lamda_a*cos(thetamax110a_2))/(2*sin(thetamax110a_2)*sin(thetamax110a_2)));
    err_a_110_2 = root8*err_d_a110_2;

% Mean error
err_d_110_2 = (sqrt((err_d_a110_2*err_d_a110_2)+(err_d_b110_2*err_d_b110_2)))/2;
err_ab_110_2 = (sqrt((err_a_110_2*err_a_110_2)+(err_b_110_2*err_b_110_2)))/2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% NaCl_111
Theta111 = NaCl111(:,1);
Theta111 = table2array(Theta111);

Imp_s_111 = NaCl111(:,2);
Imp_s_111 = table2array(Imp_s_111);

figure(7)
plot(Theta111, Imp_s_111)
ft = fittype("a1*exp(-((x-b1)/c1)^2)+d1");
xlabel('Angle (Degrees)','FontSize',18);
ylabel('Pulse Count','FontSize',18);    
title ('Number of photons Observed with Respect to Angle for NaCl(111) ','FontSize',18)

% For then gauss fit
    %n = 1 
    % alpha
    figure(8)
    Theta_a111_1 = Theta111(90:110,1);
    Imp_s_a111_1 = Imp_s_111(90:110,1);
    [F_a111_1, gof_a111_1, fit_output_a111_1] = fit( Theta_a111_1, Imp_s_a111_1, ft,"StartPoint",[140,13,1,60]);
    plot( F_a111_1,Theta_a111_1, Imp_s_a111_1);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(111), n = 1, Alpha Peak','FontSize',18)

    % beta
    figure(9)
    Theta_b111_1 = Theta111(78:95,1);
    Imp_s_b111_1 = Imp_s_111(78:95,1);
    [F_b111_1, gof_b111_1, fit_output_b111_1] = fit( Theta_b111_1, Imp_s_b111_1, ft,"StartPoint",[100,11.5,0.2,60]);
    plot( F_b111_1,Theta_b111_1, Imp_s_b111_1);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(100), n = 1, Beta Peak ','FontSize',18)

    %n = 2 
    figure(10)
    Theta_111_2 = Theta111(190:275,1);
    Imp_s_111_2 = Imp_s_111(190:275,1);
    [F_111_2,gof_111_2,fit_output_111_2] = fit( Theta_111_2, Imp_s_111_2, 'gauss2');
    plot( F_111_2,Theta_111_2, Imp_s_111_2);
    xlabel('Angle (Degrees)','FontSize',18);
    ylabel('Pulse Count','FontSize',18);    
    title ('Gauss fit for the Number of photons Observed with Respect to Angle for NaCl(111), n = 2','FontSize',18)

% Calculation spacing between two faces 'd'
% First set of peaks, n=1
% Beta
root3 = sqrt(3);
thetamax111b1 = (F_b111_1.b1*pi)/180;  %In radians
d_b111_1 = lamda_b/(2*sin(thetamax111b1)); % d = ((n*lamda)/(2*sin(theta)))
b_111_1 = d_b111_1*root3; % a = d*sqrt((h*h)+(k*k)+(l*l))
        %Error Analysis
        err_Theta_rad_111b1 = (F_b111_1.c1*pi)/180;
        err_d_b111_1 = err_Theta_rad_111b1*((lamda_b*cos(thetamax111b1))/(2*sin(thetamax111b1)*sin(thetamax111b1)));
        err_b_111_1 = root3*err_d_b111_1;

%Alpha
thetamax111a1 = (F_a111_1.b1*pi)/180;  %In radians
d_a111_1 = (lamda_a)/(2*sin(thetamax111a1));
a_111_1 = d_b111_1*root3; %a = d*sqrt((h*h)+(k*k)+(l*l))
        %Error Analysis
        err_Theta_rad_111a1 = (F_a111_1.c1*pi)/180;
        err_d_a111_1 = err_Theta_rad_111a1*((lamda_a*cos(thetamax111a1))/(2*sin(thetamax111a1)*sin(thetamax111a1)));
        err_a_111_1 = root3*err_d_a111_1;

% Mean error
err_d_111_1 = (sqrt((err_d_a111_1*err_d_a111_1)+(err_d_b111_1*err_d_b111_1)))/2;
err_ab_111_1 = (sqrt((err_a_111_1*err_a_111_1)+(err_b_111_1*err_b_111_1)))/2;

% Second set of peaks, n=2
% Beta
root3 = sqrt(3);
thetamax111b_2 = (F_111_2.b2*pi)/180;  %In radians
d_b111_2 = (2*lamda_b)/(2*sin(thetamax111b_2)); % d = ((n*lamda)/(2*sin(theta)))
b_111_2 = d_b111_2*root3; % a = d*sqrt((h*h)+(k*k)+(l*l))
        %Error Analysis
        err_Theta_rad_111b2 = (F_111_2.c2*pi)/180;
        err_d_b111_2 = err_Theta_rad_111b2*((2*lamda_b*cos(thetamax111b_2))/(2*sin(thetamax111b_2)*sin(thetamax111b_2)));
        err_b_111_2 = root3*err_d_b111_2;

%Alpha
% Second set of peaks
thetamax111a_2 = (F_111_2.b1*pi)/180;  %In radians
d_a111_2 = (2*lamda_a)/(2*sin(thetamax111a_2));
a_111_2 = d_a111_2*root3; %a = d*sqrt((h*h)+(k*k)+(l*l))
        %Error Analysis
        err_Theta_rad_111a2 = (F_111_2.c1*pi)/180;
        err_d_a111_2 = err_Theta_rad_111a2*((2*lamda_a*cos(thetamax111a_2))/(2*sin(thetamax111a_2)*sin(thetamax111a_2)));
        err_a_111_2 = root3*err_d_a111_2;

% Mean error
err_d_111_2 = (sqrt((err_d_a111_2*err_d_a111_2)+(err_d_b111_2*err_d_b111_2)))/2;
err_ab_111_2 = (sqrt((err_a_111_2*err_a_111_2)+(err_b_111_2*err_b_111_2)))/2;




