%% Trajectory Calcs

clear all
close all

%% Constants
% Target Alt = 343 km, Target Vel = 3387 m/s
GM = 4.282837e13; % m3s-2

Thrust_S1 = 5500; %N
Vel_exh_S1 = 3214.38; %ms-1
m_fuel_dot_S1 = 1.792; %kg/s
burn_1_time = 113;

Thrust_S2 = 2000; %N
Vel_exh_S2 = 2893.95; %ms-1
m_fuel_dot_S2 = 0.651428; %kg/s
burn_2_time = 90; %s

Stage_1_Inert_Mass = 93.6; %kg



% Sim Variables
% Sep_time = [210]; %s 
% Gamma_PO = deg2rad([72]);

%Sep_time = [210, 240, 360];
%Gamma_PO = deg2rad([72,69,64]);

Sep_time = [30:5:480]; %s 
Gamma_PO = deg2rad([20:5:85]);

%Sep_time = 150;
%Gamma_PO = 0.9599;
% Border Values

MECO1 = burn_1_time;
S2_Ignition = burn_1_time + Sep_time;
%% Initial Values

M_0 = 400; %kg
r_0 = 3396200; %m
V_0 = Thrust_S1/M_0 - GM/(r_0^2); %m/s
Alt = 0; %m
Theta_0 = 0; % rad
Gamma_0 =  pi/2; % rad
Ground_Pass_0 = 0;
TVC_0 = 0;

% Total_Gamma = deg2rad(85);
% Gamma_rate = Total_Gamma/(burn_1_time+burn_2_time - 10);

%% Analtics Setup
SimLength = 2000;
Results = zeros(numel(Sep_time),numel(Gamma_PO)); % 1 undershoot, 2 overshoot, 3 Success
figure(1)
hold on
grid on
xlabel('Time (s)')
ylabel('Altitude (m)')
yline(343500,':');
yline(342500,':');

Pos_mat = zeros(numel(Sep_time)*numel(Gamma_PO),4);

figure(2)
hold on
grid on
xlabel('Time (s)')
ylabel('Velocity (m/s)')
yline(3382,':');
yline(3392,':');

figure(3)
hold on
grid on
xlabel('Time (s)')
ylabel('Gamma (deg)')

figure(4)
hold on
grid on
xlabel('Velocity (m/s)')
ylabel('Altitude (m)')
yline(342000,'--');
yline(344000,'--');
xline(3377,'--');
xline(3397,'--');

figure(5)
hold on
grid on
xlabel('Altitdue (m)')
ylabel('Gamma (deg)')
xline(342000,'--');
xline(344000,'--');

figure(6)
subplot(2,1,1)
hold on
grid on
ylabel('Altitude (m)')
yline(343500,':');
yline(342500,':');

subplot(2,1,2)
hold on
grid on
ylabel('Velocity (m/s)')
yline(3382,':');
yline(3392,':');

figure(7)
hold on
grid on
ylabel('Rate of Change of Flight Path Angle (rad/s)')
xlabel('Time (s)')

figure(8)
hold on
grid on
title('Joe''s stupid graph')

it = 0;
%% Loop
for cx = 1:numel(Sep_time)
    for cc = 1:numel(Gamma_PO)
        it = it+1;
        time = zeros(SimLength,1);
        Altitude = zeros(SimLength,1);
        Velocity = zeros(SimLength,1);
        Mass = zeros(SimLength,1);
        r = zeros(SimLength,1);
        g = zeros(SimLength,1);
        Accel = zeros(SimLength,1);
        Delta_Gamma = zeros(SimLength,1);
        Gamma = zeros(SimLength,1);
        Theta = zeros(SimLength,1);
        Ground_Pass = zeros(SimLength,1);
        TVC = zeros(SimLength,1);
        Grav_Gamma = zeros(SimLength,1);
        TVC_Mag = zeros(SimLength,1);
        TVC_Mag_Lim = zeros(SimLength,1);
        TVC_Flag = zeros(SimLength,1);
        TVC_Act = zeros(SimLength,1);
        for t = 1:SimLength

            % First Burn
            if t<burn_1_time
                if t == 1
                    time(t) = t;
                    Altitude(t) = Alt;
                    Velocity(t) = V_0;
                    Mass(t) = M_0;
                    r(t) = r_0;
                    g(t) = GM/(r_0^2);
                    Accel(t) = 0;
                    Gamma(t) = Gamma_0;
                    Delta_Gamma(t) = 0;
                    Theta(t) = Theta_0;
                    Ground_Pass(t) = Ground_Pass_0;
                    TVC(t) = TVC_0;
                    Grav_Gamma(t) = 0;
                    TVC_Mag(t) = 0;
                    TVC_Mag_Lim(t) = 0;
                    TVC_Flag(t) = 0;
                    TVC_Act(t) = 0;

                % Pre-PO, no gamma, no forced gamma
                elseif t<10     
                    time(t) = t;
                    Altitude(t) = Altitude(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                    Velocity(t) = Velocity(t-1) + (m_fuel_dot_S1*Vel_exh_S1/Mass(t-1)) - (g(t-1)*sin(Gamma(t-1)));
                    Mass(t) = Mass(t-1)-m_fuel_dot_S1;
                    r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                    g(t) = GM/(r(t-1)^2);
                    Accel(t) = ((m_fuel_dot_S1*Vel_exh_S1/Mass(t-1)) - (g(t-1)*Gamma(t-1)))/9.81;
                    Delta_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                    Gamma(t) = Gamma(t-1) + Delta_Gamma(t);  
                    %Gamma(t) = Gamma(t-1) - Gamma_rate;
                    Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                    Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                    TVC(t) = 0;
                    Grav_Gamma(t) = Delta_Gamma(t);
                    TVC_Mag(t) = 0;
                    TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                    

                    %Pitch Over - Set initial pitch over angle 
                elseif t==10
                    time(t) = t;
                        Altitude(t) = Altitude(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                        Velocity(t) = Velocity(t-1) + (m_fuel_dot_S1*Vel_exh_S1/Mass(t-1)) - (g(t-1)*sin(Gamma(t-1)));
                        Mass(t) = Mass(t-1)-m_fuel_dot_S1;
                        r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                        g(t) = GM/(r(t-1)^2);
                        Accel(t) = ((m_fuel_dot_S1*Vel_exh_S1/Mass(t-1)) - (g(t-1)*Gamma(t-1)))/9.81;
                        Gamma(t) = Gamma_PO(cc);
                        %Delta_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                        Delta_Gamma(t) = Gamma(t) - Gamma(t-1);
                        %Gamma(t) = Gamma(t-1) - Gamma_rate;
                        Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                        Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                        TVC(t) = 1;
                        Grav_Gamma(t) = 0;
                        TVC_Mag(t) = abs(Delta_Gamma(t));
                        TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                        

                else
                    if Altitude(t-1)<0
                        Results(cx,cc) = 4;
                        break
                    end
                    time(t) = t;
                    Altitude(t) = Altitude(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                    Velocity(t) = Velocity(t-1) + (m_fuel_dot_S1*Vel_exh_S1/Mass(t-1)) - (g(t-1)*sin(Gamma(t-1)));
                    Mass(t) = Mass(t-1)-m_fuel_dot_S1;
                    r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                    g(t) = GM/(r(t-1)^2);
                    Accel(t) = ((m_fuel_dot_S1*Vel_exh_S1/Mass(t-1)) - (g(t-1)*Gamma(t-1)))/9.81;
                    Delta_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                    Gamma(t) = Gamma(t-1) + Delta_Gamma(t);  
                    %Gamma(t) = Gamma(t-1) - Gamma_rate;
                    Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                    Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                    TVC(t) = 0;
                    Grav_Gamma(t) = Delta_Gamma(t);
                    TVC_Mag(t) = 0;
                   TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                end

            % Stage Seperation - Stop mass flow, remove inert mass,     
            elseif t<burn_1_time + Sep_time(cx)
                tsep = burn_1_time + Sep_time(cx);
                time(t) = t;
                Altitude(t) = Altitude(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                Velocity(t) = Velocity(t-1) - (g(t-1)*sin(Gamma(t-1)));
                if t == burn_1_time + (Sep_time(cx)/2)
                    Mass(t) = Mass(t-1) - Stage_1_Inert_Mass;
                else
                    Mass(t) = Mass(t-1);
                end
                r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                g(t) = GM/(r(t-1)^2);
                %Gamma(t) = Gamma(t-1) + (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1))); 
                %Gamma(t) = deg2rad(10);
                Delta_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                Gamma(t) = Gamma(t-1) + Delta_Gamma(t);
                Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                TVC(t) = 0;
                Grav_Gamma(t) = Delta_Gamma(t);
                TVC_Mag(t) = 0;
                TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                
            % Stage 2 Burn  
            else 
                if Altitude(t-1)<0
                    Results(cx,cc) = 4;
                    break
                end
                
                
                if Altitude(t-1) < 343500 && Altitude(t-1) > 342500 && Velocity(t-1) < 3407 && Velocity(t-1) > 3367 && Gamma(t-1) == 0
                    disp('Sim Successful')
                    Results(cx,cc) = 3;
                    Pos_mat(it,1) = it;
                    Pos_mat(it,2) = Altitude(t);
                    Pos_mat(it,3) = Velocity(t);
                    Pos_mat(it,4) = Gamma(t);
                    
                    Sep_time(cx)
                    rad2deg(Gamma_PO(cc))
                    Mass(t-1) - 47.4
                    burn_2_time - (t-1 - burn_1_time - Sep_time(cx))
                    BurnCheck = t-1 - burn_1_time - Sep_time(cx)
                    
                    figure(8)
                    scatter(Sep_time(cx),BurnCheck)
                    
                    Altitude(t-1)
                    Velocity(t-1)
                    Gamma(t-1);
                    break
                elseif Altitude(t-1) > 344000 && Velocity(t-1) > 3397
                    %disp('Sim Overshot')
                    Results(cx,cc) = 2;
                    Pos_mat(it,1) = it;
                    Pos_mat(it,2) = Altitude(t-1);
                    Pos_mat(it,3) = Velocity(t-1);
                    Pos_mat(it,4) = Gamma(t-1);
                    break
                    
%                   % Not used??
                elseif Altitude(t-1) < 342700 && Altitude(t-1) > 343300 && Velocity(t-1) < 3407 && Velocity(t-1) > 3367
                       time(t) = t;
                       Altitude(t) = Altitude(t-1) + Velocity(t-1)*sin(Gamma(t-1));
                       Velocity(t) = Velocity(t-1);
                       Mass(t) = Mass(t-1);
                       r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                       g(t) = GM/(r(t-1)^2);
                       Accel(t) = ((m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1))/9.81;
                       %Gamma(t) = Gamma(t-1) + (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1))); 
                       Gamma(t) = deg2rad(0);
                       Delta_Gamma(t) = Gamma(t) - Gamma(t-1);
                       Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                       Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                       TVC(t) = 1;
                       Grav_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                       TVC_Mag(t) = abs(Delta_Gamma(t) - Grav_Gamma(t));
                      TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                elseif Altitude(t-1) > 325000 && Altitude(t-1)<342000
%                        disp('Correction created')
                       time(t) = t;
                       Altitude(t) = Altitude(t-1) + Velocity(t-1)*sin(Gamma(t-1));
                       if Velocity(t-1) > 3407
                           Velocity(t) = Velocity(t-1);
                           Mass(t) = Mass(t-1);
                       else
                        Velocity(t) = Velocity(t-1) + (m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1);
                        Mass(t) = Mass(t-1)-m_fuel_dot_S2;
                       end
                       r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                       g(t) = GM/(r(t-1)^2);
                       Accel(t) = ((m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1))/9.81;
                       %Gamma(t) = Gamma(t-1) + (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1))); 
                       Gamma(t) = deg2rad(5);
                       Delta_Gamma(t) = Gamma(t) - Gamma(t-1);
                       Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                       Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                       TVC(t) = 1;
                       Grav_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                       TVC_Mag(t) = abs(Delta_Gamma(t) - Grav_Gamma(t));
                       TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
               elseif Altitude(t-1) > 300000 && Altitude(t-1)<325000
                       time(t) = t;
                       Altitude(t) = Altitude(t-1) + Velocity(t-1)*sin(Gamma(t-1));
                       if Velocity(t-1) > 3407
                           Velocity(t) = Velocity(t-1);
                           Mass(t) = Mass(t-1);
                       else
                        Velocity(t) = Velocity(t-1) + (m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1);
                        Mass(t) = Mass(t-1)-m_fuel_dot_S2;
                       end
                       r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                       g(t) = GM/(r(t-1)^2);
                       Accel(t) = ((m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1))/9.81;
                       %Gamma(t) = Gamma(t-1) + (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1))); 
                       Gamma(t) = deg2rad(20);
                       Delta_Gamma(t) = Gamma(t) - Gamma(t-1);
                       Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                       Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                       TVC(t) = 1;
                       Grav_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                       TVC_Mag(t) = abs(Delta_Gamma(t) - Grav_Gamma(t));
                      TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                elseif Altitude(t-1) >342000 && Altitude(t-1) <343000
                     
                       time(t) = t;
                       Altitude(t) = Altitude(t-1) + Velocity(t-1)*sin(Gamma(t-1));
                       if Velocity(t-1) > 3407
                           Velocity(t) = Velocity(t-1);
                           Mass(t) = Mass(t-1);
                       else
                        Velocity(t) = Velocity(t-1) + (m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1);
                        Mass(t) = Mass(t-1)-m_fuel_dot_S2;
                       end
                       r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                       g(t) = GM/(r(t-1)^2);
                       Accel(t) = ((m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1))/9.81;
                       %Gamma(t) = Gamma(t-1) + (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1))); 
                       Gamma(t) = deg2rad(0);
                       Delta_Gamma(t) = Gamma(t) - Gamma(t-1);
                       Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                       Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                       TVC(t) = 1;
                       Grav_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                       TVC_Mag(t) = abs(Delta_Gamma(t) - Grav_Gamma(t));
                       TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                else
                    time(t) = t;
                    Altitude(t) = Altitude(t-1) + Velocity(t-1)*sin(Gamma(t-1));
                    Velocity(t) = Velocity(t-1) + (m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1);
                    Mass(t) = Mass(t-1)-m_fuel_dot_S2;
                    r(t) = r(t-1) + (Velocity(t-1)*sin(Gamma(t-1)));
                    g(t) = GM/(r(t-1)^2);
                    Accel(t) = ((m_fuel_dot_S2*Vel_exh_S2/Mass(t-1)) - g(t-1))/9.81;
                    Delta_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                    %Gamma(t) = Gamma(t-1) + (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1))); 
                    Gamma(t) = Gamma(t-1) + Delta_Gamma(t);
                    Theta(t) = Theta(t-1) + Velocity(t-1)*cos(Gamma(t-1))/r(t-1);
                    Ground_Pass(t) = Ground_Pass(t-1) + Velocity(t-1)*cos(Gamma(t-1));
                    TVC(t) = 0;
                    Grav_Gamma(t) = (Velocity(t-1)*cos(Gamma(t-1))/r(t-1)) - ((1/Velocity(t-1))*(g(t-1))*cos(Gamma(t-1)));
                    TVC_Mag(t) = 0;
                   TVC_Mag_Lim(t) = TVC_Mag_Lim(t-1) + TVC_Mag(t) - TVC_Act(t-1);
                        if TVC_Mag_Lim(t) > deg2rad(5)
                            TVC_Act(t) = deg2rad(5);
                       elseif TVC_Mag_Lim(t) > deg2rad(1.5)
                           TVC_Act(t) = deg2rad(1.5);
                       elseif TVC_Mag_Lim(t) > deg2rad(0.5)
                           TVC_Act(t) = deg2rad(0.5);
                       else
                            TVC_Act(t) = TVC_Mag_Lim(t);
                       end
                        if TVC_Act(t) > 0
                            TVC_Flag(t) = 1;
                        else
                            TVC_Flag(t) = 0;
                        end
                end
                
                if t>(burn_1_time+burn_2_time+Sep_time(cx)) %Stop sim due to empty tanks
                    %disp('Sim Undershot')
                    Results(cx,cc) = 1;   
                    Pos_mat(it,1) = it;
                    Pos_mat(it,2) = Altitude(t-1);
                    Pos_mat(it,3) = Velocity(t-1);
                    Pos_mat(it,4) = Gamma(t-1);
                    break
                end
            end
        end


        %% Organisation 

        Gamma_deg = rad2deg(Gamma);
        Gamma = [Gamma Gamma_deg];

        Theta_deg = rad2deg(Theta);
        Theta = [Theta Theta_deg];
        Output_mat = [time Altitude Velocity Mass r g Accel Delta_Gamma Gamma Theta Ground_Pass TVC Grav_Gamma TVC_Mag TVC_Mag_Lim TVC_Act TVC_Flag];

        for ct = 1:height(Output_mat)
            if Output_mat(ct,1) == 0
                Output_mat = Output_mat(1:ct-1,:);
                time = time(1:ct-1,:);
                Altitude = Altitude(1:ct-1,:);
                Velocity = Velocity(1:ct-1,:);
                Mass = Mass(1:ct-1,:);
                r = r(1:ct-1,:);
                g = g(1:ct-1,:);
                Accel = Accel(1:ct-1,:);
                Delta_Gamma = Delta_Gamma(1:ct-1,:);
                Gamma = Gamma(1:ct-1,:);
                Theta = Theta(1:ct-1,:);
                Ground_Pass = Ground_Pass(1:ct-1,:);
                TVC = TVC(1:ct-1,:);
                Grav_Gamma = Grav_Gamma(1:ct-1,:);
                TVC_Mag = TVC_Mag(1:ct-1,:);
                TVC_Mag_Lim = TVC_Mag_Lim(1:ct-1,:);
                TVC_Act = TVC_Act(1:ct-1,:);
                TVC_Flag = TVC_Flag(1:ct-1,:);
                break
            end
        end
        Output_tab = array2table(Output_mat,'VariableNames',{'time','Altitude','Velocity','Mass','r','g','Acceleration','Delta_Gamma','Gamma','Gamma_deg','Theta','Theta_deg','Ground_Pass','TVC','Grav_Gamma','TVC_Mag','TVC_Mag_Lim','TVC_Act','TVC_Flag'});
        
%           if Gamma(end,1) ==0
            Tag = Results(cx,cc);
            switch Tag
                case 1
                    figure(1)
                    plot(time,Altitude,'--','Color','#A2142F')
                    figure(2)
                    plot(time,Velocity,'--','Color','#A2142F')
                    figure(3)
                    plot(time,Gamma(:,2),':','Color','#A2142F')
                    figure(4)
                    plot(Velocity,Altitude,'--','Color','#A2142F')
                    
                    figure(6)
                    subplot(2,1,1)
                    plot(time,Altitude,':','Color','#A2142F')
                    subplot(2,1,2)
                    plot(time,Velocity,':','Color','#A2142F')
                case 3
                    figure(1)
                    plot(time,Altitude,'--','Color','#77AC30','LineWidth',2)
                    figure(2)
                    plot(time,Velocity,'--','Color','#77AC30','LineWidth',2)
                    figure(3)
                    plot(time,Gamma(:,2),'--','Color','#77AC30','LineWidth',3)
                    figure(4)
                    plot(Velocity,Altitude,'--','Color','#77AC30','LineWidth',2)
                    figure(5)
                    plot(Altitude,Gamma(:,2),'--','Color','#77AC30','LineWidth',2)
                    figure(6)
                    subplot(2,1,1)
                    plot(time,Altitude,'Color','#77AC30','LineWidth',3)
                    subplot(2,1,2)
                    plot(time,Velocity,'Color','#77AC30','LineWidth',3)
                    figure(7)
                    plot(time,Delta_Gamma,'--','Color','#77AC30','LineWidth',3)
                    
                    figure
                    subplot(2,1,1)
                    plot(time,Delta_Gamma)
                    hold on
                    plot(time,Grav_Gamma)
                    plot(time,TVC_Mag)
                    plot(time,TVC_Mag_Lim)
                    plot(time,TVC_Act)
                    grid on
                    legend('Delta_Gamma','Grav_Gamma','TVC_Mag','TVC_Mag_Lim','TVC_Act')
                    subplot(2,1,2)
                    plot(time,TVC)
                    hold on
                    plot(time,TVC_Flag)
                    legend('Flag','Real Flag')
                    
                    
                    edges = deg2rad([0 0.5 1.5 5 6]);
                    disctab(cx,:) = histcounts(TVC_Act,edges);
                    
                    S1TVC_Act = TVC_Act(1:tsep);
                    S2TVC_Act = TVC_Act(tsep+1:end);
                    
                    S1disctab(cx,:) = histcounts(S1TVC_Act,edges);
                    S2disctab(cx,:) = histcounts(S2TVC_Act,edges);
                    
                    
                    
                    
                case 2
                    figure(1)
                    plot(time,Altitude,'--','Color','#0072BD')
                    figure(2)
                    plot(time,Velocity,'--','Color','#0072BD')
                    figure(3)
                    plot(time,Gamma(:,2),':','Color','#0072BD')
                    figure(4)
                    plot(Velocity,Altitude,'--','Color','#0072BD')
                    figure(6)
                    suplot(2,1,1)
                    plot(time,Altitude,':','Color','#0072BD')
                    subplot(2,1,2)
                    plot(time,Velocity,':','Color','#0072BD')
            end 
%           end
    end
end
%% Analytics
% figure
% subplot(5,1,1)
% plot(time,Altitude)
% grid on
% xlabel('Time (s)')
% ylabel('Altitude (m)')
% xline(MECO1,'--')
% xline(S2_Ignition,'--')
% yline(343000,'--')
% subplot(5,1,2)
% plot(time,Velocity)
% grid on
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% xline(MECO1,'--')
% xline(S2_Ignition,'--')
% yline(3387,'--')
% subplot(5,1,3)
% plot(time,Mass)
% grid on
% xlabel('Time (s)')
% ylabel('Mass (kg)')
% xline(MECO1,'--')
% xline(S2_Ignition,'--')
% subplot(5,1,4)
% plot(time,Accel)
% grid on
% xlabel('Time (s)')
% ylabel('Acceleration (g''s)')
% xline(MECO1,'--')
% xline(S2_Ignition,'--')
% subplot(5,1,5)
% plot(time,Ground_Pass)
% grid on
% xlabel('Time (s)')
% ylabel('Ground Track (m)')
% xline(MECO1,'--')
% yline(S2_Ignition,'--')
% 
% figure
% plot(time,Gamma(:,2))
% grid on
% xlabel('Time (s)')
% ylabel('Gamma (deg)')
% yyaxis right
% plot(time,Theta(:,2))
% ylabel('Theta (deg)')
% 
figure
angle = 0:0.01:pi*2;
surfrad = 3396200;
xp = surfrad*cos(angle);
yp = surfrad*sin(angle);

plot(xp,yp)
grid on
hold on
plot(r,Theta)
% 
% figure
% plot(Ground_Pass,Altitude)
% hold on
% xlabel('Ground Track Distance(m)')
% ylabel('Altitude (m)')
% axis equal
% grid on
% yline(343000,'--')
% 
% figure
% plot(time,Gamma(:,2))
% grid on
% xlabel('time (s)')
% ylabel('Gamma (deg)')
% 
% figure
% plot(time,Delta_Gamma(:,1))
% grid on
% xlabel('time(s)')
% ylabel('Delta_Gamma (rad)')

figure
% scatter(Pos_mat(:,3),Pos_mat(:,2),[],Pos_mat(:,1))
scatter(Pos_mat(:,3),Pos_mat(:,2))
grid on
xlabel('Velocity (m/s)')
ylabel('Altitude (m)')
yline(342000,'--');
yline(344000,'--');
xline(3377,'--');
xline(3397,'--');
title('Final Orbit Parameters')
% colorbar

idx = Pos_mat(:,4) == 0;
Pos_mat_Gamma = Pos_mat(idx,:);

figure
% scatter(Pos_mat(:,3),Pos_mat(:,2),[],Pos_mat(:,1))
scatter(Pos_mat_Gamma(:,3),Pos_mat_Gamma(:,2))
grid on
xlabel('Velocity (m/s)')
ylabel('Altitude (m)')
yline(342000,'--');
yline(344000,'--');
xline(3377,'--');
xline(3397,'--');
title('Final Orbit Parameters')
