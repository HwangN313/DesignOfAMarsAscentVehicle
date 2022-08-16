function [Burn2] = traj_sim_obj(Stage_Sep_time)


%% Sim Conditions
Time_step = 0.1; %s
                % Time Step should be 1s or lower to ensure that the stage
                % sep time is an integer
%% Recording
Length = 2000/Time_step; % 2000 at 1s
output = zeros(Length,6);

%% Initial Launch Stage - Pre-pitch over



% Intial Conditions

GM = 4.282837e13; % m3s-2
mfuel_dot_s1 = 1.792 * Time_step;
V_exh_1 = 3214.38 * Time_step;


mfuel_dot_s2 = 0.6514 * Time_step;
V_exh_2 = 2893.95 * Time_step;

V_1 = 0;
Alt_1 = 0;
R_1 = 3396200;
Gamma_1 = deg2rad(90);
M_1 = 400;
g_1 = 3.71;
t = 0;

% Stage_Sep = Stage_Sep_time;

InitVert = 10/Time_step;

for ct = 1:InitVert
    V_2 = V_1 + (mfuel_dot_s1 * V_exh_1/M_1) - g_1;
    Alt_2 = Alt_1 + V_1;
    R_2 = R_1 + V_1;
    M_2 = M_1 - mfuel_dot_s1;
    g_2 = GM/(R_1^2);
    
    V_1 = V_2;
    Alt_1 = Alt_2;
    R_1 = R_2;
    M_1 = M_2;
    g_1 = g_2;
    
    Gamma_2 = Gamma_1;
    
    t = t+1;
    output(t,1) = t;
    output(t,2) = V_2;
    output(t,3) = Alt_2;
    output(t,4) = R_2;
    output(t,5) = M_2;
    output(t,6) = Gamma_2;

end

%% Pitch Over

% Deg - make sure to use rad for maths
% PO_Target = deg2rad(Pitch_Over_Angle);
PO_Target = deg2rad(30); %actual heading is 90 - target.
Max_PO_rate = deg2rad(5)*Time_step;
Gamma_1 = deg2rad(90);


while PO_Target > 0
    V_2 = V_1 + (mfuel_dot_s1 * V_exh_1/M_1) - (g_1*sin(Gamma_1));
    Alt_2 = Alt_1 + (V_1*sin(Gamma_1));
    R_2 = R_1 + (V_1*sin(Gamma_1));
    M_2 = M_1 - mfuel_dot_s1;
    g_2 = GM/(R_1^2);
    if PO_Target>Max_PO_rate
        Gamma_2 = Gamma_1 - Max_PO_rate;
        PO_Target = PO_Target - Max_PO_rate;
    else
        Gamma_2 = Gamma_1 - PO_Target;
        PO_Target = 0;
    end
    
    V_1 = V_2;
    Alt_1 = Alt_2;
    R_1 = R_2;
    M_1 = M_2;
    Gamma_1 = Gamma_2;
    g_1 = g_2;
    
    t = t+1;
    output(t,1) = t;
    output(t,2) = V_2;
    output(t,3) = Alt_2;
    output(t,4) = R_2;
    output(t,5) = M_2;
    output(t,6) = Gamma_2;
end

%% Main Burn 1 - Burn til tank empty

while M_2 > 198
    V_2 = V_1 + (mfuel_dot_s1 * V_exh_1/M_1) - (g_1*sin(Gamma_1));
    Alt_2 = Alt_1 + (V_1*sin(Gamma_1));
    R_2 = R_1 + (V_1*sin(Gamma_1));
    M_2 = M_1 - mfuel_dot_s1;
    g_2 = GM/(R_1^2);
%     if Gamma_2>deg2rad(40)
        Gamma_2 = Gamma_1 + (V_1*cos(Gamma_1)/R_1) - ((1/V_1)*g_1*cos(Gamma_1));
%     else 
%         Gamma_2 = deg2rad(40);
%     end
    
    
    V_1 = V_2;
    Alt_1 = Alt_2;
    R_1 = R_2;
    M_1 = M_2;
    Gamma_1 = Gamma_2;
    g_1 = g_2;
    
    t = t+1;
    output(t,1) = t;
    output(t,2) = V_2;
    output(t,3) = Alt_2;
    output(t,4) = R_2;
    output(t,5) = M_2;
    output(t,6) = Gamma_2;
    
end


%% Stage Seperation

Stage_Sep_time = Stage_Sep_time/Time_step;

for ct = 1:Stage_Sep_time
    V_2 = V_1 - (g_1*sin(Gamma_1));
    Alt_2 = Alt_1 + (V_1*sin(Gamma_1));
    R_2 = R_1 + (V_1*sin(Gamma_1));
    M_2 = M_1;
    g_2 = GM/(R_1^2);
    Gamma_2 = Gamma_1 + (V_1*cos(Gamma_1)/R_1) - ((1/V_1)*g_1*cos(Gamma_1));
    
    
    
    V_1 = V_2;
    Alt_1 = Alt_2;
    R_1 = R_2;
    M_1 = M_2;
    Gamma_1 = Gamma_2;
    g_1 = g_2;
    
    t = t+1;
    output(t,1) = t;
    output(t,2) = V_2;
    output(t,3) = Alt_2;
    output(t,4) = R_2;
    output(t,5) = M_2;
    output(t,6) = Gamma_2;
    
    
end

if Alt_2<0
    return 
end

M_1 = M_1 - 93.6;
B1PO_time = t; % Burn 1 + Pitch over time
%% Stage 2

while M_1 > 47.4
    V_2 = V_1 + (mfuel_dot_s2 * V_exh_2/M_1) - (g_1*sin(Gamma_1));
    Alt_2 = Alt_1 + (V_1*sin(Gamma_1));
    R_2 = R_1 + (V_1*sin(Gamma_1));
    M_2 = M_1 - mfuel_dot_s2;
    g_2 = GM/(R_1^2);
%     if Alt_2 > 300000 && Alt_2 < 325000
%         Gamma_2 = deg2rad(20);
%     elseif Alt_2 > 100000 && Alt_2 < 300000
%         Gamma_2 = deg2rad(60);
%     elseif Alt_2 > 325000 && Alt_2 < 340000
%         Gamma_2 = deg2rad(5);
%     else
        Gamma_2 = Gamma_1 + (V_1*cos(Gamma_1)/R_1) - ((1/V_1)*g_1*cos(Gamma_1));
%     end
    
    
    V_1 = V_2;
    Alt_1 = Alt_2;
    R_1 = R_2;
    M_1 = M_2;
    Gamma_1 = Gamma_2;
    g_1 = g_2;
    
    t = t+1;
    output(t,1) = t;
    output(t,2) = V_2;
    output(t,3) = Alt_2;
    output(t,4) = R_2;
    output(t,5) = M_2;
    output(t,6) = Gamma_2;
    
    if V_2>3377 && V_2 <3397 && Alt_2 >333000 && Alt_2 <353000 && Gamma_2 <5 && Gamma_2 > -5
        break
    end
    
    
    
end

%% Output Trim
for ct = 1:length(output)
    if output(ct,1) ==0
        output = output(1:ct-1,:);
        break
    end
end

output(:,7) = rad2deg(output(:,6));

Velocity = output(end,2);
Altitude = output(end,3);
Gamma = output(end,7);
Burn2 = (output(end,1) - B1PO_time)*Time_step;



% %% Plots
% figure
% plot(output(:,3),output(:,2))
% xlabel('Altitude')
% ylabel('Velocity')
% 
% figure
% plot(output(:,1),output(:,7))
% xlabel('Time')
% ylabel('Gamma')

end