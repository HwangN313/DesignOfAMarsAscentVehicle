%% Thrust Script

Burn_time = 120; %s
ISP = 3069.568;
m0 = 400;
MarsG = 3.71;
MarsGM = 0.042828e6; %km3/s2
m1fuel = 204.4;
m2fuel = 51.1;

Mass_Flow1 = m1fuel/Burn_time;


%Launch
Thrust = Mass_Flow1 * ISP;
Weight = MarsG * m0;
Accel_Launch = (Thrust - Weight)/m0;

%MECO1
fun = @(t) (9.3621/2)*t.^2;
displace = integral(fun,0,Burn_time)