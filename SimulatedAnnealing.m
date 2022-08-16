%% Simulated Annealing
close all
clear all

figure(1)
hold on
xlabel('Separation Time (s)')
ylabel('Burn 2 Time (s)')
grid on
title('Simulated Annealing')

figure(2)
hold on
ylabel('Separation Time (s)')
xlabel('Iterations')
grid on
title('Seperation Time')

figure(3)
hold on
ylabel('Burn 2 Time (s)')
xlabel('Iterations')
grid
title('Burn 2 Time')

figure(4)
hold on
ylabel('Guess Track - Loop')
xlabel('Iterations')
grid on


%% Probability Setup

% k = 1.380649e-23; % Boltzmann Constant
k = 1;


%% Algorithm Setup


Iterations = 2000;
PO_angles=[0:1:30];

Guess_ct = 1;
Boundary = 5;
Flag_2 =0;

LB_Hard = 20;
UB_Hard = 700;


%% Tracking
n = 0;

Guess_Limit = 0;
CurrentGuessct = 0;
Rolled = 0;
Tried = 0;
Final_Solutions = zeros(1,4);

%% Algorithm Loop

x0 = 400;
for cx = 1:numel(PO_angles)
    Track = zeros(1,2);
    GuessTrack = zeros(1,2);
    Rolled_track = zeros(1,2);
     
     T = 100; 
    for iter = 1:Iterations
   

        % Create Guesses
        while Flag_2 == 0 
            lb = x0 - Boundary;
            if lb <10
                lb = 10;
            end
            ub = x0 + Boundary;
            if ub>UB_Hard
                ub = UB_Hard;
            end
            x1 = lb + (ub - lb).*rand(1,1);
            x1 = round(x1,1);
            [Flag_2,Burn2] = traj_sim_Flag_PO(x1, PO_angles(cx));
    %         [Flag_2,Burn2] = traj_sim_Flag(x1);
            Boundary = Boundary + 1;
            GuessTrack(iter,1:2) = [x1,Burn2];
            Guess_ct = Guess_ct+1;
            CurrentGuessct = CurrentGuessct +1;
            n = n+1;
            if n >100
                T = T*0.5;
                n = 0;
            end
            if CurrentGuessct >1000
                disp('Guess Limit reached')
                Guess_Limit = 1;
                break
            end
        end

        if Guess_Limit == 1
            break
        end

    %     [Flag_1,Burn1] = traj_sim_Flag_PO(x0);
        [Flag_1,Burn1] = traj_sim_Flag_PO(x0, PO_angles(cx));
        delta_E = Burn2-Burn1;
        


        Prob = exp(-delta_E/(k*T));
        Prob_Track(iter,1) = Prob;
        Temp_Track(iter,1) = T;

        Track(iter,1) = x1;
        Track(iter,2) = Burn2;
        Track(iter,3) = PO_angles(cx);
        Track(iter,4) = Flag_2;
        
        Tried = Tried + 1;
        delta_E_Track(iter,1) = -delta_E;

        if delta_E < 0
            x0 = x1;
        elseif rand(1)<Prob
            x0 = x1;
            Rolled = Rolled +1;
            Rolled_track(iter,1) = 1;
        end   
        
        Flag_2 = 0;
        Boundary = 5;
        X0_Track(iter,1) = x0;
        CurrentGuessct = 0;
        
        try
            for ci = 2:length(Rolled_track)
                Rolled_track(ci,2) = Rolled_track(ci-1,2) + Rolled_track(ci,1);
            end
        end
    end

val = min(Track(:,2));
idx = Track(:,2) == val;

Best_Solutions = Track(idx,:);
Final_Solutions = [Final_Solutions;Best_Solutions];
    
figure(1)
scatter(Track(:,1),Track(:,2))

figure(2)
plot(Track(:,1))

figure(3)
plot(Track(:,2))

figure(4)
plot(GuessTrack(:,1))
PO_angles(cx);

figure(5)
hold on
plot(Temp_Track)
grid on
title('Temperature Reduction')

figure(6) 
hold on
plot(Prob_Track)
axis([0 inf 0 1])

figure(7)
hold on
plot(Rolled_track(:,2))
xlabel('Iterations')
axis equal
grid on
title('Cumulative Probability Success')
% figure
% subplot(3,1,1)
% plot(delta_E_Track(:,1))
% subplot(3,1,2)
% plot(Temp_Track(:,1))
% subplot(3,1,3)
% plot(Prob_Track(:,1))


end

Final_Solutions = Final_Solutions(2:end,:);

val_fin = min(Final_Solutions(:,2));
idx_fin = Final_Solutions(:,2) == val_fin;
Final_Solutions = Final_Solutions(idx_fin,:);
Final_Solutions = unique(Final_Solutions,'rows');
names = {'Separation Time (s)','Burn 2 Time (s)','Pitch-Over Angle (deg)','Validity Flag'};
Final_Solutions = array2table(Final_Solutions,'VariableNames',names);