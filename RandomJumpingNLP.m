%% Random Jumping Method
% Create random guesses, take minimum
clear all
close all
%% Sim Setup
PO_Angle = [0:80];
BestTrack = zeros(1,5);
for cx = 1:numel(PO_Angle)

%% First Stage Bounds
ub = 700;
lb = 300;
r = zeros(20,1);
Flag = zeros(size(r));
%% First Stage Guesses
while unique(Flag) ==0
    for ct = 1:numel(r) 
        r(ct,1) = rand;
        X(ct,1) = round(lb + r(ct)*(ub - lb),1);
        [Flag(ct,1),f(ct,1)] = traj_sim_Flag_PO(X(ct),PO_Angle(cx));
    end
end

Stage_1_Results = [r,X,f,Flag];
for ax = 1:length(Stage_1_Results)
    Stage_1_Results(ax,5) = PO_Angle(cx);
end

%% Filter & Locate Minimum Value

idx = Stage_1_Results(:,4) ==1;
Filt_S1_results = Stage_1_Results(idx,:);

[val_min, idx_min] = min(Filt_S1_results(:,3));
Min_S1 = Filt_S1_results(idx_min,:);
%% Stage 2 Bounds
ub2 = Min_S1(2) + 5;
lb2 = Min_S1(2) - 5;

r2 = zeros(20,1);
Flag2 = zeros(size(r));

%% Stage 2 Loop

while unique(Flag2) == 0
    for ct = 1:numel(r) 
        r2(ct,1) = rand;
        X2(ct,1) = round(lb2 + r2(ct)*(ub2 - lb2),1);
        [Flag2(ct,1),f2(ct,1)] = traj_sim_Flag_PO(X2(ct),PO_Angle(cx));
    end
end

Stage_2_Results = [r2,X2,f2,Flag2];
for ax = 1:length(Stage_2_Results)
    Stage_2_Results(ax,5) = PO_Angle(cx);
end

idx2 = Stage_2_Results(:,4) ==1;
Filt_S2_results = Stage_2_Results(idx2,:);

[val_min2, idx_min2] = min(Filt_S2_results(:,3));
Min_S2 = Filt_S2_results(idx_min2,:);

%% Final Value

if Min_S2(3) < Min_S1(3)
    Best_Run = Min_S2;
else
    Best_Run = Min_S1;
end
i = 1;
BestTrack = [BestTrack;Best_Run];
i = i+1;
%% Plots

figure(1)
xlabel('Seperation Time (s)')
ylabel('Second Stage Burn Time (s)')
title('Random Jumping Non-Linear Programming')
grid on
hold on
sz = 60;

for ct = 1:length(Stage_1_Results)
    if Stage_1_Results(ct,4) ==1
        scatter(Stage_1_Results(ct,2),Stage_1_Results(ct,3),sz,[0.4660 0.6740 0.1880])
%     else
%         scatter(Stage_1_Results(ct,2),Stage_1_Results(ct,3),sz,[0.6350 0.0780 0.1840])
    end
end

for ct = 1:length(Stage_2_Results)
    if Stage_2_Results(ct,:) == Best_Run
        scatter(Stage_2_Results(ct,2),Stage_2_Results(ct,3),sz,'green','filled')
    elseif Stage_2_Results(ct,4) ==1
        scatter(Stage_2_Results(ct,2),Stage_2_Results(ct,3),sz,[0.4660 0.6740 0.1880])
%     else
%         scatter(Stage_2_Results(ct,2),Stage_2_Results(ct,3),sz,[0.6350 0.0780 0.1840])
    end
end


end
BestTrack = BestTrack(2:end,:);
ShortestBurn = min(BestTrack(:,3));
ShortestBurnIndex = BestTrack(:,3) == ShortestBurn;


FinalSolutions = BestTrack(ShortestBurnIndex,:);
names = {'Random Guess','Corresponding Separataion Time','Burn 2 Length','Flag','Pitch-Over Angle'};
FinalSolutions = array2table(FinalSolutions,'VariableNames',names);
