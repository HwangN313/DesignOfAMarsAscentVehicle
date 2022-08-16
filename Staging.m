%% 2 Stage

clear all 
close all


%% Setup
Init_Mass           = 400; % kg, initial mass
S2_1_Frac           = [0.1:0.01:0.9]; % stage 1 fraction vector
% S2_1_Frac         = 0.3;      % stage 1 fraction test value
S2_2_Frac           = 1 - S2_1_Frac; % Stage 2 frac value
Payload_Mass        = 16;

%Inert_Mass          = 72.5; % kg, Inert mass not prop sys
Inert_Mass = 111;
Prop_System_Mass    = 7;    % kg, inert prop sys mass

Init_2Stage_Prop_Mass       = Init_Mass - Inert_Mass - (2*Prop_System_Mass) - Payload_Mass;

%% Staging

S2_ISP1                 = 3068.568;
S2_Prop_1_Mass          = Init_2Stage_Prop_Mass*S2_1_Frac;
S2_Final_1_Mass         = Init_Mass - S2_Prop_1_Mass;

S2_DV1                  = S2_ISP1*log(Init_Mass./S2_Final_1_Mass);


ISP2                    = S2_ISP1;
S2_Init_Mass_2          = S2_Final_1_Mass - Prop_System_Mass - (Inert_Mass*S2_1_Frac);
S2_Prop_Mass_2          = Init_2Stage_Prop_Mass*S2_2_Frac;
S2_Final_Mass_2         = S2_Init_Mass_2-S2_Prop_Mass_2;


S2_DV2                  = ISP2*log(S2_Init_Mass_2./S2_Final_Mass_2);

S2_DV                   = S2_DV1 + S2_DV2;

% Analysis
figure
plot(S2_1_Frac,S2_DV)
[val_max_S2_DV,idx_max_S2_DV]     = max(S2_DV);
Best_Frac                         = S2_1_Frac(idx_max_S2_DV)
grid on
xlabel('Stage 1 Fraction')
ylabel('Delta-V (m/s)')
title('2 Stage Delta-V')

%% 3 Stage

S3_1_Frac = [0.1:0.01:0.9]';
S3_2_Frac = [0.1:0.01:0.9];
S3_Frac_Mat = NaN(height(S3_1_Frac),width(S3_2_Frac));
for ct = 1:height(S3_1_Frac)
    for cx = 1:width(S3_2_Frac)
        Val = 1 - S3_1_Frac(ct) - S3_2_Frac(cx);
        if Val>0
            S3_Frac_Mat(ct,cx) = Val;
        end
    end
end

S3_DV_Mat = NaN(size(S3_Frac_Mat));

for ct = 1:height(S3_Frac_Mat)
    for cx = 1:width(S3_Frac_Mat)
        
        S3_Frac1 = S3_1_Frac(cx);
        S3_Frac2 = S3_2_Frac(ct);
        S3_Frac3 = S3_Frac_Mat(cx,ct);
        
        if isnan(S3_Frac1)||isnan(S3_Frac2)||isnan(S3_Frac3)
            continue
        end

        S3_Init_Prop_Mass = Init_Mass - Payload_Mass - Inert_Mass - 3*Prop_System_Mass;

        S3_ISP1 = 3068.58;
        S3_ISP2 = S3_ISP1;
        S3_ISP3 = S3_ISP1;

        S3_Final_1_Mass = Init_Mass - (S3_Init_Prop_Mass*S3_Frac1);
        S3_DV1 = S3_ISP1 * log(Init_Mass/S3_Final_1_Mass);

        S3_Init_2_Mass = S3_Final_1_Mass - Prop_System_Mass - (S3_Frac1*Inert_Mass);
        S3_Final_2_Mass = S3_Init_2_Mass - (S3_Init_Prop_Mass*S3_Frac2);
        S3_DV2 = S3_ISP2 * log(S3_Init_2_Mass/S3_Final_2_Mass);

        S3_Init_3_Mass = S3_Final_2_Mass - Prop_System_Mass - (S3_Frac2*Inert_Mass);
        S3_Final_3_Mass = S3_Init_3_Mass - (S3_Init_Prop_Mass*S3_Frac3);
        S3_DV3 = S3_ISP3*log(S3_Init_3_Mass/S3_Final_3_Mass);

        S3_DV_Mat(ct,cx) = S3_DV1 + S3_DV2 + S3_DV3;
    end
end

%% Analytics
figure
surf(S3_1_Frac,S3_2_Frac,S3_DV_Mat)
xlabel('First Stage Fraction')
ylabel('Second Stage Fraction')
zlabel('Delta-V (m/s)')
%[val3,idx3] = max(max(S3_DV_Mat)) % if matrix, then max(A) returns a row vector containing the maximum value of each column of A. 
MaxDV3 = max(max(S3_DV_Mat));
[f2_idx,f1_idx] = find(S3_DV_Mat==MaxDV3);
MaxFracs = [S3_1_Frac(f1_idx),S3_2_Frac(f2_idx),S3_Frac_Mat(f1_idx,f2_idx)]
title('3 Stage Delta-V Outputs')


