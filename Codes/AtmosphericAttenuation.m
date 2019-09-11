% Tue 09/08/2018
% Scr-F11
% fog/smoke attenuation
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculae fog/smoke attenuation.


%% atmospheric attenuation
switch logical(true)  % q value based on Kim model
    case (Vis_FS > 50)  % Visibility > 50 km
        q_FS = 1.6;
    case ((Vis_FS <= 50) && (Vis_FS > 6))  % 6 km < Visibility <= 50 km
        q_FS = 1.3;
    case ((Vis_FS <= 6) && (Vis_FS > 1))  % 1 km < Visibility <= 6 km
        q_FS = 0.67*Vis_FS + 0.34;
    case ((Vis_FS <= 1) && (Vis_FS > 0.5))  % 0.5 km < Visibility <= 1 km
        q_FS = Vis_FS + 0.5;
    case (Vis_FS <= 0.5)  % Visibility <= 0.5 km
        q_FS = 0;
    otherwise  % otherwise
        q_FS = 0;
end
beta_l = (-log(T_th_FS)/Vis_FS)*(lambda/lambda_0)^(-q_FS);  % beta_lambda
Atm_Att_Sim = exp(-beta_l*Link_Len/1e3);  % fog/smoke attenuation

fprintf(['Fog/smoke:\tVisibility\t= %.2f km\n'...
    'Fog/smoke:\tLoss\t\t= %.3f dB\n'], Vis_FS, -10*log10(Atm_Att_Sim));  % print out the loss