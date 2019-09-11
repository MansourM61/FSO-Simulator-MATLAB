% Tue 06/08/2018
% Scr-F08
% calculating Q-factor
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate Q-factor from received signal


%% Q-factor calculation
V_1_Sim = mean(RX_Sim_OOK(TX_OOK > 0));  % Mean(V_1)
V_0_Sim = mean(RX_Sim_OOK(TX_OOK <= 0));  % Mean(V_0)
sig_1_Sim = std(RX_Sim_OOK(TX_OOK > 0));  % Std(V_1)
sig_0_Sim = std(RX_Sim_OOK(TX_OOK <= 0));  % Std(V_0)
Q_Fac = abs(V_1_Sim - V_0_Sim)/(sig_1_Sim + sig_0_Sim);  % Q-factor