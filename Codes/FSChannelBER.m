% Tue 09/08/2018
% Scr-F12
% analytical BER over fog/smoke chanel
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER over fog/smoke channel


%% BER calculation
if((Turb_EN == false) && (PE_EN == false) && (FS_EN == true))
    BER_Anl_fs = qfunc(sqrt(SNR_Anl*Atm_Att_Sim^2));  % analytical BER over fog/smoke channel
end
