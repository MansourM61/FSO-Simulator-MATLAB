% Tue 04/08/2018
% Scr-F07
% analytical BER over clear channel
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER over clear channel


%% BER calculation
if((Turb_EN == false) && (PE_EN == false) && (FS_EN == false))
    BER_Anl_cl = qfunc(sqrt(SNR_Anl));  % analytical BER over clear channel
end