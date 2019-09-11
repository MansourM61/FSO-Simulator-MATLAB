% Tue 04/08/2018
% Scr-F07
% analytical BER based on Log-Normal model
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER based on
% Log-Normal model. The expression is adopted from Pap-F07-P01


%% BER calculation
if((Turb_EN == true) && (PE_EN == false) && strcmp(Turb_Mod, 'LN'))
    BER_Anl_LN_turb = zeros(1, size(I0, 2));  % empty BER array
    [zi, wi] = GHQRGen(k_ord_GH, 1e-14);  % series expansion coefficients
    % SNR value loop
    for Index_I = 1:length(I0)
        U_SISO = eta*I0(Index_I)*exp(-2*sig_x_turb^2 + zi*sqrt(8*sig_x_turb^2))/sqrt(2*N0);  % calculate Q function argument
        Par_SISO = sum(wi.*qfunc(U_SISO));  %  caclulate the series expansion
        BER_Anl_LN_turb(Index_I) = Par_SISO/sqrt(pi);  % calculate BER
    end
end