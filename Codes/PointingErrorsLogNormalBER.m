% Tue 08/08/2018
% Scr-F10
% analytical BER based on Pointing Errors + Log-Normal model
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER based on
% Pointing Errors and Log-Normal model. The expression is adopted from Pap-F07-P03.


%% BER calculation
if((PE_EN == true) && (Turb_EN == true) && strcmp(Turb_Mod, 'LN'))
    u_a = (s_PE/sig_j_PE)^2 + 2*sig_x_turb^2*gam_PE^2 + 2*sig_x_turb^2*gam_PE^4;  % u_a axulliary parameter
    C = 2^(gam_PE^2 - 1)*gamma(gam_PE^2/2 + 1/2)*exp(u_a)/(sqrt(pi)*(A0_PE_v)^(gam_PE^2));  % constant coefficient of BER
    Par = ((SNR_Anl)*2).^(-gam_PE^2/2);  % varying part of BER
    BER_Anl_PE_LN_turb = C*Par;  % calculate the asymptotic BER
end
