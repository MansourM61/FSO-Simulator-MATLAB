% Tue 08/08/2018
% Scr-F10
% analytical BER based on Pointing Errors + Log-Normal model
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER based on
% Pointing Errors and Log-Normal model. The expression is adopted from Pap-F07-P03.


%% BER calculation
alpha = PHI;  % parameter assignment
beta = SAI;  % parameter assignment
if((PE_EN == true) && (Turb_EN == true) && strcmp(Turb_Mod, 'GG'))
    Par_1 = 2^(beta - 1)*gamma(beta/2 + 1/2)*exp(-s_PE^2/(2*sig_j_PE^2) +...
        -(s_PE^2*gam_PE^2/sig_j_PE^2)/(2*beta - 2*gam_PE^2))*...
        (alpha*beta/A0_PE_v)^beta;  % numinator section
    Par_2 = gamma(alpha)*gamma(beta)*sin((alpha - beta)*pi)*...
        gamma(-(alpha - beta) + 1)*...
        abs(gam_PE^2 - beta)*beta;  % denominator section
    Par_3 = gamma((beta + 1)/2)*sqrt(pi)*gam_PE^2;  % constant section
    Par_v = ((SNR_Anl)*2).^(-beta/2);  % varying part of BER
    BER_Anl_PE_GG_turb = (Par_1/Par_2)*Par_3*Par_v;  % calculate the asymptotic BER
end