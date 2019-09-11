% Tue 08/08/2018
% Scr-F10
% analytical BER based on pointing errors model
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER based on
% pointing errors model. The expression is adopted from Pap-F07-P03 and
% based on Gauss–Laguerre quadrature series expansion described in 
% Appendix A from Thesis by Mojtaba Mansour Abadi.


%% BER calculation
if((Turb_EN == false) && (PE_EN == true))
    BER_Anl_PE = zeros(1, size(I0, 2));  % empty BER array
    [x_i, w_i] = GLQRGen(k_ord_GL, 1e-14);  % series expansion coefficients
    % SNR value loop
    for Index_I = 1:length(I0)
        A_par = eta*I0(Index_I)*A0_PE_v/(sqrt(2*N0));  % Q function constant parameter
        U_SISO = exp(-4*sig_j_PE^2*x_i/W2_eq_PE_v);  %  Q function varying parameter
        Par_SISO = qfunc(A_par*U_SISO).*besseli(0, s_PE*sqrt(2*x_i)/sig_j_PE);  % F(x)
        BER_Anl_PE(Index_I) = exp(-s_PE^2/(2*sig_j_PE^2))*sum(w_i.*Par_SISO);  % calculate BER
    end
end