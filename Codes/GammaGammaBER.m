% Tue 04/08/2018
% Scr-F07
% analytical BER based on Gamma-Gamma model
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to calculate analytical BER based on
% Gamma-Gamma model. The expression is adopted from Pap-F07-P02.


%% BER calculation
if((Turb_EN == true) && (PE_EN == false) && strcmp(Turb_Mod, 'GG'))
    C = 2^(PHI + SAI - 3)/(sqrt(pi^3)*gamma(PHI)*gamma(SAI));  % constant parameter
    z = SNR_Anl*(2/(PHI*SAI))^2;  % meijer G function argument
    a1 = (1 - PHI)/2;  % meijer G function parameter
    a2 = (2 - PHI)/2;  % meijer G function parameter
    a3 = (1 - SAI)/2;  % meijer G function parameter
    a4 = (2 - SAI)/2;  % meijer G function parameter
    a5 = 1;  % meijer G function parameter
    b1 = 0;  % meijer G function parameter
    b2 = 1/2;  % meijer G function parameter
    Par_a = [a1, a2, a3, a4, a5];  % meijer G function parameter vector
    Par_b = [b1, b2];  % meijer G function parameter vector
    Str = 'meijerG(2, 4, [%f, %f, %f, %f, %f], [%f, %f], %f)';  % meijer G function engine string
    BER_Anl_GG_turb = zeros(size(z));  % BER array
    % SNR value loop
    for Index_I = 1:length(I0)
        MGF_EGC = evalin(symengine, sprintf(Str, Par_a, Par_b, 2*z(Index_I)));  % use engine to evaluate meijer G function; 
                                                                                % THE COEFFICIENT 2 IN THE ARGUMENT IS ADDED TO MATCH THEORY AND SIMULATION; 
                                                                                % In defferent literatures, the SNR value has a different definitions. The coefficient 2
                                                                                % is needed to compensate for the difference.
                                                                                % In Pap-F07-P01, SNR = (eta*I)^2/(2*N0)
                                                                                % In Pap-F07-P02, SNR = (eta*I)^2/(1*N0)
        BER_Anl_GG_turb(Index_I) = C*MGF_EGC;  % calculate BER
    end
end
