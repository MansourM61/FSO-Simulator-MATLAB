% Tue 07/08/2018
% Scr-F08
% applying pointing errors effect on the signal
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to apply the effect of the pointing errors


%% pointing errors effect

A0_PE_v = 1;  % set the geometrical loss to 0 dB; it is considered in another way
Mean_anl = A0_PE_v*gam_PE^2/(1 + gam_PE^2)*exp(-s_PE^2/(2*(1 + gam_PE^2)*sig_j_PE^2));  % theoritical mean value of the intensity after pointing errors
Var_anl = A0_PE_v^2*gam_PE^2/(2 + gam_PE^2)*exp(-s_PE^2/((2 + gam_PE^2)*sig_j_PE^2)) -...
    A0_PE_v^2*gam_PE^4/(1 + gam_PE^2)^2*exp(-s_PE^2/((1 + gam_PE^2)*sig_j_PE^2));  % theoritical variance value of the intensity after pointing errors

U = randn(1, 2*NoPE);  % generate a normal random sequence
V = randn(1, 2*NoPE);  % generate a normal random sequence
r_h = sig_j_PE*U + mu_h_PE;  % horizontal jitter random sequence
r_v = sig_j_PE*V + mu_v_PE;  % vertical jitter random sequence
clear U;  % release the memory
clear V;  % release the memory
PE_Sim_0 = A0_PE_v*A0_PE_h*exp(-2*r_v.^2/W2_eq_PE_v).*exp(-2*r_h.^2/W2_eq_PE_h);  % pointing errors coefficient
clear r_h;  % release the memory
clear r_v;  % release the memory
if(strcmp(Resamp_PE, 'RES'))
    PE_Sim_1 = resample(PE_Sim_0, NoP, NoPE);  % resample generated numbers based on signal frequency to fading frequency ratio
elseif(strcmp(Resamp_PE, 'RECT'))
    PE_Sim_1 = rectpulse(PE_Sim_0, ceil(NoP/NoPE));  % resample generated numbers based on signal frequency to fading frequency ratio
else
    error('pick a resampling method');  % print out an error message and exit
end
clear PE_Sim_0;  % release the memory
PE_Sim_2 = PE_Sim_1(ceil(length(PE_Sim_1)/4) + (1:NoP));  % cut the extra numbers
clear PE_Sim_1;  % release the memory
if(min(PE_Sim_2) < 0)  % make sure all the coefficients are positive
    PE_Sim_Offset = -min(PE_Sim_2) + eps;  % set the offset value to proper level if minimum coeeficient is negative
else
    PE_Sim_Offset = eps;  % set the offset value to proper level if minimum coeeficient is negative
end
PE_Sim = PE_Sim_2 + PE_Sim_Offset;  % add the offset to turbulence coefficients
clear PE_Sim_2;  % release the memory

Mean_sim = mean(PE_Sim);  % simulated mean value of the intensity after pointing errors
Var_sim = var(PE_Sim);  % simulated variance value of the intensity after pointing errors

fprintf(['Pointing Errors:\t\tMean Value\t(Theory = %.3f,\t'...
    'Simulation = %.3f)\n'], Mean_anl, Mean_sim);  % print out the mean value
fprintf(['Pointing Errors:\t\tVar Value\t(Theory = %.3f,\t'...
    'Simulation = %.3f)\n'], Var_anl, Var_sim);  % print out the mean value