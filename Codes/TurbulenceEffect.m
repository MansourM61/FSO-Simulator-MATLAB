% Tue 04/08/2018
% Scr-F07
% applying turbulence effect on the signal
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to apply the effect of the turbulence


%% turbulence effect
if(strcmp(Turb_Mod, 'LN'))  % Log-Normal model

    Mean_anl = exp(mu_x_turb + 0.5*sig_x_turb^2);  % theoritical mean value of the intensity after Log-Normal turbulence
    Var_anl = 4*(exp(sig_x_turb^2) - 1)*exp(2*mu_x_turb + sig_x_turb^2);  % theoritical variance of the intesity after Log-Normal turbulence; the coefficient '4' comes from '2' in the exponent 'exp(2*y1)'
    
    % generating random coefficinets based on Log-Normal distribution
    U = randn(1, 2*NoT);  % generate a normal random sequence
    y1 = sig_x_turb*U + mu_x_turb;  % calculate exponent argument
    clear U;  % release the memory
    Turb_Sim_0 = exp(2*y1);  % calculate Log-Normal random sequence; THE COEFFICIENT 2 IS NEEDED, THE REASON sig1 = sqrt(sig2_R_Sim/4) DID NOT LEAD TO CORRECT VALUE.
    clear y1;  % release the memory
    if(strcmp(Resamp_Turb,'RES'))
        Turb_Sim_1 = resample(Turb_Sim_0, NoP, NoT);  % resample generated numbers based on signal frequency to fading frequency ratio
    elseif(strcmp(Resamp_Turb,'RECT'))
        Turb_Sim_1 = rectpulse(Turb_Sim_0, ceil(NoP/NoT));  % resample generated numbers based on signal frequency to fading frequency ratio
    else
        error('pick a resampling method');  % print out an error message and exit
    end
    clear Turb_Sim_0;  % release the memory
    Turb_Sim_2 = Turb_Sim_1(ceil(length(Turb_Sim_1)/4) + (1:NoP));  % cut the extra numbers
    clear Turb_Sim_1;  % release the memory
    if(min(Turb_Sim_2) < 0)  % make sure all the coefficients are positive
        Turb_Sim_Offset = -min(Turb_Sim_2) + eps;  % set the offset value to proper level if minimum coeeficient is negative
    else
        Turb_Sim_Offset = eps;  % set the offset value to proper level if minimum coeeficient is negative
    end
    Turb_Sim_3 = Turb_Sim_2 + Turb_Sim_Offset;  % add the offset to turbulence coefficients
    clear Turb_Sim_2;  % release the memory
    Dummy = sum(Turb_Sim_3);  % sum up turbulence coefficinets
    Turb_Sim = Turb_Sim_3/Dummy*NoP;  % normalize turbulence coefficients
    clear Turb_Sim_3;  % release the memory
    lognD_Sim = fitdist(Turb_Sim(1:(floor(NoT/2) + 1):end)', 'lognormal');  % estimate turbulence parameters
    sig2_I_anl = mean(Turb_Sim.^2)/mean(Turb_Sim)^2 - 1;  % simulated scintillation index for Log-Normal model
elseif(strcmp(Turb_Mod, 'GG'))  % Gamma-Gamma model

    Mean_anl = 1;  % theoritical mean value of the intensity after Log-Normal turbulence
    Var_anl = 1/PHI + 1/SAI + 1/(PHI*SAI);  % theoritical variance of the intesity after Log-Normal turbulence

    % generating random numbers based on Gamma-Gamma distribution
    D1 = gamrnd(PHI, 1, [1, 2*NoT]);  % generate a Gammma random sequence
    X1 = (1/PHI)*(D1);  % scale the distribution
    clear D1;  % release the memory
    D1 = gamrnd(SAI, 1, [1, 2*NoT]);  % generate a Gammma random sequence
    Y1 = (1/SAI)*(D1);  % scale the distribution
    clear D1;  % release the memory
    Turb_Sim_0 = X1.*Y1;  % calculate Gamma-Gammma random sequence
    clear X1 Y1;  % release the memory
    if(strcmp(Resamp_Turb,'RES'))
        Turb_Sim_1 = resample(Turb_Sim_0, NoP, NoT);  % resample generated numbers based on signal frequency to fading frequency ratio
    elseif(strcmp(Resamp_Turb,'RECT'))
        Turb_Sim_1 = rectpulse(Turb_Sim_0, ceil(NoP/NoT));  % resample generated numbers based on signal frequency to fading frequency ratio
    else
        error('pick a resampling method');  % print out an error message and exit
    end
    clear Turb_Sim_0;  % release the memory
    Turb_Sim_2 = Turb_Sim_1(ceil(length(Turb_Sim_1)/4) + (1:NoP));  % cut the extra numbers
    clear Turb_Sim_1;  % release the memory
    if(min(Turb_Sim_2) < 0)  % make sure all the coefficients are positive
        Turb_Sim_Offset = -min(Turb_Sim_2) + eps;  % set the offset value to proper level if minimum coeeficient is negative
    else
        Turb_Sim_Offset = eps;  % set the offset value to proper level if minimum coeeficient is negative
    end
    Turb_Sim_3 = Turb_Sim_2 + Turb_Sim_Offset;  % add the offset to turbulence coefficients
    clear Turb_Sim_2;  % release the memory
    Dummy = sum(Turb_Sim_3);  % sum up turbulence coefficinets
    Turb_Sim = Turb_Sim_3/Dummy*NoP;  % normalize turbulence coefficients
    clear Turb_Sim_3;  % release the memory
    sig2_I_anl = mean(Turb_Sim.^2)/mean(Turb_Sim)^2 - 1;  % simulated scintillation index for Gamma-Gamma model
else
    error('pick a turbulence model');  % print out an error message and exit
end

Mean_sim = mean(Turb_Sim);  % simulated mean value of the intensity after turbulence
Var_sim = var(Turb_Sim);  % simulated variance value of the intensity after turbulence
sig2_I_sim = mean(Turb_Sim.^2)/mean(Turb_Sim)^2 - 1;  % simulated scintillation index

fprintf('Turbulence:\tRytov Variance\t\t\t\t= %.3f\n', sig2_R_Sim);  % print out Rytov variance
fprintf('Turbulence:\tAperture Averaging Factor\t= %.3f\n', AAF);  % print out aperture averaging factor
fprintf(['Turbulence:\tMean Value\t\t\t(Theory = %.3f,\t'...
    'Simulation = %.3f)\n'], Mean_anl, Mean_sim);  % print out the mean value
fprintf(['Turbulence:\tVar Value\t\t\t(Theory = %.3f,\t'...
    'Simulation = %.3f)\n'], Var_anl, Var_sim);  % print out the mean value
fprintf(['Turbulence:\tScintillation Index\t(Theory = %.3f,\t'...
    'Simulation = %.3f)\n'], sig2_I_anl, sig2_I_sim);  % print out the scintillation index