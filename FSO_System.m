% Mon 03/09/2018
%
% FSO system simulation in MATLAB
%
% Writer: Mojtaba Mansour Abadi
% Description: This code is simulating a free-space optical communication system.
% The simulation parameters are defined in 'GlobalParameters.m' file.
% The simulation can be done for various parameter sweeps.
% In the 'Main Loop' section you can follow the steps from generating
% random bits to extracting the bits from the received signal and
% calculating BER.
% Due to an unkown memory error, instead of using function in the code, I
% have used a mechanism simula to MACRO definition in C/C++ language.
% Therefore, whenever needed, a sperate file containing the codes to
% perform a task is called with no argument/return values. This way all the
% variables defined in this file will be available to the code in the
% seperate file.
%

clc;
clearvars;
close all;


%% global parameters
run('.\\GlobalParameters.m');  % load the simulation parameters


%% initialisation
run('.\\Codes\\InitialisingParameters.m');  % initialise the simulation parameters


SNR_dB_Sim = zeros(1, Loop_Order);  % SNR array (dB)
BER_Sim = zeros(1, Loop_Order);  % BER array
Q_Sim = zeros(1, Loop_Order);  % Q-factor array
NoE_Sim = zeros(1, Loop_Order);  % no of erronous bits at each simulation step
Threshold_Sim = zeros(NoT*Thresh_Len_Coeff + 1, Loop_Order);  % threshold array; including values for adaptive thresholding


%% Main Loop
% calculate BER for each given parameter
for Index = 1:Loop_Order
    
    
    % ---------------------------------------------------------------------
    % assign the simulation parameters
    Avg_Opt_Pow = Avg_Opt_Pow_Array(Index);  % set the average optical power
    ext_ratio = ext_ratio_Array(Index);  % set the extinction ratio
    
    % ---------------------------------------------------------------------
    % Section 2
    % laser modulation part 1
    % calculating power levels for each SNR value
    P_Sim_avg = Avg_Opt_Pow;   % calculate average optical power
    P_Sim_0 = 2*P_Sim_avg/(1 + ext_ratio);  % output power for bit 0
    P_Sim_1 = 2*P_Sim_avg/(1 + 1/ext_ratio);  % output power for bit 1
    DelP_Sim = P_Sim_1 - P_Sim_0;  % amplityde of optical signal
    
    
    % to Increase the accurary, repeat each BER step "RepOrd" times
    for Index_Rep = 1:RepOrd
        
        % print out information regarding current iteration
        fprintf(['Index = %d out of %d\nAverage optical power\t= %5.3f dBm,\t %5.3f mW',...
            '\nExtinction ratio\t\t= %5.2f\n'],...
            Index, Loop_Order, 10*log10(Avg_Opt_Pow) + 30, Avg_Opt_Pow*1e3, ext_ratio);  % print out the SNR index
        fprintf('Repetition\t\t\t\t= %d out of %d\n', Index_Rep, RepOrd);  % print out the repetition index
        fprintf([repmat('*\t', 1, 18), '\n']);  % print out the separater
        
        
        % -----------------------------------------------------------------
        % Section 1
        % Section 2
        % pseudorandom binary sequence (PRBS) generation and NRZ-OOK bits
        Raw_Bit_Sim = randi([0, 1], 1, NoB);  % generate random bits
        Signal_Sim_Out = rectpulse(Raw_Bit_Sim, NoS);  % resample bits
        
        
        % -----------------------------------------------------------------
        % Section 2
        % Section 3
        % laser modulation part 2
        % generating output optical signal based on generated OOK signal
        % and calculated power levels
        Laser_Sim_Out = DelP_Sim*(Signal_Sim_Out - Sig_Level);  % generate laser output (W)
        Sig_Sim_Power = var(Laser_Sim_Out);  % measure the signal power
        clear Signal_Sim_Out;  % relesse the memory
        
        
        % -----------------------------------------------------------------
        % finding the 1st rising edge of signal based on transmit signal
        % finding rising edge
        run('.\\Codes\\FindFirstRisingEdge.m');  % find the first rising edge
        
        
        % -----------------------------------------------------------------
        % Section 11
        % geometrical loss effect
        if (GL_EN == true)
            run('.\\Codes\\GeometricalLoss.m');  % apply the geometrical loss
        else
            Geo_Loss_Sim = 1;  % geometrical loss is ignored
        end
        
        
        % -----------------------------------------------------------------
        % Section 11
        % fog/smoke channel effect
        if (FS_EN == true)
            run('.\\Codes\\AtmosphericAttenuation.m');  % apply the atmospheric atteaution effect            
        else
            Atm_Att_Sim = 1;  % atmospheric attenuationm effet is ignored
        end
        
        
        % -----------------------------------------------------------------
        % Section 9
        % Section 10
        % turbulence channel effect
        if (Turb_EN == true)
            run('.\\Codes\\TurbulenceEffect.m');  % apply the turbulence effect
        else
            Turb_Sim = 1;  % turbulence effet is ignored
        end
        
        
        % -----------------------------------------------------------------
        % Section 9
        % Section 10
        % pointing error channel effect
        if (PE_EN == true)
            run('.\\Codes\\PointingErrosEffect.m');  % apply the pointing errors effect
        else
            PE_Sim = 1;  % turbulence effet is ignored
        end
        
        
        % -----------------------------------------------------------------
        % Section 4
        % FSO channel effect
        Misc_Loss = 10^(-MiscLoss/10);  % total FSO channel loss
        Fading_Sim = Misc_Loss*Geo_Loss_Sim*Atm_Att_Sim*Turb_Sim.*PE_Sim;  % channel fading coefficient
        Mean_Fading_Sim = mean(Fading_Sim);  % mean value of the simulated fading coefficient
        Const_Loss = Atm_Att_Sim*Geo_Loss_Sim*Misc_Loss;  % constant loss due to atmosphere, propagation, misc loss; this loss only uses onstant attenuation rather than varying fading
        Var_Fading_Sim = var(Fading_Sim);  % variance value of the simulated fading coefficient
        clear Turb_Sim;  % release the memory
        clear PE_Sim;  % release the memory
        if(strcmp(Fading_Add, 'M1'))
            Rec_Sim_Opt = Fading_Sim.*P_Sim_avg + Laser_Sim_Out;  % received optical signal after loss and turbulence effect (Method 1)
        elseif(strcmp(Fading_Add, 'M2'))
            Rec_Sim_Opt = Fading_Sim.*Laser_Sim_Out + P_Sim_avg;  % received optical signal after loss and turbulence effect (Method 2)
        elseif(strcmp(Fading_Add, 'M3'))
            Rec_Sim_Opt = Fading_Sim.*(Laser_Sim_Out + P_Sim_avg);  % received optical signal after loss and turbulence effect (Method 3)
        else
            error('pick a fading implementation method');  % print out an error message and exit
        end
        fprintf('Channel Coefficient:\tMean Value\t\t= %.3f\n',...
            Mean_Fading_Sim);  % print out the mean value
        fprintf('Channel Coefficient:\tVar Value\t\t= %.3f\n',...
            Var_Fading_Sim);  % print out the mean value
        fprintf('Channel Coefficient:\tConstant Loss\t= %.3f dB\n',...
            -10*log10(Const_Loss));  % print out the mean value
        clear Laser_Sim_Out  % release the memory
        clear Fading_Sim;  % release the memory
        
        
        % -----------------------------------------------------------------
        % Section 5
        % optical to electical coversion (photodetector)
        Rec_Sim_Sig = PD_Resp*PD_Gain*Rec_Sim_Opt;  % PD conversion
        clear Rec_Sim_Opt;  % release the memory
        
        
        % -----------------------------------------------------------------
        % Section 6
        % applying SNR to the received signal
        P_noise_SHN = (PD_Gain^2)*(2*q_ch*PD_Resp*Avg_Opt_Pow*Const_Loss*BW)/PD_RL;  % total power of noise (W) due to shot noise
        P_noise = P_noise_NEP + P_noise_SHN + P_noise_BGD;  % total power of noise (W)        
        Rec_Sim_Sig_Power = ((Const_Loss*PD_Resp*PD_Gain)^2*Sig_Sim_Power)/PD_RL;  % received signal power (W); normalised to 50 Ohms
        SNR = Rec_Sim_Sig_Power/P_noise;  % signal-to-noise ratio (SNR)
        SNR_dB_Sim(Index) = 10*log10(SNR);  % update SNR array with SNR (dB)
        AdditiveNoise_Sim = randn(1, NoP)*sqrt(P_noise*PD_RL);  % generating white Gaussian noise
        Det_Sim_Sig = Rec_Sim_Sig + AdditiveNoise_Sim;  % adding white Gaussian noise to the detected signal
        clear Rec_Sim_Sig;  % release the memory
        clear AdditiveNoise_Sim;  % release the memory
        
        
        % -----------------------------------------------------------------
        % performing signall processing on detected electical signal
        LPF_RX_Sim_1 = Det_Sim_Sig - mean(Det_Sim_Sig);  % remove DC level
        P2P_RX_Sim = std(LPF_RX_Sim_1)*2;  % calculate peak-to-peak of received signal
        LPF_RX_Sim_2 = LPF_RX_Sim_1/P2P_RX_Sim;  % normalise received signal
        clear Det_Sim_Sig;  % release the memory
        RX_Sim_OOK = LPF_RX_Sim_2((Index_RE + round(NoS/2)):NoS:end);  % sample and hold the detected signal
        clear LPF_RX_Sim_1;  % release the memory
        clear LPF_RX_Sim_2;  % release the memory
        
        
        % -----------------------------------------------------------------
        % Section 8
        % analysing detected signal
        % calculating Q-factor parameters
        SigLength = length(RX_Sim_OOK);  % sampled signal length
        run('.\\Codes\\QFactorExtractor.m');  % calculate Q-factor
        
        
        % -----------------------------------------------------------------
        % Section 7
        % thresholding and extract the bits
        % generating transmit bits matrix
        run('.\\Codes\\Quantisation.m');  % extract bits from the received signal
        
        
        % -----------------------------------------------------------------
        % Section 8
        % BER calculation
        [NoE_Sim_Temp, BER_Sim_Temp] = biterr(TX_OOK_bit, RX_Sim_OOK_bit);  % calculate BER and no of erronous bits
        clear TX_OOK_bit;  % release the memory
        clear RX_Sim_OOK_bit;  % release the memory
        BER_Sim(Index) = BER_Sim(Index) + BER_Sim_Temp;  % update the BER array for each iteration
        NoE_Sim(Index) = NoE_Sim(Index) + NoE_Sim_Temp;  % update the no of errors array for each iteration
        Q_Sim(Index) = Q_Sim(Index) + Q_Fac;  % update the Q factor array for each iteration
        
        
        % -----------------------------------------------------------------
        % signal detection presentaion for each iteration
        fprintf(['SNR\t\t\t= %4.2f dB,\tBER\t\t\t= %5.3g\n',....
            'Threshold\t= %5.3f,\tQ-factor\t= %5.3f\n\n'],...
            SNR_dB_Sim(Index), BER_Sim_Temp, Threshold_Sim_Temp, Q_Fac);  % print out the received signal properties for each iteration
        
    end
    
    % -----------------------------------------------------------------
    % BER calculation for each SNR
    BER_Sim(Index) = BER_Sim(Index)/RepOrd;  % update the BER array
    Threshold_Sim(Index) = Threshold_Sim(Index)/RepOrd;  % update the threshold array
    Q_Sim(Index) = Q_Sim(Index)/RepOrd;  % update the Q-factor array
    
    % -----------------------------------------------------------------
    % signal detection presentaion for each SNR
    fprintf('SNR\t\t\t\t\t= %4.2f\nAverage BER\t\t\t= %5.3g\nAverage Threshold\t= %5.3f\nAverage Q-factor\t= %5.3f\n',...
        SNR_dB_Sim(Index), BER_Sim(Index), Threshold_Sim(1, Index), Q_Sim(Index));  % print out the received signal properties for each SNR
    fprintf([repmat('-', 1, 70), '\n']);  % print out the separater
    
end


%% analytical initialisation
eta = 1;  % responsivity for theoritical BER
N0 = 1;  % noise spectral density for theoritical BER
k_ord_GH = 400;  % no of Gauss�Hermite quadrature series expansion for theoritical BER
k_ord_GL = 50;  % no of Gauss�Laguerre quadrature series expansion for theoritical BER
SNR_dB_Anl = linspace(SNR_dB_Sim(1), SNR_dB_Sim(end), 100);  % SNR (dB) array for theoritical BER
SNR_Anl = 10.^(SNR_dB_Anl/10);  % SNR array for theoritical BER
I0 = sqrt(N0*2*SNR_Anl)/eta;  % intensity for theoritical BER


%% theoritical BER calculation
% theoritical analysis - pointing errors
run('.\\Codes\\PointingErrorsBER.m');  % calculating analytical BER based on Log-Normal model
  
% theoritical analysis - weak turbulence
run('.\\Codes\\LogNormalBER.m');  % calculating analytical BER based on Log-Normal model

% theoritical analysis - strong turbulence
run('.\\Codes\\GammaGammaBER.m');  % calculating analytical BER based on Gamma-Gamma model

% theoritical analysis - pointing errors + weak turbulence
run('.\\Codes\\PointingErrorsLogNormalBER.m');  % calculating analytical BER based on Gamma-Gamma model

% theoritical analysis - pointing errors + strong turbulence
run('.\\Codes\\PointingErrorsGammaGammaBER.m');  % calculating analytical BER based on Gamma-Gamma model

% theoritical analysis - clear channel
run('.\\Codes\\ClearChannelBER.m');  % calculating analytical BER over clear channel

% theoritical analysis - fog/smoke channel
run('.\\Codes\\FSChannelBER.m');  % calculating analytical BER over clear channel


%% plotting Results
Style = {'', 'o', '*', 's', '^', 'h', 'x', '+', 'd', 'v', '<', '>', 'p'};

figure;
hold on;
box on;

Sim_Type = '';  % simulation type; used for changing the name of the output graph

if((Turb_EN == false) && (PE_EN == false) && (FS_EN == false))
    MarkerPlot(SNR_dB_Anl, BER_Anl_cl, 'b', '-', Style{1}, 10);
    str_par = 'clear channel';
    Sim_Type = '-Clear_Channel';  % adjust simulation type variable
elseif((Turb_EN == false) && (PE_EN == false) && (FS_EN == true))
    MarkerPlot(SNR_dB_Anl, BER_Anl_fs, 'b', '-', Style{1}, 10);
    str_par = 'fog/smoke channel';
    Sim_Type = '-Fog_Channel';  % adjust simulation type variable
elseif((Turb_EN == false) && (PE_EN == true))
    MarkerPlot(SNR_dB_Anl, BER_Anl_PE, 'b', '-', Style{1}, 10);
    str_par = sprintf('pointing error channel with \\sigma_{j} = %4.2f', sig_j_PE);
    Sim_Type = '-PE_Channel';  % adjust simulation type variable
elseif((Turb_EN == true) && (PE_EN == false) && strcmp(Turb_Mod, 'LN'))
    MarkerPlot(SNR_dB_Anl, BER_Anl_LN_turb, 'b', '-', Style{1}, 10);
    str_par = sprintf('turbulence channel (Log-Normal model) with \\sigma_{R}^2 = %4.2f', sig2_R_Sim);
    Sim_Type = '-LN_Channel';  % adjust simulation type variable
elseif((Turb_EN == true) && (PE_EN == false) && strcmp(Turb_Mod, 'GG'))
    MarkerPlot(SNR_dB_Anl, BER_Anl_GG_turb, 'b', '-', Style{1}, 10);
    str_par = sprintf('turbulence channel (Gamma-Gamma model) with \\sigma_{R}^2 = %4.2f', sig2_R_Sim);
    Sim_Type = '-GG_Channel';  % adjust simulation type variable
elseif((Turb_EN == true) && (PE_EN == true) && strcmp(Turb_Mod, 'LN'))
    MarkerPlot(SNR_dB_Anl, BER_Anl_PE_LN_turb, 'b', '-', Style{1}, 10);
    str_par = sprintf(['turbulence channel (Log-Normal model) with \\sigma_{R}^2 = %4.2f\n'...
        'and pointing error channel with \\sigma_{j} = %4.2f'], sig2_R_Sim, sig_j_PE);
    Sim_Type = '-LN_PE_Channel';  % adjust simulation type variable
elseif((Turb_EN == true) && (PE_EN == true) && strcmp(Turb_Mod, 'GG'))
    MarkerPlot(SNR_dB_Anl, BER_Anl_PE_GG_turb, 'b', '-', Style{1}, 10);
    str_par = sprintf(['turbulence channel (Gamma-Gamma model) with \\sigma_{R}^2 = %4.2f\n'...
        'and pointing error channel with \\sigma_{j} = %4.2f'], sig2_R_Sim, sig_j_PE);
    Sim_Type = '-GG_PE_Channel';  % adjust simulation type variable
end

MarkerPlot(SNR_dB_Sim, BER_Sim, 'r', '--', Style{2}, 10);
Dummy = BER_Sim(BER_Sim > 0);
axis([SNR_dB_Sim(1), SNR_dB_Sim(end), min(Dummy), 1]);
xlabel('SNR (dB)');
ylabel('BER');

str_title = sprintf('BER vs. SNR for %s', str_par);
title(str_title);

MakeitPretty(gcf, [10, 9], ['L', 'G'], [12, 1.5, 5, 10], ['.\\Graphs\\BER_SNR', Sim_Type]);