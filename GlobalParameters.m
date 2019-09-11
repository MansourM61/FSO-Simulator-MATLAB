% Tue 04/08/2018
% Scr-F07
% parameters for FSO system
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to define all parameters.


%% simulation
Sim_Par = 'POW';  % simulation parameter; 'POW' = average optical power is varid, 'EXR' = extention ratio is varied.
Avg_Opt_Pow_start = 0.1e-3;  % optical average power - start value (W)
Avg_Opt_Pow_end = 5e-3;  % optical average power - stop value (W)
Avg_Opt_Pow_step = 5;  % optical average power - step value (W)
Avg_Opt_Pow_const = 100e-3;  % onstant optical average power (W)
ext_ratio_start = 2;  % extinction ratio = P_max/P_min - start value
ext_ratio_end = 200;  % extinction ratio = P_max/P_min - stop value
ext_ratio_step = 10;  % extinction ratio = P_max/P_min - step value
ext_ratio_const = 20;  % constant extinction ratio = P_max/P_min
RepOrd = 2;  % repeatation Order


%% bit operations
Sync_EN = false;  % if synchronisation is needed
NoB = 5e6;  % no of bits for each burst transmission
BR = 1e6;  % Data rate (bps)
NoS = 5;  % No of samples per bit


%% baseband modulation
lambda = 1550e-9; % laser wavelength (m)
Sig_Step = 5e-3;  % signal step for each signal level (A)
BW = BR*1.25;  % bandwidth of baseband signal


%% laser
Las_Eff = 0.5;  % laser internal modulation efficiency (W/A)
Sig_Level = 0.5;  % a calibration constant - the DC level of optical signal for NRZ-OOK modulation


%% background light
P_background = 0.01e-3;  % background light power (W); radiation from diffused from sun and sky


%% geometrical loss
GL_EN = false;  % if the geometrical loss is taken into account
Prop_Mod = 'UNI';  % beam propagation/illumination model; 'UNI' = uniform propagation/illumination, 'GUS' = Gaussian propagation/illumination


%% lossy channel
MiscLoss = 0;  % miscellaneous channel loss (dB)
Link_Len = 500;  % link Length (m)
Fading_Add = 'M2';  % fading effect; 'M1' = fading affects average power, 'M2' = fading affects signal, 'M3' = fading affects average power + signal


%% fog/smoke channel
FS_EN = false;  % if the effect of fog/smoke is present
Vis_FS = 2.25;  % fog/smoke visibility (km)
lambda_0 = 550e-9;  % reference wavelength (green light) (m)
T_th_FS = 2/100;  % contrast threshold (typical value = 2%)


%% turbulence channel
Turb_EN = true;  % if the effect of turbulence is present
Cn2 = 4e-13;  % the refractive index structure coefficient (m^-2/3)
F_t = 500;  % turbulence maximum frequency (Hz); inverse of tepmoral coherence(1 - 10 mS)
Resamp_Turb = 'RECT';  % method of resampling the turbulence samples; 'RES' = uses 'resample' function, 'RECT' = uses 'rectpulse' function
Turb_Mod = 'GG';  % turbulence model; 'LN' = Log-Normal model, 'GG' = Gamma-Gamma model


%% pointing error channel
PE_EN = true;  % if the effect of pointing errors is present
sig_j_PE = 0.5;  % pointing errors horizontal jitter (m)
mu_h_PE = 0;  % horizontal displacement
mu_v_PE = 0;  % vertical displacement
F_p = 500;  % pointing errors maximum frequency (Hz); inverse of tepmoral coherence(1 - 10 mS)
Resamp_PE = 'RECT';  % method of resampling the pointing errors samples; 'RES' = uses 'resample' function, 'RECT' = uses 'rectpulse' function


%% transmitter optics
TX_AP_Type = 'ANG';  % transmitter beam parameter type; 'ANG' = divergence angle is given, 'DIA' = Tx aperture diameter is given
Prop_Type = 'GUS';  % laser propagation model; 'GUS' = Gaussain propagaion model, 'UNI' = uniform propagation model
theta_d_v = 10;  % full vertical divergence angle (Deg)
theta_d_h = 10;  % full horizontal divergence angle (Deg)
w_tx_v = 5e-3;  % vertical beam size (m)
w_tx_h = 5e-3;  % horizontal beam size (m)
Tx_Pos = [0; 0; 0];  % laser source cartesian position in the global coordinate (m, m, m)
Tx_Dir = [cosd(90); sind(90); 0];  % direction of the source propagation; cannot be 0
Tx_Ori = 0;  % orientation of the source (Deg); around local z axis


%% receiver optics
Rx_Ap_dia = 0.005;  % receiver aperture diameter (m)
Rx_Pos = [sqrt(Link_Len^2 - 2^2); 2; 0];  % receiver aperture cartesian position in the global coordinate (m, m, m)
Rx_Dir = [cosd(45 + 180); sind(45 + 180); 0];  % direction of the receiver aperture facing; normal to the detector face; cannot be 0
Rx_Ori = 0;  % orientation of the aperture (Deg); around local z axis
Rx_Trn = 85;  % receiver aperture transmittance (%)
AFOV = 1;  % receiver aperture full-angle angular field-of-view (Deg)


%% photodetector
PD_Resp = 0.5;  % responsivity of photodiode (A/W)
PD_Gain = 1e1;  % transimpedance amplifier gain (V/A)
PD_NEP = 1e-12;  % noise equivalent power (NEP) of optical receiver (W/sqrt(Hz))
PD_RL = 50;  % receiver load impedance (Ohms)
BW_BR_r = 1.25;  % bandwidth to bit rate ratio (Hz/bps)
q_ch = 1.60217662e-19;  % electron charge (C)


%% detection
Thresh_Len_Coeff = 1;  % this coefficient is used to shorten/widen the length of adaptive threshold estimation