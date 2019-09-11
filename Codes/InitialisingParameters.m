% Tue 14/08/2018
% Scr-F16
% initialising the parameters
% Writer: Mojtaba Mansour Abadi
% Description: This code is used to perform initialisations and calcutaion
% of required parameter.


%% initialisation
F_s = BR*NoS;  % sampling frequency (Hz)
NoP = NoB*NoS;  % total no of points
NoT = ceil((F_t/F_s)*NoP);  % total no of turbulence samples
NoPE = ceil((F_p/F_s)*NoP);  % total no of pointing error samples

if(strcmp(TX_AP_Type, 'ANG'))  % angle of divergence is given
    theta_0_v = (theta_d_v*pi/180)/2;  % half vertical divergence angle (Rad)
    theta_0_h = (theta_d_h*pi/180)/2;  % half horizontal divergence angle (Rad)
    w_0_v = lambda/(pi*theta_0_v);  % aperture vertical radius (m)
    w_0_h = lambda/(pi*theta_0_h);  % aperture horizontal radius (m)
elseif(strcmp(TX_AP_Type, 'DIA'))  % aperture diameter is given
    w_0_v = w_tx_v/2;  % aperture vertical radius (m)
    w_0_h = w_tx_h/2;  % aperture horizontal radius (m)
else
    error('pick an aperture type for transmitter');  % print out an error message and exit
end
z_0_v = pi*w_0_v^2/lambda;  % vertical Rayleigh range (m)
z_0_h = pi*w_0_h^2/lambda;  % horizontal Rayleigh range (m)
W_Rx_v = w_0_v*sqrt(1 + (Link_Len/z_0_v)^2);  % vertical beam radius at receiver (m)
W_Rx_h = w_0_h*sqrt(1 + (Link_Len/z_0_h)^2);  % horiontal beam radius at receiver (m)

k_Sim = 2*pi/lambda;  % wave number (rad/m)
sig2_R_Sim = 1.23*Cn2*k_Sim^(7/6)*Link_Len^(11/6);  % Rytov variance
d_AAF = sqrt((k_Sim*Rx_Ap_dia^2)/(4*Link_Len));  % auxilary parameter for aperture averaging effect
AAF = (1 + 1.062*d_AAF^2)^(-7/6);  % aperture averagin factor

sig_x_turb = sqrt(AAF*sig2_R_Sim/4);  % sigma parameter weak turbulence
mu_x_turb = -AAF*sig2_R_Sim/4;  % mu parameter weak turbulence
sig2_I_LN_anl = AAF*sig2_R_Sim;  % analytical scintillation index for Log-Normal model

Mean_LN_Sim_ANL = exp(2*mu_x_turb + 2*sig_x_turb^2);  % theoritical mean value of the intensity after Log-Normal turbulence
Var_LN_Sim_ANL = (exp(4*sig_x_turb^2) - 1)*exp(4*mu_x_turb + 4*sig_x_turb^2);  % theoritical variance of the intesity after Log-Normal turbulence; cofficient 4 comes from exp(2*X)

slnX2 = 0.49*sig2_R_Sim/(1 + 0.653*d_AAF^2 + 1.11*sig2_R_Sim^(6/5))^(7/6);  % logarithmic large-scale scintillation
slnY2 = 0.51*sig2_R_Sim*(1 + 0.69*sig2_R_Sim^(6/5))^(-5/6)/(1 + 0.9*d_AAF^2 + 0.621*d_AAF^2*sig2_R_Sim^(6/5));  % logarithmic small-scale scintillation
PHI = 1/(exp(slnX2) - 1);  % large-scale scintillation
SAI = 1/(exp(slnY2) - 1);  % small-scale scintillation
sig2_I_GG_anl = 1/PHI + 1/SAI + 1/(PHI*SAI);  % analytical scintillation index for Gamma-Gamma model

Mean_GG_Sim_ANL = 1;  % theoritical mean value of the intensity after Gamma-Gamma turbulence
Var_GG_Sim_ANL = 1/PHI + 1/SAI + 1/(PHI*SAI);  % theoritical variance of the intesity after Gamma-Gamma turbulence

v_PE_v = sqrt(pi/2)*(Rx_Ap_dia/(2*W_Rx_v));  % aperture to vertical beam size ratio
W2_eq_PE_v = W_Rx_v^2*sqrt(pi)*erf(v_PE_v)/(2*v_PE_v*exp(-v_PE_v^2));  % equivalent vertical beam size at receiver side

v_PE_h = sqrt(pi/2)*(Rx_Ap_dia/(2*W_Rx_h));  % aperture to horizontal beam size ratio
W2_eq_PE_h = W_Rx_h^2*sqrt(pi)*erf(v_PE_h)/(2*v_PE_h*exp(-v_PE_h^2));  % equivalent horizontal beam size at receiver side

A0_PE_v = 1;  % Geometrical loss for vertical direction; The value should be (erf(v_PE_v))^2, but we consider the geometrical loss in another form
A0_PE_h = 1;  % Geometrical loss for horizontal direction; The value should be (erf(v_PE_h))^2, but we consider the geometrical loss in another form

gam_PE = sqrt(W2_eq_PE_v)/(2*sig_j_PE);  % pointing errors gamma ratio
s_PE = sqrt(mu_h_PE^2 + mu_v_PE^2);  % boresight displacement
Mean_PE_Sim_Anl = A0_PE_v*gam_PE^2/(1 + gam_PE^2)*exp(-s_PE^2/(2*(1 + gam_PE^2)*sig_j_PE^2));  % theoritical mean value of the intensity after pointing errors
Var_PE_Sim_Anl = A0_PE_v^2*gam_PE^2/(2 + gam_PE^2)*exp(-s_PE^2/((2 + gam_PE^2)*sig_j_PE^2)) -...
    A0_PE_v^2*gam_PE^4/(1 + gam_PE^2)^2*exp(-s_PE^2/((1 + gam_PE^2)*sig_j_PE^2));  % theoritical variance value of the intensity after pointing errors

MEAN_PE_LN_Sim_Anl = A0_PE_v*gam_PE^2/(1 + gam_PE^2)*exp(-s_PE^2/(2*(1 + gam_PE^2)*sig_j_PE^2));  % theoritical mean value of the intensity after Log-Normal and pointing errors
Var_PE_LN_Sim_Anl = A0_PE_v^2*gam_PE^2/(2 + gam_PE^2)*exp(4*sig_x_turb^2 - s_PE^2/((2 + gam_PE^2)*sig_j_PE^2)) -...
    A0_PE_v^2*gam_PE^4/(1 + gam_PE^2)^2*exp(-s_PE^2/((1 + gam_PE^2)*sig_j_PE^2));  % theoritical variance value of the intensity after Log-Normal and pointing errors

Mean_PE_GG_Sim_Anl = A0_PE_v*gam_PE^2/(1 + gam_PE^2)*exp(-s_PE^2/(2*(1 + gam_PE^2)*sig_j_PE^2));  % theoritical mean value of the intensity after Gamma-Gamma and pointing errors
Var_PE_GG_Sim_Anl = A0_PE_v^2*gam_PE^2/(2 + gam_PE^2)*exp(-s_PE^2/((2 + gam_PE^2)*sig_j_PE^2))*...
    (1 + 1/PHI + 1/SAI + 1/(PHI*SAI)) - ...
    A0_PE_v^2*gam_PE^4/(1 + gam_PE^2)^2*exp(-s_PE^2/((1 + gam_PE^2)*sig_j_PE^2));  % theoritical variance value of the intensity after Gamma-Gamma and pointing errors

if(strcmp(Sim_Par, 'POW'))  % check if the simulation parameter is average optical power
    Loop_Order = Avg_Opt_Pow_step;  % set the loop no of iteration
    Avg_Opt_Pow_Array = linspace(Avg_Opt_Pow_start, Avg_Opt_Pow_end, Avg_Opt_Pow_step);  % create average optical power array
    ext_ratio_Array = ones(1, Loop_Order)*ext_ratio_const;  % create extinction ratio array
elseif(strcmp(Sim_Par, 'EXR'))  % check if the simulation parameter is extinction ratio
    Loop_Order = ext_ratio_step;  % set the loop no of iteration
    Avg_Opt_Pow_Array =  ones(1, Loop_Order)*Avg_Opt_Pow_const;  % create average optical power array
    ext_ratio_Array = linspace(ext_ratio_start, ext_ratio_end, ext_ratio_step);  % create extinction ratio array
else
    error('pick a simulation type');  % print out an error message and exit
end

BW = BW_BR_r*BR;  % required bandwidth (Hz)
P_noise_NEP = PD_NEP*sqrt(BW);  % power of noise (W) due to NEP
P_noise_BGD = (PD_Gain^2)*(2*q_ch*PD_Resp*P_background*BW)/PD_RL;  % total power of noise (W) due to background noise