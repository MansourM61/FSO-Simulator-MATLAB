% Tue 04/08/2018
% Scr-F07
% quantising the received signal
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to extract bits from received signal


%% quantising the signal
TX_OOK_bit = zeros(1, SigLength);  % fill transmit bit array with bit 0
TX_OOK_bit(TX_OOK > 0) = 1;  % change the bits to 1 based on threshold
clear TX_OOK;  % release the memory
% generating received bits matrix and threshold measurement
Threshold_Sim_Temp = MiscLoss*PD_Resp*PD_Gain*P_Sim_avg;  % fixed threshold based on the CSI
Threshold_Sim(1, Index) = Threshold_Sim(Index) + Threshold_Sim_Temp;  % update the threshold array
Index_Step = floor(length(RX_Sim_OOK)/(NoT*Thresh_Len_Coeff));  % size of bit chunck to perform adaptive thresholding
% adaptive thresholding loop
for Index_T = 1:(NoT*Thresh_Len_Coeff)
    Index_S = 1 + (Index_T - 1)*Index_Step;  % index of the lower boundary of bit chunck
    Index_E = Index_T*Index_Step;  % index of the higher boundary of bit chunck
    Threshold_Sim(Index_T + 1, Index) =...
        Threshold_Sim(Index_T + 1, Index) +...
        (max(RX_Sim_OOK(Index_S:Index_E)) +...
        min(RX_Sim_OOK(Index_S:Index_E)))/2;  % adaptive threshold calculation
end
RX_Sim_OOK_bit = zeros(1, SigLength);  % fill received bit array with bit 0
RX_Sim_OOK_bit(RX_Sim_OOK > Threshold_Sim_Temp) = 1;  % change the bits to 1 based on threshold
clear RX_Sim_OOK;  % release the memory