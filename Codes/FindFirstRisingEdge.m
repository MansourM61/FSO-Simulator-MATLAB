% Tue 04/08/2018
% Scr-F07
% finding the first rise transition
% Writer: Mojtaba Mansour Abadi
% Description: This file is used to find the first low to high transision


%% finding the first rising edge
Index_RE = 1;  % initial the rise edge index
if (Sync_EN == true)  % if synchronasation is needed
    Signal_Temp = Laser_Sim_Out - mean(Laser_Sim_Out);  % remove DC level from the signal
    P2P_Sim = std(Signal_Temp)*2;  % calculate peak-to-peak of the signal
    % loop through all points
    for Index_R = 1:(NoP - 1)
        F1 = Signal_Temp(Index_R)/P2P_Sim;  % difference of bit(i) and mean value of sequence
        F2 = Signal_Temp(Index_R + 1)/P2P_Sim;  % difference of bit(i+1) and mean value of sequence
        if ((F2 > 0) && (F1 < 0))  % if bit(i) == 0 and bit(i+1) == 1
            Index_RE = Index_R + 1;  % store the index
            break;  % leave the loop
        end
    end
    
    fprintf('The first rising edge index = %d\n', Index_RE);  % print out the first rising edge index
    
    TX_OOK =...
        (Laser_Sim_Out((Index_RE + round(NoS/2)):NoS:end)  - mean(Laser_Sim_Out))/...
        DelP_Sim - mean(Laser_Sim_Out);  % generate the transmit synchronised bits matrix
else
    TX_OOK = Raw_Bit_Sim;   % generate the transmit synchronised bits matrix
end

clear Raw_Bit_Sim;  % relesse the memory
