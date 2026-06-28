clc
clear
close all;
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels

load("SNR_data_MIMO_SDMA_RSMA.mat")
x_index = 1:20;
SNR = 5:5:40;

for i_snr = 1:length(Avg_SDMA)
    for i_max = x_index

        common(i_max,:,i_snr) = diag(RSMA_C_Precoder_all(:,:,i_max,i_snr)'*RSMA_C_Precoder_all(:,:,i_max,i_snr));
        private(i_max,:,i_snr) = diag(RSMA_P_Precoder_all(:,:,i_max,i_snr)'*RSMA_P_Precoder_all(:,:,i_max,i_snr));
        private_SDMA(i_max,:) = diag(SDMA_Precoder_all(:,:,i_max,i_snr)'*SDMA_Precoder_all(:,:,i_max,i_snr));
    end
%     figure()
%         subplot(3,1,1)
%         bar(x_index,common(:,:,i_snr));
%         xlabel('Channel index number');
%         ylabel('Power [W]');
%         grid on;
%         title(['Power distribution on RSMA common stream ($Q_{c}=2$) precoder at ' num2str(SNR(i_snr)) ' dB SNR'])
%         legend('stream 1','stream 2')
%         
%         subplot(3,1,2)
%         bar(x_index,private(:,:,i_snr));
%         xlabel('Channel index number');
%         ylabel('Power [W]');
%         grid on;
%         legend('stream 1','stream 2','stream 3','stream 4')
%         title(['Power distribution on RSMA private stream precoder at ' num2str(SNR(i_snr)) ' dB SNR'])
%         
%         subplot(3,1,3)
%         bar(x_index,private_SDMA);
%         xlabel('Channel index number');
%         ylabel('Power [W]');
%         grid on;
%         legend('stream 1','stream 2','stream 3','stream 4')
%        title(['Power distribution on SDMA precoder at ' num2str(SNR(i_snr)) ' dB SNR'])
end

%%
% try to plot the power distribution of the user antenna for the private
% stream precoder design.
 figure()
 for i = 1:8
        subplot(4,2,i)
        bar(x_index,common(:,:,i));
        xlabel('Channel index number');
        ylabel('Power [W]');
        grid on;
        title([ num2str(SNR(i)) ' dB SNR'])
 end
legend('stream 1','stream 2')


 figure()
 for i = 1:8
        subplot(4,2,i)
        bar(x_index,private(:,:,i))
        xlabel('Channel index number');
        ylabel('Power [W]');
        grid on;
        title([num2str(SNR(i)) ' dB SNR'])
 end
legend('stream 1','stream 2','stream 3','stream 4')
