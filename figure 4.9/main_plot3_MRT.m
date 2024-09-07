clc
clear
close all;
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
load("SNR_result_with_SVD_MRT_initilisation.mat")
Avg_RSMA_WMMSE = Avg_RSMA;
Avg_SDMA_WMMSE = Avg_SDMA;


load('SNR_RSMA_data_adaptive.mat')

figure()
plot(SNR_dB,Avg_R1_RSMA,"--x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB,Avg_R2_RSMA,"--o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB,Avg_selective_RSMA,"-*",'LineWidth',1,'MarkerSize',8,'Color',"#EDB120"); % select the highest rate among all streams
hold on;
plot(SNR_dB,Avg_SDMA_BD,"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
hold on;
plot(SNR_dB,Avg_RSMA_WMMSE,'-diamond','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB,Avg_SDMA_WMMSE,'--*','LineWidth',1,'MarkerSize',8)
grid on;
xlabel("SNR [dB]")
ylabel("Sum rate [bps/Hz]");
legend("RSMA SVD-BD with $Q_{c} = 1$","RSMA SVD-BD with $Q_{c} = 2$",'RSMA SVD-BD with Adaptive $Q_{c}$','SDMA BD low complexity','RSMA SVD-MRT WMMSE','SDMA MRT WMMSE','Location','northwest');

% small figure in plot
axes("Position",[0.65,0.15,0.2,0.2]) % location of the small figure and size
box on; % show the box
index_x = 5:8;
% area where to magnify
plot(SNR_dB(index_x),Avg_R1_RSMA(index_x),"--x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB(index_x),Avg_R2_RSMA(index_x),"--o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB(index_x),Avg_selective_RSMA(index_x),"-*",'LineWidth',1,'MarkerSize',8,'Color',"#EDB120"); % select the highest rate among all streams
hold on;
plot(SNR_dB(index_x),Avg_SDMA_BD(index_x),"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
hold on;
plot(SNR_dB(index_x),Avg_RSMA_WMMSE(index_x),'-diamond','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB(index_x),Avg_SDMA_WMMSE(index_x),'--*','LineWidth',1,'MarkerSize',8);
axis tight;
grid on;