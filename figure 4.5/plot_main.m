clc
clear
close all;
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
load('SNR_data_MIMO_SDMA_RSMA.mat')
SNR_dB = 5:5:40;
Avg_RSMA_WMMSE = Avg_RSMA_max;
Avg_SDMA_WMMSE = Avg_SDMA;
load('Low_complexity_sumrate.mat')


%% low complexity plot with optimised result
figure()
plot(SNR_dB,Avg_RSMA_1,"--x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB,Avg_RSMA_2,"--o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB,Avg_RSMA_max,"-*",'LineWidth',1,'MarkerSize',8,'Color',"#EDB120"); % select the highest rate among all streams
hold on;
plot(SNR_dB,Avg_SDMA,"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
hold on;
plot(SNR_dB,Avg_RSMA_WMMSE,'-diamond','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB,Avg_SDMA_WMMSE,'--*','LineWidth',1,'MarkerSize',8)
grid on;
xlabel("SNR [dB]")
ylabel("Sum rate [bps/Hz]");
legend("RSMA SVD-BD $Q_{c} = 1$ with proposed","RSMA SVD-BD $Q_{c} = 2$ with proposed",'RSMA SVD-BD Adaptive $Q_{c}$ with proposed','SDMA BD low complexity','RSMA SVD-BD/MRT $Q_{c}=2$ WMMSE','SDMA MRT/BD WMMSE','Location','northwest');

% small figure in plot
axes("Position",[0.65,0.15,0.2,0.2]) % location of the small figure and size
box on; % show the box
index_x = 5:8;
% area where to magnify
plot(SNR_dB(index_x),Avg_RSMA_1(index_x),"--x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB(index_x),Avg_RSMA_2(index_x),"--o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB(index_x),Avg_RSMA_max(index_x),"-*",'LineWidth',1,'MarkerSize',8,'Color',"#EDB120"); % select the highest rate among all streams
hold on;
plot(SNR_dB(index_x),Avg_SDMA(index_x),"--v",'LineWidth',1,'MarkerSize',8,'Color',"#7E2F8E"); % SDMA BD rate
hold on;
plot(SNR_dB(index_x),Avg_RSMA_WMMSE(index_x),'-diamond','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB(index_x),Avg_SDMA_WMMSE(index_x),'--*','LineWidth',1,'MarkerSize',8);
axis tight;
grid on;

