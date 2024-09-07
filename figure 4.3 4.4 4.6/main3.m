clc
clear
close all
%% Extra path and function
addpath('function\')
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% System Parameter
Nt = 40;                             % Tx number
N_k = 4;                            % Rx number per user
Num_U = 2;                          % number of user
Nr = N_k*Num_U;                     % total number of Rx in the system    
SNR_dB = 5:5:40;                    % SNR in dB
SNR = 10.^(SNR_dB./10);             % snr in linear
Max = 100;                          % number of realisation
tolerance =1e-5;                    % tolerance for iteative algorithm
Rnn = eye(N_k);                     % noise covariance matrix
max_count = 1e3;                    % maximum number of loop for iterative algorithm

% predefine the variable size
weight = [1;1];
sigma = [1,1];

% generate theb channel pool
H_all = (randn(Nr,Nt,Max) + 1i* randn(Nr,Nt,Max))./sqrt(2);


All_t = zeros(Max,length(SNR_dB));
R_max_SDMA = zeros(Max,length(SNR_dB));
R_1_stream = zeros(Max,length(SNR_dB));
R_2_stream = zeros(Max,length(SNR_dB));
R_stream_max = zeros(Max,length(SNR_dB));
R_SDMA_optimal = zeros(Max,length(SNR_dB));
R_RSMA_optimal = zeros(Max,length(SNR_dB));

%%  Start the loop

for i_snr = 1:length(SNR)
    i_snr
     for i_max = 1:Max
        H = H_all(:,:,i_max);
        % BD precoder design
        [P_k]=BD_design(N_k,H); % do not use the combine
      
        %% optimal t by newton raphson method for private part
        [t11] = netwon_raphson_method2(H,P_k,Num_U,tolerance,SNR(i_snr),1);
        [t22] = netwon_raphson_method2(H,P_k,Num_U,tolerance,SNR(i_snr),2);
        [t33] = netwon_raphson_method2(H,P_k,Num_U,tolerance,SNR(i_snr),3);
        [t44] = netwon_raphson_method2(H,P_k,Num_U,tolerance,SNR(i_snr),4);
        
        %% using newton raphson to determine the number of common stream (min 1, max 2)

        % SDMA Rate
        [R_k_SDMA_BD,~]=SDMA_BD_ZF(SNR(i_snr),H,P_k,N_k);

        % RSMA Sum rate with different stream
        R_max_SDMA(i_max,i_snr) = sum(R_k_SDMA_BD);

        %  one stream
        [R_1_stream(i_max,i_snr),~,~]=adaptive_common_stream(H,1,P_k,SNR(i_snr),t11,Rnn,Num_U);
        % two streams
        [R_2_stream(i_max,i_snr),~,~]=adaptive_common_stream(H,2,P_k,SNR(i_snr),t22,Rnn,Num_U);
        % three streams
        [R_3_stream(i_max,i_snr),~,~]=adaptive_common_stream(H,3,P_k,SNR(i_snr),t33,Rnn,Num_U);
        % four streams
        [R_4_stream(i_max,i_snr),~,~]=adaptive_common_stream(H,4,P_k,SNR(i_snr),t44,Rnn,Num_U);

        % exhaust result
        [Optimal_R_1_stream(i_max,i_snr),~]=RSMA_exhaust_search1(H,Num_U,P_k,SNR(i_snr),Rnn,1);
        [Optimal_R_2_stream(i_max,i_snr),~]=RSMA_exhaust_search1(H,Num_U,P_k,SNR(i_snr),Rnn,2);
        [Optimal_R_3_stream(i_max,i_snr),~]=RSMA_exhaust_search1(H,Num_U,P_k,SNR(i_snr),Rnn,3);
        [Optimal_R_4_stream(i_max,i_snr),~]=RSMA_exhaust_search1(H,Num_U,P_k,SNR(i_snr),Rnn,4);
   

    end
end 
%% select the highest one sum rate
[R_max_RSMA]=selection_max(SNR,Max,R_1_stream,R_2_stream,R_3_stream,R_4_stream);

%%
Avg_RSMA_4 = mean(R_4_stream,1);
Avg_RSMA_3 = mean(R_3_stream,1);
Avg_RSMA_2 = mean(R_2_stream,1);
Avg_RSMA_1 = mean(R_1_stream,1);
Avg_RSMA_max = mean(R_max_RSMA,1);

Avg_RSMA_4_exhaust = mean(Optimal_R_4_stream,1);
Avg_RSMA_3_exhaust = mean(Optimal_R_3_stream,1);
Avg_RSMA_2_exhaust = mean(Optimal_R_2_stream,1);
Avg_RSMA_1_exhaust = mean(Optimal_R_1_stream,1);

Avg_SDMA = mean(R_max_SDMA,1);

%% figure()
% figure 43
figure()
plot(SNR_dB,Avg_RSMA_1,'--x','LineWidth',1,'MarkerSize',8,'Color',"#D95319");
hold on;
plot(SNR_dB,Avg_RSMA_1_exhaust,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB,Avg_RSMA_2,'--o','LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB,Avg_RSMA_2_exhaust,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#D95319");
hold on;
plot(SNR_dB,Avg_SDMA,'--v');
hold on;
legend("RSMA SVD-BD $Q_{c}=1$ with proposed","RSMA SVD-BD $Q_{c}=1$ with exhaust","RSMA SVD-BD $Q_{c}=2$ with proposed","RSMA SVD-BD $Q_{c}=2$ with exhaust","SDMA BD low complexity")
grid on;
xlabel("SNR [dB]");ylabel("Sum Rate [bps/Hz]")

axes("Position",[0.65,0.15,0.2,0.2]) % location of the small figure and size
box on; % show the box
% area where to magnify
index_x = 5:8;
plot(SNR_dB(index_x),Avg_RSMA_1(index_x),"-x",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD"); % 1 stream
hold on;
plot(SNR_dB(index_x),Avg_RSMA_2(index_x),"-o",'LineWidth',1,'MarkerSize',8,'Color',"#D95319"); % 2 streams
hold on;
plot(SNR_dB(index_x),Avg_RSMA_1_exhaust(index_x),"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#0072BD");
hold on;
plot(SNR_dB(index_x),Avg_RSMA_2_exhaust(index_x),"--diamond",'LineWidth',1,'MarkerSize',8,'Color',"#D95319");
axis tight;
grid on;

%% error
error_4 = (Avg_RSMA_4_exhaust - Avg_RSMA_4)./Avg_RSMA_4_exhaust;
error_3 = (Avg_RSMA_3_exhaust - Avg_RSMA_3)./Avg_RSMA_3_exhaust;
error_2 = (Avg_RSMA_2_exhaust - Avg_RSMA_2)./Avg_RSMA_2_exhaust;
error_1 = (Avg_RSMA_1_exhaust - Avg_RSMA_1)./Avg_RSMA_1_exhaust;

figure() 
% figure 4 4
plot(SNR_dB,error_1*100,'--o','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB,error_2*100,'--x','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB,error_3*100,'--*','LineWidth',1,'MarkerSize',8);
hold on;
plot(SNR_dB,error_4*100,'--v','LineWidth',1,'MarkerSize',8);
legend("$Q_{c}=1$","$Q_{c}=2$","$Q_{c}=3$","$Q_{c}=4$");
grid on;
xlabel("SNR [dB]");ylabel("Relative Error [$\%$]")

% figure 4.6 in the report
% figure()
% max_index = 1:100;
% plot(max_index,R_1_stream(max_index,4),'--x','LineWidth',1,'MarkerSize',8);
% hold on;
% plot(max_index,R_2_stream(max_index,4),'--*','LineWidth',1,'MarkerSize',8);
% hold on;
% plot(max_index,R_max_RSMA(max_index,4),'-o','LineWidth',1,'MarkerSize',8);
% xlabel("Channel Index");
% ylabel("Sum Rate [bps/Hz]");
% grid on;

% legend("RSMA low complexity $Q_{c} = 1$",'RSMA low complexity $Q_{c}=2$','RSMA low complexity Adaptive $Q_{c}$')
