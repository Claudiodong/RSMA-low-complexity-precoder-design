clc
clear
close all
%% Extra path and function
addpath('function\')
addpath("SDMA_MIMO_SNR_function\");
addpath('RSMA_MIMO_SNR_function\');
% Set the default for all text to use LaTeX interpreter
set(0, 'defaultTextInterpreter', 'latex');          % For text
set(0, 'defaultLegendInterpreter', 'latex');        % For legends
set(0, 'defaultAxesTickLabelInterpreter', 'latex'); % For tick labels
%% System Parameter
Nt = 4;                             % Tx number
N_k = 2;                            % Rx number per user
Num_U = 2;                          % number of user
Nr = N_k*Num_U;                     % total number of Rx in the system    
SNR_dB = 5:5:40;                    % SNR in dB
SNR = 10.^(SNR_dB./10);             % snr in linear
Max = 1;                          % number of realisation
tolerance =1e-5;                    % tolerance for iteative algorithm
Rnn = eye(N_k);                     % noise covariance matrix
Q_c = 2;                            % number of common stream (initilise)
max_count = 1e3;                    % maximum number of loop for iterative algorithm
alpha = 0.2;
weight = [1;1];
sigma = [1,1];

% channel
H_all = (randn(N_k*Num_U,Nt,Max) + 1i* randn(N_k*Num_U,Nt,Max))./sqrt(2);

% predefine the variable
All_t = zeros(Max,length(SNR_dB));
R_max_SDMA = zeros(Max,length(SNR_dB));
R_1_stream = zeros(Max,length(SNR_dB));
R_2_stream = zeros(Max,length(SNR_dB));
R_stream_max = zeros(Max,length(SNR_dB));
R_SDMA_optimal = zeros(Max,length(SNR_dB));
R_RSMA_optimal = zeros(Max,length(SNR_dB));


%%  Start the loop

for i_snr = 1:length(SNR)
     R = zeros(1,Max);
     for i_max = 1:Max
        H = H_all(:,:,i_max);

        % BD precoder design
        [P_k,~]=BD_design(N_k,H); % do not use the combine
      
        %% optimal t by newton raphson method for private part
        [t] = netwon_raphson_method(H,P_k,Num_U,tolerance);% irrevelant to SNR

        %% using newton raphson to determine the number of common stream (min 1, max 2)       
        % two streams
        [R_2_stream(i_max,i_snr),Pri_precoder,Com_precoder]=adaptive_common_stream(H,2,Nr,P_k,SNR(i_snr),t,Rnn,Num_U);
       

        % adaptive stream
        [R_stream_max(i_max,i_snr),~,~]=adaptive_common_stream(H,Q_c,Nr,P_k,SNR(i_snr),t,Rnn,Num_U);

        

        % SDMA
        [Rate_p_k,Pri_precoder_optimised_SDMA]=SDMA_MIMO_rate(H,N_k,max_count,tolerance,sigma,SNR(i_snr),weight);
        [Rate_p_k1,Pri_precoder_optimised_SDMA1]=SDMA_MIMO_rate1(H,N_k,max_count,tolerance,sigma,SNR(i_snr),weight,P_k);
        R_SDMA(i_max,i_snr) = sum(Rate_p_k);
        R_SDMA1 = sum(Rate_p_k1);

        % RSMA
        [R1_1,R2_1,Common_precoder_optimsied1,Pri_precoder_optimised1]=RSMA_MIMO_rate1(H,N_k,SNR(i_snr),sigma,weight,tolerance,max_count,Pri_precoder,Com_precoder);
        [R1,R2,Common_precoder_optimsied,Pri_precoder_optimised]=RSMA_MIMO_rate(H,N_k,alpha,SNR(i_snr),sigma,weight,tolerance,max_count);
        R_RSMA(i_max,i_snr)= R1+R2;
        R_RSMA1= R1_1+R2_1;


        % store all the precoder design information
        if (R_SDMA(i_max,i_snr) > R_SDMA1)
            R_SDMA_optimal(i_max,i_snr) = R_SDMA(i_max,i_snr);
            SDMA_Precoder_all(:,:,i_max,i_snr) = Pri_precoder_optimised_SDMA;
        else
            R_SDMA_optimal(i_max,i_snr) = R_SDMA1;
            SDMA_Precoder_all(:,:,i_max,i_snr) = Pri_precoder_optimised_SDMA1;
        end

        if ( R_RSMA(i_max,i_snr) >  R_RSMA1)
            R_RSMA_optimal(i_max,i_snr) = R_RSMA(i_max,i_snr);
            RSMA_C_Precoder_all(:,:,i_max,i_snr) = Common_precoder_optimsied;
            RSMA_P_Precoder_all(:,:,i_max,i_snr) = Pri_precoder_optimised;
        else
            R_RSMA_optimal(i_max,i_snr) = R_RSMA1;
            RSMA_C_Precoder_all(:,:,i_max,i_snr) = Common_precoder_optimsied1;
            RSMA_P_Precoder_all(:,:,i_max,i_snr) = Pri_precoder_optimised1;
        end
        
    end
end 

%% average rate
Avg_SDMA_BD = mean(R_SDMA_optimal,1);
Avg_RSMA_max = mean(R_RSMA_optimal,1);
Avg_RSMA_MRT = mean(R_RSMA,1);
Avg_SDMA_MRT = mean(R_SDMA,1);

figure()
plot(SNR_dB,Avg_RSMA_MRT,'--x','LineWidth',1,'MarkerSize',8,'Color',"#0072BD")
hold on;
plot(SNR_dB,Avg_RSMA_max,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#0072BD")
hold on;
plot(SNR_dB,Avg_SDMA_BD,'--diamond','LineWidth',1,'MarkerSize',8,'Color',"#D95319")
hold on;
plot(SNR_dB,Avg_SDMA_MRT,'--x','LineWidth',1,'MarkerSize',8,'Color',"#D95319")
legend("RSMA SVD-MRT initialisation",'RSMA SVD-BD/MRT initialisation','SDMA BD/MRT initialisation','SDMA MRT initialisation');
grid on;
% title("Optimised precoder sum rate with different initilisation")
xlabel('SNR [dB]')
ylabel("Sum rate [bps/Hz]")



