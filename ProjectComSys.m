% In The Name Of God
clc
clear 
warning off 
%% Transmitter
% Last digit of student number is 8

% Part 1
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 4000;                 % Number of bits to process (4*100,000)
OV_factor = 2;              % Oversampling factor
Rs = 4e7;
Fs = 2 * Rs;

% Generating bits
rng default                        % Use default random number generator
Sample_vector = randi([0 1],n,1);  % Generate vector of binary data

% Converting Bits To Symbols :
L_sample = length(Sample_vector);
reshaped_vector = reshape(Sample_vector,L_sample/k,k);   % Reshape data into binary k-tuples, k = log2(M)
symbols_vector = bi2de(reshaped_vector);                 % Convert to integers

% Drawing Plot for First 25 Symb. 
figure
stem(symbols_vector(1:25));
title('Symbols');
xlabel('Symbol Index');
ylabel('Integer Value');

% Modulating :

Qam_vector = qammod(symbols_vector,M,0);

% Drawing Plot
scatterplot(Qam_vector);
title('Modulated Symbols');
% part 2
% We'll upsample in next part.

% Part 3

% Difining Tx Filter parameters
alpha = 0.35;
sps = 4;             % Rate of sampling
Nd = 10;             % Nd

% Impulse Response of Tx
Tx_filter = rcosdesign(alpha , Nd , sps);

% Upsampling & Applying Tx Filter
Output_Tx_filter = upfirdn(Qam_vector , Tx_filter , sps );


%%%%%
% 
% % Correcting Delay
% Delay_Number1 = floor(abs(length(Output_Tx_filter)/2 - length(Qam_vector)*2))+1;
% 
% % Synchronizing Received Signal
% uuuu = Output_Tx_filter(Delay_Number1-1:length(Output_Tx_filter)-Delay_Number1-1) ;
% 
% scatterplot(uuuu);
% 

%%%%%%

% Drawing plot
scatterplot(Output_Tx_filter)
title('Modulated Signal After Tx Filter');

% Eye diagram
eyediagram(real(Output_Tx_filter),2*sps);
title('Eye Diagram Of The Output_Tx_filter')

% Plotting the Amplitude Response Of the Tx Filter :

fft_Tx = my_Complex_fft(Tx_filter);
fft_Output_Tx_filter = my_Complex_fft(Output_Tx_filter);

% Comparing Amplitude Response Of The Tx Filter & Transmitted Signal
L_Tx_fft = length(Tx_filter);
L_Out_Tx = length(Output_Tx_filter);
f = Fs*(-L_Tx_fft/2:L_Tx_fft/2-1)/L_Tx_fft;
f1 = Fs*(-L_Out_Tx/2:L_Out_Tx/2-1)/L_Out_Tx;
figure
subplot(2,1,1)
plot(f,abs(fft_Tx))
title('|Tx(f)|')
xlabel('f');
subplot(2,1,2)
plot(f1,abs(fft_Output_Tx_filter))
title('|Transmitted Signal|')
xlabel('f');



%% Transmission Through Channel

i = complex(0,1);

% Coefficients Of Impulse Response Of Channel:
C_coff = [-0.4863 - 0.4696*i , 0.3339 + 0.8264*i , 0.0278 - 0.5302*i , 0.2233 - 0.3868*i , 0.0221 + 0.3155*i] ;

% After Applying Channel Model:
Output_Channel = filter(C_coff,1,Output_Tx_filter) ;

% Drawing plot
scatterplot(Output_Channel)
title('Modulated Signal After Getting through Channel')

% Adding WGN to signal
SNR = 15 ;
Noisy_Output_Channel = awgn(Output_Channel,SNR,'measured') ;

% Drawing plot
scatterplot(Noisy_Output_Channel)
title('Modulated Signal After Adding Noise')


% Computing H(f) of Channel
Hf_Channel = my_Complex_fft(impz(C_coff , 1 , length(Noisy_Output_Channel)));

% Computing and plotting H(f) of Channel
L_Channel = length(Hf_Channel);
f = Fs*(-L_Channel/2:L_Channel/2-1)/L_Channel;

figure
subplot(2,1,1)
plot(f,abs(Hf_Channel))
title('|Hc(f)|')
xlabel('f');
subplot(2,1,2)
plot(f,angle(Hf_Channel))
title('<Hc(f)')
xlabel('f');

%% Receiver

 
%%
% LMS-algorithm

N = 1000;               % Number of samples
         
data = randi([0 1],1,N);        % Random signal

d = qammod(data,16);    % 16qam Modulated signal (desired/output)
r = filter(C_coff,1,data);          % Signal after passing through channel
x = awgn(r, SNR);              % Noisy Signal after channel (given/input)

%LMS-parameters

iterations = 1000; % Number of iterations
eta = 0.001;  % Step size
order = 20;    % Order of the equalizer
U = zeros(1,order); % Input frame
W = zeros(1,order); % Initial Weigths
% Algorithm

for k = 1 : iterations
    for n = 1 : 1000
        U(1,2:end) = U(1,1:end-1);  % Sliding window
        U(1,1) = x(n);   % Present Input
     
        y = (W)*U';             % Calculating output of LMS
        e = d(n) - y;           % Instantaneous error 
        W = W +  eta * e * U ;  % Weight update rule of LMS
        J(k,n) = e * e';        % Instantaneous square error
    end
end

%Calculation of performance parameters

MJ = mean(J,2);     % Mean square error
CS = freqz(C_coff);    % Channel Spectrum
NF = (0:length(CS)-1)./(length(CS));          % Normalized Frequencies
IMR = -10*log10(real(CS).^2 + imag(CS).^2);   % Inverse channel magnitude response (desired)
IPR = -imag(CS)./real(CS);                    % Inverse channel phase response (desired)

ES = freqz(W);        % Equalizer Spectrum
EMR = 10*log10(real(ES).^2 + imag(ES).^2);    % Equalizer magnitude response
EPR = imag(ES)./real(ES);                     % Equalizer phase response

%Plots

fsize=12;   % Font size of figure
lw=2;       % linewidth of plot
figure % MSE
plot(10*log10(MJ),'->k','linewidth',lw)
hg=legend('MSE','Location','Best');
grid minor
xlabel('Epochs iterations','FontSize',fsize);
ylabel('Mean squared error (dB)','FontSize',fsize);
title('Cost function','FontSize',2*fsize);
set(hg,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)

figure % Magnitude reponse
subplot(2,1,1)
plot(NF,IMR,'b','linewidth',lw)
hold on
plot(NF,EMR,'--r','linewidth',lw)
hg=legend('Inverse Channel','Equalizer','Location','Best');
grid minor
xlabel('Normalized Frequency','FontSize',fsize);
ylabel('Magnitude (dB)','FontSize',fsize);
title('Magnitude response','FontSize',2*fsize);
set(hg,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)

subplot(2,1,2)
plot(NF, IPR,'g','linewidth',lw)
hold on
plot(NF, EPR,'--b','linewidth',lw)
hg=legend('Inverse Channel','Equalizer','Location','Best');
grid minor
xlabel('Normalized Frequency','FontSize',fsize);
ylabel('Phase shift (rad)','FontSize',fsize);
title('Phase response','FontSize',2*fsize);
set(hg,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)
%%

% Channel Equalization
h_channel = impz(C_coff , 1 , length(Noisy_Output_Channel));
Output_equalizer = ifft( fft(Noisy_Output_Channel) ./ fft(h_channel) );

% Drawing Plot
scatterplot(Output_equalizer);
title('Modulated Signal After Equalization')

% Applying Rx Filter And Down-sampling
Output_Rx_Filter = upfirdn(Output_equalizer, Tx_filter , 1 , sps );

% Drawing Plot
scatterplot(Output_Rx_Filter);
title('Modulated Signal After Rx Filter')

% Eye- Diagram


OUTPUT = Output_Rx_Filter ;
Received_Demod_Signal = qamdemod(OUTPUT,M,0) ;

% Correcting Delay
Delay_Number = floor(abs(length(symbols_vector)/2 - length(Received_Demod_Signal)/2))+1;

% Synchronizing Received Signal
Received_Signal_Sync = Received_Demod_Signal(Delay_Number:length(Received_Demod_Signal)-21+Delay_Number) ;


% Drawing Plot Of First 50 Samples To Check Recieved Signal And Transmitted Signal
figure
stem(symbols_vector(1:50));
hold on
stem(Received_Signal_Sync(1:50));
title('Comparing First 50 Symbols Of Transmitted And Received Signals')
legend('Trasmitted','Received')


% As SNR rises the error bits decrease which is obvious since the power of
% the signal increases.
