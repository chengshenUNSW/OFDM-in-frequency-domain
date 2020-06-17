% OFDM in frequency domain


clear all
nFFT        = 64; % fft size
nDSC        = 64; % number of data subcarriers
nBitPerSym  = 64; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym        = 10^4; % number of symbols

EbN0dB      = [0:50]; % bit to noise ratio
EsN0dB      = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio

for ii = 1:length(EbN0dB)

   % Transmitter
   ipBit = rand(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod = 2*ipBit-1; % BPSK modulation 0 --> -1, 1 --> +1
   ipMod = reshape(ipMod,nBitPerSym,nSym).'; % grouping into multiple symbolsa; time by frequency

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF = [zeros(nSym,(nFFT-nDSC)/2) ipMod(:,[1:nBitPerSym/2]) ipMod(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,(nFFT-nDSC)/2)] ;
   

   % multipath channel
   nTap = 1;
   ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + j*randn(nSym,nTap)); 
   % each taps are iid--is that the reality? 
   % each series of taps are independent--is that the reality?
   % But, for each OFDM symbol, the channel taps stays constant
   
   % computing and storing the frequency response of the channel, for use at recevier
   hF = fftshift(fft(ht,64,2));

   
   % Gaussian noise of unit variance, 0 mean
   %nt = 1/sqrt(2)*[randn(1,nSym*(80+nTap-1)) + j*randn(1,nSym*(80+nTap-1))];
   nt = 10^(-EsN0dB(ii)/20)*1/sqrt(2)*[randn(nSym,64) + j*randn(nSym,64)];
   %nF = fft(nt,64,2);
   nF = nt;
   
   yFF = xF.*hF + nF;
   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   %ytt = sqrt(80/64)*xt; 
   %yt = ytt + 10^(-EsN0dB(ii)/20)*nt;
  
   % equalization by the known channel frequency response
   yF = yFF./hF;

   % extracting the required data subcarriers
   yMod = yF(:,[(nFFT-nDSC)/2+[1:nBitPerSym/2] (nFFT-nDSC)/2+[nBitPerSym/2+1:nBitPerSym] ]); 

   % BPSK demodulation
   % +ve value --> 1, -ve value --> -1
   ipModHat = 2*floor(real(yMod/2)) + 1;
   ipModHat(find(ipModHat>1)) = +1;
   ipModHat(find(ipModHat<-1)) = -1;

   % converting modulated values into bits
   ipBitHat = (ipModHat+1)/2;
   ipBitHat = reshape(ipBitHat.',nBitPerSym*nSym,1).';

   % counting the errors
   nErr(ii) = size(find(ipBitHat - ipBit),2);

end

simBer = nErr/(nSym*nBitPerSym);
EbN0Lin = 10.^(EbN0dB/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));

close all; figure(2);
semilogy(EbN0dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(EbN0dB,simBer,'mx-','LineWidth',2);
axis([0 50 10^-6 1])
grid on
legend('Rayleigh-Theory', 'Rayleigh-Simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
