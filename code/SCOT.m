%% Smoothed Coherence Transform (SCOT)
% Input data must not be filtered

function [C,lag] = SCOT(data_n,data_m,freq,fs)

[Pxy,f] = cpsd(data_n,data_m,[],[],2*length(data_m)-1,fs,'onesided');
Pxx = pwelch(data_n,[],[],2*length(data_n)-1,'onesided');
Pyy = pwelch(data_m,[],[],2*length(data_m)-1,'onesided');

ind1 = find(f <= freq(1),1,'last'); ind2 = find(f >= freq(end),1,'first');

Y = Pxy(ind1:ind2)./sqrt(Pxx(ind1:ind2).*Pyy(ind1:ind2));

C = ifftshift(ifft(Y,2*length(data_m))); % SCOT
lag = linspace(-1*length(C)/2,length(C)/2,length(C));