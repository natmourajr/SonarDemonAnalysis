%Calcula Espectrograma de um determinado sinal

function[xDemon,fDemon,tDemon]=calcula_espectro(valor,novr3,nfft2,fs2)
[Y,f,t] = spectrogram(valor-mean(valor),hanning(nfft2,'periodic'),novr3,nfft2,fs2);
Y=abs(Y);						% Modulo da FFT
ind=(1:8);Y(ind,:)=repmat(Y(length(ind),:),[length(ind) 1]); % Descarta 8 primeiros bins
Y=Y./tpsw(Y);% Normaliza usando TPSW
xDemon=Y';
fDemon=f*60;
tDemon=t;