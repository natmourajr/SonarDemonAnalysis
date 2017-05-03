
function [y0,y1,y2,Fmax,fs2,nfft2,novr3]=reamostra1(y0,y1,y2)

filtro0='pf1k_2k.biq';	% Filtro para banda do sinal
filtro1='1de25_3.biq';	% Primeiro filtro de decimacao
filtro2='pb700.biq';		% Segundo filtro de decimacao

fs=31250;
novr1=25; %Valor da decimação a ser realizada nos dados.
nfft2=2048/2;
Fmax=1500/60; %Faixa de frequencia a ser apresentada na tela
fs2=2*Fmax;
novr3=floor(nfft2-2*Fmax*0.5);	% Calcula overlap para calculo da FFT
fs1=fs; %Frequencia de amostragem do sinal original
novr2=round((fs1/novr1)/(2*Fmax));	% Calcula segunda decimacao.
[b0,a0]=le_filtro(filtro0);	% Le filtro passa banda
[b1,a1]=le_filtro(filtro1);	% Le filtro primeira decimação
[b2,a2]=le_filtro(filtro2);	% Le filtro segunda decimação

%y0=filter(b0,a0,y0);			% Aplica filtro de seleção da banda do sinal
%y1=filter(b0,a0,y1);			% Aplica filtro de seleção da banda do sinal
        
y0=abs(y0);%y=y.^2;			% Demodula sinal
y1=abs(y1);%y=y.^2;			% Demodula sinal
y2=abs(y2);%y=y.^2;			% Demodula sinal
        
%b=ones(1,novr1)/novr1;a=1;	% Media movel (Procedimento do Bossan no MAGS)
y0=filter(b1,a1,y0);			% Aplica filtro para primeira decimacao
y0=y0(novr1:novr1:end);		% Decima sinal
y1=filter(b1,a1,y1);			% Aplica filtro para primeira decimacao
y1=y1(novr1:novr1:end);		% Decima sinal
y2=filter(b1,a1,y2);			% Aplica filtro para primeira decimacao
y2=y2(novr1:novr1:end);		% Decima sinal

y0=filter(b2,a2,y0);		% Aplica filtro para primeira decimacao
y0=y0(novr2:novr2:end);	% Decima sinal
y1=filter(b2,a2,y1);		% Aplica filtro para primeira decimacao
y1=y1(novr2:novr2:end);	% Decima sinal
y2=filter(b2,a2,y2);		% Aplica filtro para primeira decimacao
y2=y2(novr2:novr2:end);	% Decima sinal


        
