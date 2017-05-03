%Programa realiza a analise DEMON com os dados sinteticos em toda a janela
%temporal de 60 segundos com os dados adquiridos no CIAMA, usando dois contatos cavitando.

clear all
close all

scrsz = get(0,'ScreenSize');
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize',14);


addpath('functions/comuns');
addpath('functions/demon');

%Carregando as marca??es com o tamanho da janela temporal
load an0.mat;
y0_total=y;
load an1.mat;
y1_total=y;


tamanho_arq=min([length(y0_total) length(y1_total)]);

y0_total=y0_total(1:tamanho_arq);
y1_total=y1_total(1:tamanho_arq);



%Inicio da inicializacao das variaveis
janela=60; %Tamanho da janela temporal a ser analisada
fs=31250; %Frequ?ncia de amostragem do sistema
num_janela=floor(tamanho_arq/fs/janela); %N?mero de janelas a serem analisadas
inicio=1;
fim=janela;
valor=janela*num_janela; %Valor para realizar a parada das janelas
npontos=128; %N?mero de pontos da distribui??o
aux=1; %Vari?vel auxiliar, ponteiro dos elementos dos vetores
ss=0; %Vari?vel que seleciona se vai usar os algoritmos de ICA ou NMF.
      %%ss=0, usar os algoritmos de ICA. Se ss=1, usar os algoritmos de NMF.
r=2; %Valor que deve ser usado na separa??o NMF para determinar o n?mero de componentes.
novr1=25; %Valor da decima??o a ser realizada nos dados.
nfft2=2048/2;
Fmax=1500/60; %Faixa de frequencia a ser apresentada na tela
fs2=2*Fmax;
bin_inicial=1;
bin_final=120;
novr3=floor(nfft2-2*Fmax*0.5);	% Calcula overlap para calculo da FFT
fs1=fs; %Frequencia de amostragem do sinal original
novr2=round((fs1/novr1)/(2*Fmax));	% Calcula segunda decimacao.
filtro0='pf1k_2k.biq';	% Filtro para banda do sinal
filtro1='1de25_3.biq';	% Primeiro filtro de decimacao
filtro2='pb700.biq';		% Segundo filtro de decimacao
var_armazena=0;



%Leitura dos filtros para a realiza??o da an?lise DEMON
[b0,a0]=le_filtro(filtro0);	% Le filtro passa banda
[b1,a1]=le_filtro(filtro1);	% Le filtro primeira decima??o
[b2,a2]=le_filtro(filtro2);	% Le filtro segunda decima??o

opt.ncomp=3;
opt.approach='defl';
opt.nonlin='tanh';
opt.thresh=0.09;
opt.norm_h=0; %valor de referencia=0
opt.norm_w=1; %Valor de referencia=1
opt.alpha=0.99;
opt.niter=100;
opt.verbose=3;
opt.stabilization='on';
opt.sampleSize=0;
opt.finetune='on';
% BSS
alg = 'jade';


while(aux <= num_janela)
    if aux==1
        y0=y0_total(inicio:fs*fim);
        y1=y1_total(inicio:fs*fim);
    else
        y0=y0_total(inicio:fim);
        y1=y1_total(inicio:fim);
    end
    
    MMt1 = [y0;y1];
    [wt1,Stempo1,errst1]=m_nmf(MMt1,alg,opt);
    
    Matriz_tempo(aux,:)=y0;%mat_input(1,:);
    Matriz_tempo(aux,:)=y1;%mat_input(2,:);
    Matriz_comp_tempo1(aux,:)=Stempo1(1,:);%mat_comp(1,:);
    Matriz_comp_tempo1(aux,:)=Stempo1(2,:);%mat_comp(2,:);
    
    
    [y0,y1,Fmax,fs2,nfft2,novr3]=reamostra(y0,y1);
    MMt=[y0;y1];
    [wt,Stempo,errst]=m_nmf(MMt,alg,opt);
    Matriz_ream1(aux,:)=y0;%mat_input(1,:);
    Matriz_ream2(aux,:)=y1;%mat_input(2,:);
    Matriz_comp1(aux,:)=Stempo(1,:);%mat_comp(1,:);
    Matriz_comp2(aux,:)=Stempo(2,:);%mat_comp(2,:);
    [xDemon1,fDemon]=calcula_espectro(Matriz_ream1(aux,:),novr3,nfft2,fs2);
    [xDemon2]=calcula_espectro(Matriz_ream2(aux,:),novr3,nfft2,fs2);
    
    contato1=mean(xDemon1);
    contato2=mean(xDemon2);
    rotacao=fDemon;
        
    MMf=[contato1;contato2];
    
    [wf,Sfreq,errsf]=m_nmf(MMf,alg,opt);
    
    Matriz_comp1_freq(aux,:)=Sfreq(1,:);%mat_comp_freq(1,:);
    Matriz_comp2_freq(aux,:)=Sfreq(2,:);%mat_comp_freq(2,:);
    Matriz_contato1_freq(aux,:)=contato1;%mat_input_freq(1,:);
    Matriz_contato2_freq(aux,:)=contato2;%mat_input_freq(2,:);
    
    fig=figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    subplot(2,2,1),plot(rotacao,10*log10(normaliza(abs(contato1),0)),'r'),xlabel('Rotation (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 1 ',aux),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(contato1),0)),rotacao, 120.1, 146.5 );
    
    subplot(2,2,3);plot(rotacao,10*log10(normaliza(abs(contato2),0)),'r'),xlabel('Rotation (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 2 ',aux),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(contato2),0)),rotacao, 120.1, 146.5 );
    
    subplot(2,2,2),plot(rotacao,10*log10(normaliza(abs(Matriz_comp1_freq(aux,:)),0)),'b'),xlabel('Rotation (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do componente 1 alg. %s (freq.)',aux,alg),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_comp1_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    
    subplot(2,2,4);plot(rotacao,10*log10(normaliza(abs(Matriz_comp2_freq(aux,:)),0)),'b'),xlabel('Rotation (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do componente 2 alg. %s (freq.)',aux,alg),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_comp2_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    inicio=fim+1;
    fim=fim+janela*fs;
    aux=aux+1;

end
