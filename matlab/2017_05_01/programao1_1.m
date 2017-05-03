%Junta todos os programas para realizar todas as análises dos sinais
%simulados

clear all
close all

addpath('functions/comuns');
addpath('functions/demon');

load Matriz_190_205.mat %Carrega as marcações já em segmentos de 60s

scrsz = get(0,'ScreenSize');
set(0,'DefaultAxesFontName', 'Helvetica')
set(0,'DefaultAxesFontSize',14);
fs=31250; %Frequemcia de amostragem do sinal
tamanho_janela=size(Matriz1,2)/fs; %Tamanho da janela em segundos
taxa_r=625;
indice_janela=size(Matriz1,1);
inicio=1;
fim=120;
bins_freq=fim-inicio+1;

Matriz_comp1=zeros(size(Matriz1,1),tamanho_janela*fs/taxa_r); %Matriz de armazenamento do componente 1 no domínio do tempo
Matriz_comp2=zeros(size(Matriz2,1),tamanho_janela*fs/taxa_r); %Matriz de armazenamento do componente 2 no domínio do tempo
Matriz_comp3=zeros(size(Matriz2,1),tamanho_janela*fs/taxa_r); %Matriz de armazenamento do componente 3 no domínio do tempo
Matriz_ream1=zeros(size(Matriz1,1),tamanho_janela*fs/taxa_r); %Matriz de armazenamento do contato 1 no domínio do tempo
Matriz_ream2=zeros(size(Matriz1,1),tamanho_janela*fs/taxa_r); %Matriz de armazenamento do contato 2 no domínio do tempo
Matriz_ream3=zeros(size(Matriz1,1),tamanho_janela*fs/taxa_r); %Matriz de armazenamento do contato 3 no domínio do tempo
Matriz_comp1_freq=zeros(size(Matriz1,1),bins_freq); %Matriz de armazenamento do componente 1 no domínio da frequencia
Matriz_comp2_freq=zeros(size(Matriz2,1),bins_freq); %Matriz de armazenamento do componente 2 no domínio da frequencia
Matriz_comp3_freq=zeros(size(Matriz2,1),bins_freq); %Matriz de armazenamento do componente 3 no domínio da frequencia
Matriz_contato1_freq=zeros(size(Matriz1,1),bins_freq); %Matriz de armazenamento do contato 1 no domínio da frequencia
Matriz_contato2_freq=zeros(size(Matriz1,1),bins_freq); %Matriz de armazenamento do contato 2 no domínio da frequencia
Matriz_contato3_freq=zeros(size(Matriz1,1),bins_freq); %Matriz de armazenamento do contato 3 no domínio da frequencia
vet_kl_ica=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre os componentes no dominio do tempo
vet_kl_mist=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre as observações no dominio do tempo
vet_kl_ica_freq=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre os componentes no dominio da frequencia
vet_kl_mist_freq=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre as observações no dominio da frequencia
vet_kl_comp_total=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre os componentes das janelas totais no domínio do tempo 
vet_kl_padrao_ouro1=zeros(1,indice_janela); %Vetoer de armazenamento da div. KL entre o padrao ouro e o componente 1 no dominio do tempo
vet_kl_padrao_ouro2=zeros(1,indice_janela); %Vetoer de armazenamento da div. KL entre o padrao ouro e o componente 2 no dominio do tempo
vet_kl_padrao_ouro1_freq=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre o padrao ouro e o componente 1 no dominio da frequencia
vet_kl_padrao_ouro2_freq=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre o padrao ouro e o componente 2 no dominio da frequencia
vet_ems_comp1=zeros(1,indice_janela); %Vetor de aramzenamento do erro medio quadratico entre o padrão ouro e o componente 1 no domminio do tempo
vet_ems_comp2=zeros(1,indice_janela); %Vetor de aramzenamento do erro medio quadratico entre o padrão ouro e o componente 2 no domminio do tempo
vet_ems_comp1_freq=zeros(1,indice_janela); %Vetor de aramzenamento do erro medio quadratico entre o padrão ouro e o componente 1 no domminio da frequencia
vet_ems_comp2_freq=zeros(1,indice_janela); %Vetor de aramzenamento do erro medio quadratico entre o padrão ouro e o componente 2 no domminio da frequencia
vet_ems_mist1_freq=zeros(1,indice_janela); %Vetor de aramzenamento do erro medio quadratico entre o padrão ouro e a observacao 1 no domminio da frequencia
vet_ems_mist2_freq=zeros(1,indice_janela); %Vetor de aramzenamento do erro medio quadratico entre o padrão ouro e a observacao 2 no domminio da frequencia
vet_kl_padrao_ouro1_obs1=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre o padrao ouro e 0 contato 1 no dominio do tempo
vet_kl_padrao_ouro2_obs2=zeros(1,indice_janela); %Vetor de armazenamento da div. KL entre o padrao ouro e 0 contato 2 no dominio da tempo

r=3;


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

    
for aux1=1:size(Matriz1,1)
    
    [y0,y1,y2,Fmax,fs2,nfft2,novr3]=reamostra1(Matriz1(aux1,:),Matriz2(aux1,:),Matriz3(aux1,:));
    MMt=[y0;y1;y2];
    [wt,Stempo,errst]=m_nmf(MMt,alg,opt);
    Matriz_ream1(aux1,:)=y0;%mat_input(1,:);
    Matriz_ream2(aux1,:)=y1;%mat_input(2,:);
    Matriz_ream3(aux1,:)=y2;
    mat_comp(1,:)=Stempo(1,:);
    mat_comp(2,:)=Stempo(2,:);
    mat_comp(3,:)=Stempo(3,:);
    Matriz_comp1(aux1,:)=Stempo(1,:);%mat_comp(1,:);
    Matriz_comp2(aux1,:)=Stempo(2,:);%mat_comp(2,:);
    Matriz_comp3(aux1,:)=Stempo(3,:);%mat_comp(3,:);
    [xDemon1,fDemon]=calcula_espectro(Matriz_ream1(aux1,:),novr3,nfft2,fs2);
    [xDemon2]=calcula_espectro(Matriz_ream2(aux1,:),novr3,nfft2,fs2);
    [xDemon3]=calcula_espectro(Matriz_ream3(aux1,:),novr3,nfft2,fs2);

    contato1=mean(xDemon1);
    contato2=mean(xDemon2);
    contato3=mean(xDemon3);
    contato1=contato1(inicio:fim);
    contato2=contato2(inicio:fim);
    contato3=contato3(inicio:fim);
    rotacao=fDemon(inicio:fim);
    
    [xDemon_comp1]=calcula_espectro(Matriz_comp1(aux1,:),novr3,nfft2,fs2);
    [xDemon_comp2]=calcula_espectro(Matriz_comp2(aux1,:),novr3,nfft2,fs2);
    [xDemon_comp3]=calcula_espectro(Matriz_comp3(aux1,:),novr3,nfft2,fs2);
    contato_comp1=mean(xDemon_comp1);
    contato_comp2=mean(xDemon_comp2);
    contato_comp3=mean(xDemon_comp3);
    contato_comp1=contato_comp1(inicio:fim);
    contato_comp2=contato_comp2(inicio:fim);
    contato_comp3=contato_comp3(inicio:fim);
    %vet_kl_ica(aux1)=calcula_kl_freq(contato_comp1,contato_comp2);
    
    MMf=[contato1;contato2;contato3];
    
    [wf,Sfreq,errsf]=m_nmf(MMf,alg,opt);
    %[W,Sfreq,errs,varout] = nmf_alg(MMf,r,'niter',1000,'verb',3,'thresh',10e-8,'H0',H0, 'W0',W0,'norm_w',1); %Usando a div KL
    %[W,Sfreq,errs,vout]=nmf_amari(MMf,r,'alpha',2,'niter',1000,'verb',3,'thresh',1e-8,'H0',H0, 'W0',W0,'norm_w',1);%, 'norm_w',0,'norm_h',0);     %Usando a div alpha

    %[mat_input_freq, mat_comp_freq,coef_freq]=corr_ord(contato1,contato2,Sfreq(1,:),Sfreq(2,:));
    Matriz_comp1_freq(aux1,:)=Sfreq(1,:);%mat_comp_freq(1,:);
    Matriz_comp2_freq(aux1,:)=Sfreq(2,:);%mat_comp_freq(2,:);
    Matriz_comp3_freq(aux1,:)=Sfreq(3,:);%mat_comp_freq(2,:);
    Matriz_contato1_freq(aux1,:)=contato1;%mat_input_freq(1,:);
    Matriz_contato2_freq(aux1,:)=contato2;%mat_input_freq(2,:);
    Matriz_contato3_freq(aux1,:)=contato3;%mat_input_freq(2,:);
    
    
    fig=figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    subplot(3,2,1),plot(rotacao,10*log10(normaliza(abs(contato1),0)),'r'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 1 ',aux1),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(contato1),0)),rotacao, 120.1, 146.5 );
    
    subplot(3,2,3);plot(rotacao,10*log10(normaliza(abs(contato2),0)),'r'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 2 ',aux1),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(contato2),0)),rotacao, 120.1, 146.5 );
    
    subplot(3,2,5);plot(rotacao,10*log10(normaliza(abs(contato3),0)),'r'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 3 ',aux1),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(contato2),0)),rotacao, 120.1, 146.5 );
    
    for i = 1:size(Stempo,1)
        sub_ind = 2*i;
        aux = calcula_espectro(mat_comp(i,:),novr3,nfft2,fs2);
        componente=mean(aux);
        componente=componente(inicio:fim);
        subplot(size(Stempo,1),2,sub_ind),plot(rotacao,10*log10(normaliza(abs(componente),0))),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do componente %i alg. %s (tempo)',aux1,i,alg),'FONTSIZE',14),grid
        %m_func_write_text1(10*log10(normaliza(abs(componente),0)),rotacao, 120.1, 146.5);
        
    end
    xlabel('Rotação (rpm)','FONTSIZE',14);
    %saveas(fig,sprintf('./figuras_final/espectros_exp1_fig/espectro_%st_janela%i.fig',alg,aux1),'fig');
    %saveas(fig,sprintf('./figuras_final/espectros_exp_eps/espectro_%st_janela%i.eps',alg,aux1),'eps');
    %saveas(fig,sprintf('./figuras_final/espectros_exp_jpg/espectro_%st_janela%i.jpg',alg,aux1),'jpg');
end
%close all

for aux1=1:size(Matriz1,1)
    fig=figure('Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)]);
    subplot(3,2,1),plot(rotacao,10*log10(normaliza(abs(Matriz_contato1_freq(aux1,:)),0)),'r'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 1 ',aux1),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_contato1_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    subplot(3,2,3);plot(rotacao,10*log10(normaliza(abs(Matriz_contato2_freq(aux1,:)),0)),'r'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 2 ',aux1),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_contato2_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    subplot(3,2,5);plot(rotacao,10*log10(normaliza(abs(Matriz_contato3_freq(aux1,:)),0)),'r'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do contato 3 ',aux1),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_contato2_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    subplot(3,2,2),plot(rotacao,10*log10(normaliza(abs(Matriz_comp1_freq(aux1,:)),0)),'b'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do componente 1 alg. %s (freq.)',aux1,alg),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_comp1_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    
    subplot(3,2,4);plot(rotacao,10*log10(normaliza(abs(Matriz_comp2_freq(aux1,:)),0)),'b'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do componente 2 alg. %s (freq.)',aux1,alg),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_comp2_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    subplot(3,2,6);plot(rotacao,10*log10(normaliza(abs(Matriz_comp3_freq(aux1,:)),0)),'b'),xlabel('Rotação (rpm)','FONTSIZE',14),ylabel('Amplitude (dB)','FONTSIZE',14),title(sprintf('Janela %i do componente 3 alg. %s (freq.)',aux1,alg),'FONTSIZE',14),grid
    %m_func_write_text1(10*log10(normaliza(abs(Matriz_comp2_freq(aux1,:)),0)),rotacao, 120.1, 146.5 );
    
    xlabel('Rotação (rpm)','FONTSIZE',14);
    %saveas(fig,sprintf('./figuras_final/espectros_exp1_fig/espectro_%sf_janela%i.fig',alg,aux1),'fig');
    %saveas(fig,sprintf('./figuras_final/espectros_exp_eps/espectro_%sf_janela%i.eps',alg,aux1),'eps');
    %saveas(fig,sprintf('./figuras_final/espectros_exp_jpg/espectro_%sf_janela%i.jpg',alg,aux1),'jpg');
end


