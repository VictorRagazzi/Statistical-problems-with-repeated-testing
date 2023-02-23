%Aplica��o do protocolo ao banco de dados 
clear all, close all, clc 

caminho = 'C:\Users\Crohma\Documents\PNV\Simulation\Testes_Seq\Sinais_EEG\';

%vetor dos volunt�rios 
Vvoluntario = {'Ab';'An';'Bb';'Er';'Lu';...
    'So';'Qu';'Vi';'Sa';'Ti';'Wr'}; %vetor dos volunt�rio 

% Vvoluntario = Vvoluntario([1,2]);

%Intensidade ----------------------
%intensidade = {'70dB';'60dB';'50dB';'40dB';'30dB';'ESP'}; %quais intensidade analisadas 
%sugest�o
%vetor_Mmax = [50;50;240;440;440;20]; %n�mero m�ximo de janela para cada intensidade
Intensidade = {'50dB'};
% Intensidade = {'ESP'};
% Intensidade = {'30dB'};

Mmax = 240; %valor m�ximo
alfa = 0.05;
FP_desejado = 0.05; 
%% Parametros do protocolo de n�o detec��o. 
load(['VC_not_Unitario_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa) '.mat'],'VC_not')
% load(['VC_Unitario_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '.mat'],'VC_not')
%% Parametros do protocolo de detec��o. 

%% Dados 440 janelas
% 
% load(['NDC2_Unitario_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
%      'P', 'nRuns')
% NDC_minimo = ones (size(alfa_corrigido));
% parametros = [P, NDC_minimo,alfa_corrigido];

%% Dados para 240 janelas 
load(['NDC_Unitario_MSC_Mmax' num2str(Mmax) '.mat'],'parametros')

%%
load([caminho 'eletrodos.mat'])
pos_ele = 1; 

ganho  = 200;
alpha = 0.05; 

%% --------------------------------------------------
remoc = [.1]/ganho; 


%% 
%******poder fazer por intensidade aqui -------

for cont_vol = 1:size(Vvoluntario,1) %fazer por volunt�rio 
    
    voluntario = cell2mat(Vvoluntario(cont_vol,:)); %carregar o volunt�rio 
    intensidade = cell2mat(Intensidade); %intensidadde 
    load([caminho voluntario intensidade], 'x','Fs','binsM','freqEstim')   
  
    x = x(:,:,pos_ele);
    
     nfft = Fs;%1segundo de sinal 
         
     %retirar componente DC por janela (fiz isso pq no processamento em
     %tempo real � por janela)
     x = x - repmat(mean(x),nfft,1); %tirar a m�dia de cada de cada trecho - devido a remo��o
        
     %scluir os dois primeiros segundos do inicio da coleta 
     x(:,1:2,:) =[]; 
        
        
     %encontrar o valor m�ximo por canal 
      Vmax = max(abs(x),[],1);
      ind = Vmax>remoc;
      [sum(ind) cont_vol ];
      x = x(:,~ind); %removor o ru�do amplitude 
      x = x(:,1:Mmax);%limitar o tamanho para o valor m�ximo. 
     
      %******** fazer por canal diferente ----for nCanal = 1:16 %
      
      %%% INSERIR VC_not %%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
%       [dr,time] = protocolo_deteccao(x, parametros);
      [dr,time] = protocolo_deteccao(x, parametros, VC_not);
      
      Tdr(:,:,cont_vol) = dr;
      Ttime(:,:,cont_vol) = time;
            
end

%% An�lise de desempenho 

%TXD - analisar as freq. estimula��o 
%binsM = [82    84    86    88    90    92    94    96]
%freq. 81Hz,83,85,87,89,91,93,95Hz
TXD = mean(mean(Tdr(binsM,:,:),3),1)';

% binsR = binsM+1;
% binsR = 1:100; 
% binsR(binsM) = []; 
% binsR(1:2) = []; 
binsR = 70:104;
binsR(binsM - 71) = [];

FP = mean(mean(Tdr(binsR,:,:),3),1)';

%% mostrar resultados 
%clc
%1 - Taxa de detec��o 
figure 
plot(TXD*100,'.k','MarkerSize',10)
hold on 
%% Muda com Banco de DADOS
%Dados 30dB 440 Janelas descomentar : Mmax=Min na last line
plot([0 size(TXD,1)],[TXD(end)*100 TXD(end)*100], ':r','LineWidth',2)  %single shot 30dB
% Dados 50dB 240 Janelas : Mmax=Min na first line
% plot([0 size(TXD,1)],[TXD(1)*100 TXD(1)*100], ':r','LineWidth',2)  %single shot 50dB



% ylabel('Taxa de Detec��o','fontsize',12)
xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
ylabel('Taxa de Detec��o (%)','fontsize',12)
title(['Est�mulo ', (Intensidade{1}), ' SPL: M�todo 04'])
% box off
grid on
%% 2 - Falsos Positivo  
figure 
plot(FP*100,'.k','MarkerSize',10)
hold on 

lsup = 6.73;
linf = 3.37;

plot([0 size(FP,1)],[lsup lsup], ':r','LineWidth',2) %Limite do SUP nivel de significancia
% plot([0 size(FP,1)],[FP_desejado FP_desejado], ':r','LineWidth',2) %Limite do nivel de significancia
plot([0 size(FP,1)],[linf linf], ':r','LineWidth',2) %Limite do INF nivel de significancia

%ylabel('Falso Positivo','fontsize',12)
xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
ylabel('Taxa de Falso Positivo (%)','fontsize',12)
box off
title(['Est�mulo ', (Intensidade{1}), ' SPL: M�todo 04'])
ylim([0  9])
%figure
% boxplot(FP)
grid on


%%
% taxa de detec��o x tempo 
timeM = time(binsM,:); 
timeM(timeM==-1) = Mmax;
timeM = mean(timeM,1)'*1; %1segundo por janela
TXD = mean(mean(Tdr(binsM,:,:),3),1)' *100;

%% Single Shot - Muda com Banco de DADOS

figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');  
plot([0 1]*100, [Mmax Mmax],'-.b','linewidth',1)
%Dados 30dB 440 Janelas descomentar
% plot([TXD(end) TXD(end)], [min(timeM) max(timeM)],'-.b','linewidth',1) %30dB Last line
% Dados 50dB 240 Janelas
plot([TXD(1) TXD(1)], [min(timeM) max(timeM)],'-.b','linewidth',1) %50dB 1�Line


for ii = 1:size(parametros,1)
    plot(TXD(ii),timeM(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
end

[ p, idxs] = paretoFront([TXD,(-timeM)] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];

[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];

plot(TXD(idxs),timeM(idxs),'-or','Markersize',8,'linewidth',1.2) 

% ylim([min(-p(:,2))*.9 max(-p(:,2))*1.1])
% xlim([min(p(:,1))*80 max(p(:,1))*104])

set(axes1,'XMinorTick','on');
set(axes1,'YMinorTick','on');
box(axes1,'off');
hold off
xlabel('Taxa de Detec��o (%)','fontsize',12); 
ylabel('Tempo M�dio de Exame (s)','fontsize',12);
% title('INTENSIDADE 30dB SPL: M�todo 04','Fontsize',10)
title(['Est�mulo ', (Intensidade{1}), ' SPL: M�todo 04'])

for ii = 1:size(idxs,1)

    [I] = find((TXD == TXD(idxs(ii))) & (timeM==timeM(idxs(ii))));
    
    fprintf('\nPD = %f Tempo = %f ',TXD(idxs(ii)),timeM(idxs(ii))); 
    fprintf(' NI = %d ',length(I)); 
     I = I(1); 
    
    for jj = I  
        
        fprintf(' - Buffer:%d, M_step:%d', parametros(jj,1),parametros(jj,2)); 
    %    text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ','  num2str(parametros(jj,3)) '\}' ]);
%            text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ',100,'  num2str(parametros(jj,3)) '\}' ]);
%        text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ',100\}' ]);
       text(TXD(jj),timeM(idxs(ii))*.975,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ]);
    end    
end
fprintf('\n'); 
xlim([min(TXD(idxs))*.95,max(TXD(idxs))*1.05])
ylim([min(timeM(idxs))*.95 Mmax*1.05])
%%
% save(['timeM_30db_M4' '.mat'],'timeM');
% save(['pareto_30db_M4' '.mat'],'TXD','Mmax','parametros');

% save(['timeM_50db_M4' '.mat'],'timeM');
% save(['pareto_50db_M4' '.mat'],'TXD','Mmax','parametros');
%% 
% save(['Workspace_NND5_M4_50' '.mat']); %50dB
% save(['Workspace_M4_30' '.mat']); %30dB

