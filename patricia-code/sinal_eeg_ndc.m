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
%Intensidade = {'30dB'};
% Intensidade = {'60dB'};

Mmax = 240; %valor m�ximo 

%% Parametros do protocolo de detec��o. 

%parametros = [Min Mstep Mmax alfa_corrigido]
%parametros = [200 1 200 .05];
% nRuns  = 10000;
% Mmax = 240; %n�mero m�ximo de janela 
% alfa_teste = 0.05;
% FP_desejado =0.05;
% [alfa_corrigido,NDC_minimo,cost_alfa, P] = funcao_NDC_alfaCorrigido_Mmax(nRuns,Mmax ,alfa_teste,FP_desejado);
alfa = 0.05;
FP_desejado = 0.05; 


% load(['NDC_AlfaCorrigido_MSC_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
%     'NDC_minimo','P', 'nRuns')
load(['NDC_AlfaCorrigido_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
    'NDC_minimo','P', 'nRuns')

parametros = [P, NDC_minimo,alfa_corrigido];
% parametros(:,5) = 0.05;


%% deixar o alfa corrigido em 0.05 nesse caso n�o h� corre�ao
for ii=1:size(parametros,1)
    parametros(ii,5)=0.05;
end

%%
load([caminho 'eletrodos.mat'])
%aux=[1:10];
pos_ele =1; 

ganho  = 200;
alpha = 0.05; 

%% --------------------------------------------------
remoc = [.1]/ganho; 


%% 
%******poder fazer por intensidade aqui -------
for ncanal=1:16 %eu que fiz
    for cont_vol = 1:size(Vvoluntario,1) %fazer por volunt�rio 
    
    voluntario = cell2mat(Vvoluntario(cont_vol,:)); %carregar o volunt�rio 
    intensidade = cell2mat(Intensidade); %intensidadde 
    load([caminho voluntario intensidade], 'x','Fs','binsM','freqEstim')   
  
     x = x(:,:,pos_ele); % todas as linhas e colunas e s� o eletrodo utilizado %esse � o verdadeiro
     %%%%%%x=x(:,:,ncanal); %eu que alterei
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
      
      [dr,time] = protocolo_deteccao(x, parametros);
      
      Tdr(:,:,cont_vol) = dr; %dr �=1 se atingiu ndc e zero caso contr�rio(dr � uma matriz onde
%as linhas representam as frequ�ncias e colunas representa o indice dos conjuntos de parametros)     
      Ttime(:,:,cont_vol) = time; %salva dr em fun��o do volunt�rio
    end           
     %%fazendo a an�lise de desemprenho por canal
     %da qui para baixo eu que fiz a modifica��o
%TXD(:,ncanal) = mean(mean(Tdr(binsM,:,:),3),1)';

%Tdr(bimsm,:,:) representa o indice da frequencias de estimula��o para todos os
%parametros e volunt�rios
TXD = mean(mean(Tdr(binsM,:,:),3),1)'; %calcula a media da detec��o em fun��o dos voluntarios 
%e em seguida em fun��o dos parametros

binsR = binsM+1;
binsR = 1:100; 
binsR(binsM) = []; 
binsR(1:2) = []; 

FP(:,ncanal) = mean(mean(Tdr(binsR,:,:),3),1)';

timeM = time(binsM,:); 
timeM(timeM==-1) = Mmax;
timeMedio(:,ncanal) = mean(timeM,1)'*1; %1segundo por janela

end
%% An�lise de desempenho 

%TXD - analisar as freq. estimula��o 
%binsM = [82    84    86    88    90    92    94    96]
%freq. 81Hz,83,85,87,89,91,93,95Hz
TXD = mean(mean(Tdr(binsM,:,:),3),1)';

binsR = binsM+1;
binsR = 1:100; 
binsR(binsM) = []; 
binsR(1:2) = []; 

FP = mean(mean(Tdr(binsR,:,:),3),1)';




%% mostrar resultados 
%clc
%1 - Taxa de detec��o 

% figure 
% plot(TXD_semCorrecao(:,2),'.k','MarkerSize',6)
% hold on 
% plot([0 size(TXD_corrigido,1)],[TXD_corrigido(1,2) TXD_corrigido(1,2)], ':r','LineWidth',1.5)

% plot(TXD_corrigido(:,3),'*','Color',[0 0.4470 0.7410]','Markersize',6)
% plot(TXD_corrigido(:,2),'.b','MarkerSize',6)
% plot((TXD),'*','Color',[0 0.4470 0.7410]','Markersize',6)
plot((TXD),'.b','MarkerSize',6)
xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
ylabel('Taxa de Detec��o(%)','fontsize',12)
title('INTENSIDADE 30dB SPL','Fontsize',12)
% box off
grid on

% 
% TXD_semCorrecao(:,3)=TXD;
% TXD_corrigido(:,3)=TXD;
%%

%2 - Falsos Positivo  
% figure 
% plot(FP_semCorrecao(:,2)*100,'.k','MarkerSize',8)
% hold on 
% plot(FP_corrigido(:,2)*100,'.b','Markersize',8)
% %hold on 
% plot([0 size(FP_corrigido,1)],[5 5], ':r','LineWidth',2)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',12)
% box off
% title('40 db dois resultados na mesma figura')

figure
plot(FP)
figure
% boxplot(fpboxplot)
xlabel('NDC sem corre��o')

FP_semCorrecao(:,3)=FP;
FP_corrigido(:,3)=FP;

%come�a aqui os subplots
% figure
% subplot(2,1,1) 
% plot(FP_semCorrecao(:,1)*100,'.k','MarkerSize',5)
% hold on 
% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
% plot(FP_corrigido(:,1)*100,'.b','Markersize',5)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',10)
% box off
% title('Intensidade de 50 dB SPL','Fontsize',12)

% subplot(2,1,1)
% plot(FP_semCorrecao(:,3)*100,'.k','MarkerSize',5)
% hold on 
% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
% plot(FP_corrigido(:,3)*100,'.b','Markersize',5)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',10)
% box off
% title('Intensidade de 30 dB SPL','Fontsize',12)

plot(FP_semCorrecao(:,3)*100,'.k','MarkerSize',5)
hold on 
plot([0 size(FP_semCorrecao,1)],[6.2 6.2], ':r','LineWidth',2) %Limite do SUP nivel de significancia

% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2) %Limite do nivel de significancia

plot([0 size(FP_semCorrecao,1)],[3.9 3.9], ':r','LineWidth',2) %Limite do INF nivel de significancia

plot(FP_corrigido(:,3)*100,'.b','Markersize',5)
xlabel('�ndice dos conjuntos de Par�metros','fontsize',14)
ylabel('Taxa de Falso Positivo(%)','fontsize',14)
box off
title('INTENSIDADE 30dB SPL','Fontsize',12)


% subplot(3,1,3)
% plot(FP_semCorrecao(:,3)*100,'.k','MarkerSize',7)
% hold on 
% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
% plot(FP_corrigido(:,3)*100,'.b','Markersize',7)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',10)
% box off
% title('30 db')

% %termina aqui os subplots
% figure
% subplot(3,1,1) 
% plot(FP_corrigido(:,1)*100,'.k','MarkerSize',9)
% hold on 
% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',10)
% box off

% subplot(3,1,2)
% plot(FP_corrigido(:,2)*100,'.k','MarkerSize',9)
% hold on 
% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',10)
% box off
% 
% subplot(3,1,3)
% plot(FP_corrigido(:,3)*100,'.k','MarkerSize',9)
% hold on 
% plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
% xlabel('�ndice dos conjuntos de Par�metros','fontsize',12)
% ylabel('Taxa de Falso Positivo(%)','fontsize',10)
% box off


%%
% taxa de detec��o x tempo 
timeM = time(binsM,:); 
timeM(timeM==-1) = Mmax;
timeM = mean(timeM,1)'*1; %1segundo por janela
TXD = mean(mean(Tdr(binsM,:,:),3),1)' *100;

%plot da fronteira parero
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');  %estudar
plot([0 1]*100, [Mmax Mmax],'-.k','linewidth',1) %linha do mmax 240
plot([TXD(1) TXD(1)], [min(timeM) max(timeM)],'-.b','linewidth',1) %linha do teste unico

figure
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

plot(TXD(idxs),timeM(idxs),'-or','Markersize',8,'linewidth',1.2) %plot da fronteira

% ylim([min(-p(:,2))*.9 max(-p(:,2))*1.1])
% xlim([min(p(:,1))*80 max(p(:,1))*104])

set(axes1,'XMinorTick','on');
set(axes1,'YMinorTick','on');
box(axes1,'off');
hold off
xlabel('Detection Rate (%)','fontsize',12); 
ylabel('Mean Exam Time (s)','fontsize',12);
fprintf('\n'); 

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


% % %%
% % %testando a redu��o
% % sem_correcao=TXD;
% % % clear TXD
% % % correcao=TXD;
% % 
% % %%%
% % sem_correcao=timeM;
% % % clear timeM
% % % correcao=timeM;
% % 
% % %redu��o do tempo  (maior-menor/maior) ou (tmax-tdetc/tmax)
% % redeucao=sem_correcao-correcao;
% % aumento=correcao-sem_correcao;
% % for ii=1:size(sem_correcao,1)
% %    %%reduztempo(ii)=redeucao(ii)/sem_correcao(ii); 
% %   % if sem_correcao(ii)~= 0
% %         aumentapd(ii)=aumento(ii)/sem_correcao(ii);
% %    %end
% % end   
% % aumentapd=aumentapd*100;
% % 
% % mean(reduztempo)
% % mean(aumentapd)
% % 
% % median(aumentapd)
% % 
% % x=[1:1:1342];
% % figure
% % plot(x,aumentapd,'Linewidth',1.2)
% % hold on
% % plot(x,reduztempo,'-r','Linewidth',1.2)
% % plot([1 1350],[0  0],'-k','Linewidth',1.5)
% % xlabel('Indice do Protocolo','Fontsize',12)
% % ylabel('Aumento da Pd(%)','Fontsize',12)
% % title('Ndc corrigido vs Ajuste do alfa 40db','Fontsize',12)
% % 
% % figure
% % hist(aumentapd,120)
% % title('Histogram da distribui��o do aumento da Taxa de detec��o','Fontsize',12)
% % 
% % %plotar os conjugos de parametros em contorno 3d
% % 
% % figure
% % histogram(aumentapd,100,'BinLimits',[0,80]) 
% %  
% % figure
% % histfit(aumentapd,100)
% % 
% % matriz_dados(:,4)= reduztempo;%(P(:,2))%,P(:,2),aumentapd,reduztempo);
% % 

%% plotando a distribui��o em 3d
X=matriz_dados(:,1);
Y=matriz_dados(:,2);
Z=matriz_dados(:,3);
figure
scatter3(X,Y,Z)

figure
%mapshow(X,Y,Z)
geoshow(X,Y,Z,'DisplayType','surface');

figure
usamap(X,Y)
geoshow(X,Y, Z, 'DisplayType', 'surface')
demcmap(Z)
daspectm('m',1)
view(3)

%% tentando plotar
matriz(:,:,1)=P(:,1:2);
matriz(:,:,2)=[aumentapd aumentapd];
x_grafico=P(:,1);
y_grafico=P(:,2);
z_grafico=aumentapd;

[X,Y,Z]=xyz2grid(parametros(:,1),(parametros(:,2)),aumentapd);
surf(X,Y,Z)


%% exemplo de superficie 
x=[-2:0.2:2];
y=[-2:0.2:2];
[X,Y]=meshgrid(x,y);

Z=X.*exp(-X.^2-Y.^2);
figure(1)
mesh(X,Y,Z)

figure(2)
surf(X,Y,Z)

figure(3)
mesh(Xa,Ya,ZAmem)

%% MINHA SUPERFICIE
xaux=P(:,1);
xaux=xaux';

yaux=P(:,2);
yaux=yaux';

[Xa,Ya]=meshgrid(xaux,yaux);

for linha=1:1342
    for coluna=1:1342
        Mmin=Xa(linha,coluna);
        Mstep=Ya(linha,coluna);
        
        n=(240-Mmin)/Mstep;
        inteiro=round(n);
        if n==inteiro
            [hh,gg]=find(P==Mstep);
        else
            
        end
   
    end
end