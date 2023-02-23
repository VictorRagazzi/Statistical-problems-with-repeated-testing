%% Tra�ando graficos
close all;clearvars;clc;
% Carrega variavies (t�m que estar na pasta local do caminho)
load('Workspace_M1','FP');
fp1 = FP;
load('Workspace_M2','FP');
fp2 = FP;
load('Workspace_M3','FP');
fp3 = FP;
load('Workspace_M4','FP');
fp4 = FP;
load('Workspace_M1','linf');
load('Workspace_M1','lsup');
%%
% cria legendas
leg0 = 'FP Desejado = 5%';
legc1 = ['Confian�a = ',num2str(linf),'%'];
legc2 = ['Confian�a = ',num2str(lsup),'%'];
legC = ['Confian�a = [',num2str(linf),' - ',num2str(lsup),'] %'];
leg1 = ['M�todo 1, m�dia = ',num2str(round(mean(fp1)*100,2)),'%'];
leg2 = ['M�todo 2, m�dia = ',num2str(round(mean(fp2)*100,2)),'%'];
leg3 = ['M�todo 3, m�dia = ',num2str(round(mean(fp3)*100,2)),'%'];
leg4 = ['M�todo 4, m�dia = ',num2str(round(mean(fp4)*100,2)),'%'];
tamanho = 20;

%% plota 1
figure

%Tra�os m�dia FP e  de Limites de confian�a
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym1 = mean(fp1)*100*ones(1,length(x)); %M�dia M�todo 01
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Tra�os
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym1,'Color','k','LineStyle',':','LineWidth',2); %M�dia M�todo 01
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP
%M�todo 01
scatter(1:size(fp1,1),fp1*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);

%Legendas, t�tulos e Eixos
legend(leg0,leg1,legC)                  %M�todo 01
ylabel('Falso Positivo','fontsize',12);
xlabel('Ind�ce dos conjuntos de Par�metros','fontsize',14);
title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])
%% 2
figure

%Tra�os m�dia FP e  de Limites de confian�a
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym2 = mean(fp2)*100*ones(1,length(x)); %M�dia M�todo 02
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Tra�os
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym2,'Color','k','LineStyle',':','LineWidth',2);  %M�dia M�todo 02
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP
%M�todo 02
scatter(1:size(fp2,1),fp2*100,tamanho,'filled','CData',[0.8500, 0.3250, 0.0980]);

%Legendas, t�tulos e Eixos
legend(leg0,leg2,legC)                  %M�todo 02
ylabel('Falso Positivo','fontsize',12);
xlabel('Ind�ce dos conjuntos de Par�metros','fontsize',14);
title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])
%% 3
figure

%Tra�os m�dia FP e  de Limites de confian�a
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym3 = mean(fp3)*100*ones(1,length(x)); %M�dia M�todo 03
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Tra�os
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym3,'Color','k','LineStyle',':','LineWidth',2); %M�dia M�todo 03
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP
%M�todo 03
scatter(1:size(fp3,1),fp3*100,tamanho,'filled','CData',[0.9290, 0.6940, 0.1250]);

%Legendas, t�tulos e Eixos
legend(leg0,leg3,legC)                   %M�todo 03
ylabel('Falso Positivo','fontsize',12);
xlabel('Ind�ce dos conjuntos de Par�metros','fontsize',14);
title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])


%% 4
figure

%Tra�os m�dia FP e  de Limites de confian�a
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym3 = mean(fp4)*100*ones(1,length(x)); %M�dia M�todo 04
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Tra�os
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym3,'Color','k','LineStyle',':','LineWidth',2);%M�dia M�todo 04
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP 
%M�todo 04
scatter(1:size(fp4,1),fp4*100,tamanho,'filled','CData',[0.4940, 0.1840, 0.5560]);

%Legendas, t�tulos e Eixos
legend(leg0,leg4,legC)                  %M�todo 04
ylabel('Falso Positivo','fontsize',12);
xlabel('Ind�ce dos conjuntos de Par�metros','fontsize',14);
title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])