%% Tra�ando graficos 30dB
% close all;
clearvars;clc;
%% Carrega variavies (t�m que estar na pasta local do caminho)
%  01
load('Workspace_M1_30','FP','TXD', 'timeM', 'Tdr', 'Ttime');
fp1 = FP; txd1 = TXD; timeM1 = timeM; 
Tdr1 = Tdr;
Ttime1 = Ttime;
%% 02
load('Workspace_M2_30','FP','TXD', 'timeM', 'Tdr', 'Ttime');
fp2 = FP; txd2 = TXD; timeM2 = timeM;
Tdr2 = Tdr; Ttime2 = Ttime;
%% 03
load('Workspace_M3_30','FP','TXD', 'timeM', 'Tdr', 'Ttime');
fp3 = FP; txd3 = TXD; timeM3 = timeM;
Tdr3 = Tdr; Ttime3 = Ttime;
%% 04
load('Workspace_M4_30','FP','TXD', 'timeM', 'Tdr', 'Ttime');
fp4 = FP; txd4 = TXD; timeM4 = timeM;
Tdr4 = Tdr; Ttime4 = Ttime;
load('Workspace_M1_30','linf', 'binsM');
load('Workspace_M1_30','lsup');

%% Teste McNemar
% h = testcholdout(YHat1,YHat2,Y) returns the test decision, by conducting the mid-p-value McNemar test, 
% from testing the null hypothesis that the predicted class labels YHat1 and YHat2 have equal accuracy for
%  predicting the true class labels Y. The alternative hypothesis is that the labels have unequal accuracy.
% h = 1 indicates to reject the null hypothesis at the 5% significance level. 
% h = 0 indicates to not reject the null hypothesis at 5% level.
% 30dB Single test last line
%Comparar M�todo 04 x Standart
%Comparar M�todo 04 x M�todo 02
% freqEstim = [81 83 85 87 89 91 93 95]; %frequencia de estimula��o fixa 
% binsM = [freqEstim] +1;%bins de modula��o
%Tdr = 120 x 1342 x 11
% [valores, indices] = min(Tdr4); % devolve valor e indices em Tdr dos menores tempos de detec��o
% Ppareto = P(indices);

% Tdr2 = reshape( Tdr1, [88,1]);
%% Descobrir par�metros de Pareto
%Indices paretos: function: [p, idxs] = paretoFront([TXD,(-timeM)] );
%M�todo 1 [1306;1420;1501; 2142 ;2293;2516;2694] - 4� 2142 
%M�todo 2 [731;735;1014; 1602 ;2516;2694] - 4� 1602
%M�todo 3 [1606;2367;2376;2406; 2412 ;2444;2448;2466] - 5� 2412
%M�todo 4 [460;924;2183;2444;2376;2449; 2383; 2596;912] -  7� 2383 com 8� = 2596
[p, idxs] = paretoFront([TXD,(-timeM)] );
%%
p1 =  2694; p2 =  2694; p3 =  2466; p4= 912; % �ltimas posi��es 
%  p1 =  2516; p2 =  2516; p3 =  2448; p4 =  912; % Pen�ltimas posi��es
% p1 =  2142; p2 =  1602; p3 =  2412; p4 =  2383/2596; % 28,4% posi��es
%%
Vdr4 = Tdr4 ( binsM, 460 ,:);  % Selecionando freq Est. Melhor P volunt�rios
Vdr4 = Vdr4(:);             % Vetorizando
Vdr3 = Tdr3 ( binsM, p3 ,:);  % Selecionando freq Est. Melhor P volunt�rios
Vdr3 = Vdr3(:);             % Vetorizando
Vdr2 = Tdr2 ( binsM, p2 ,:);  % Selecionando freq Est. Melhor P volunt�rios
Vdr2 = Vdr2(:);             % Vetorizando
Vdr1 = Tdr1 ( binsM, p1 ,:);  % Selecionando freq Est. Melhor P volunt�rios
Vdr1 = Vdr1(:);             % Vetorizando
Vdrs = Tdr4 (binsM, end, :); % Vetor single test last line
Vdrs = Vdrs(:);             % Vetorizando
V = ones (size(Vdr1));

% h13 = testcholdout (Vdr1, Vdr3, Vdrs);
% h12 = testcholdout (Vdr1, Vdr2, Vdrs);
% h23 = testcholdout (Vdr2, Vdr3, Vdrs);
% h14 = testcholdout (Vdr1, Vdr4, Vdrs);

h = testcholdout (Vdr4, Vdr1, V);
%% Teste de Normalidade
Vtime4 = Ttime4 ( binsM, 912 ,:);  % Selecionando freq Est. Melhor P volunt�rios
Vtime4 = Vtime4(:);             % Vetorizando
%N�o � normal
%% 
% Vtime4 = Ttime4 ( binsM, 912 ,:);  % Selecionando freq Est. Melhor P volunt�rios
% Vtime4 = Vtime4(:);             % Vetorizando
Vtime3 = Ttime3 ( binsM, p3 ,:);  % Selecionando freq Est. Melhor P volunt�rios
Vtime3 = Vtime3(:);             % Vetorizando
% Vtime2 = Ttime2 ( binsM, p2 ,:);  % Selecionando freq Est. Melhor P volunt�rios
% Vtime2 = Vtime2(:);             % Vetorizando
% Vtime1 = Ttime1 ( binsM, p1 ,:);  % Selecionando freq Est. Melhor P volunt�rios
% Vtime1 = Vtime1(:);             % Vetorizando
Vtimes = 440*ones (size(Vtime4));
% V = ones (size(Vtime1));
[p,h]= ranksum(Vtime4, Vtimes);

%% Diferen�a significativa tempo


%% Tempo
t4 = min(timeM4); tm4 = mean(timeM4);
t3 = min(timeM3); tm3 = mean(timeM3);
t2 = min(timeM2); tm2 = mean(timeM2);
t1 = min(timeM1); tm1 = mean(timeM1);
mediat4 = (tm4-tm1)/tm1;
t = timeM4(end); %single test last line
menort4 = (t4-t)/t;
menort3 = (t3-t)/t;
menort12 = (t2-t)/t; %tempo min 1 e 2 s�o =
%% Taxa de Detec��o
dr4 = max(txd4);drm4 = mean(txd4);
dr3 = max(txd3);drm3 = mean(txd3);
dr2 = max(txd2);drm2 = mean(txd2);
dr1 = max(txd1);drm1 = mean(txd1);
ds = txd4(end);

dr = txd4(end); %single test last line
maiordr4 = (dr4-dr)/dr;
% maiordr3 = (dr3-dr)/dr;
% maiordr2 = (dr2-dr)/dr;
maiordr1 = (dr1-dr)/dr;
%% 

