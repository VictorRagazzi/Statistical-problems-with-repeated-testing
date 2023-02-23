%% Modified Rayleigh test (St�rzebecher, E., Cebulla, M., & Elberling, C. (2005))
clear
clc
M = 240;    % N�mero de janelas
L = 1000;   % Tamanho da janela
B = 65;     % Bin da frequ�ncia de interesse
% Sinal simulado no tempo
s = 0.0*cos(2*pi*(B-1)*(1:M*L)'/L);   % Sinal a ser detectado
y = s + 1*randn(M*L,1);                 % Sinal mais ru�do
% Modified Rayleigh Test
y = reshape(y,L,M); % Dividindo o sinal em janelas
Y = fft(y);         % Aplicando a transformada discreta de Fourier
Y = Y(B,:);         % Selecionando as componenetes espectrais de cada janela na frequ�ncia de interesse.
Y = sort(Y);%,'ComparisonMethod','abs'); % Colocar em ordem crescente em rela��o ao m�dulo
Y = (Y./abs(Y)).*(1:M); % Trocando o m�dulo pelo rank
Cm = mean(real(Y));
Sm = mean(imag(Y));
Rm = (sqrt(Cm^2 + Sm^2))/(M*sqrt(M));
% Para aplicar Rayleigh test � necess�rio o valor cr�tico para o qual o
% valor de Rm ser� comparado.
%% Simula��o na frequ�ncia para a hip�tese nula
clear
clc
M = 240;    % N�mero de janelas
% Modified Rayleigh Test
Y = randn(1,M) + 1j*randn(1,M);       % Simulando as componenetes espectrais de cada janela na frequ�ncia de interesse.
Y = sort(Y);%,'ComparisonMethod','abs'); % Colocar em ordem crescente em rela��o ao m�dulo
Y = (Y./abs(Y)).*(1:M); % Trocando o m�dulo pelo rank
Cm = mean(real(Y));
Sm = mean(imag(Y));
Rm = (sqrt(Cm^2 + Sm^2))/(sqrt(M));
%% Valor cr�tico por simula��o de Monte Carlo
clear
clc
M = 100;    % N�mero de janelas
a = 0.01;   % N�vel de signific�ncia
nRuns = 10000;
Rm = zeros(nRuns,1);
for k = 1:nRuns
    % Modified Rayleigh Test
    Y = randn(1,M) + 1j*randn(1,M);       % Simulando as componenetes espectrais de cada janela na frequ�ncia de interesse.
    Y = sort(Y);%,'ComparisonMethod','abs'); % Colocar em ordem crescente em rela��o ao m�dulo
    Y = (Y./abs(Y)).*(1:M); % Trocando o m�dulo pelo rank
    Cm = mean(real(Y));
    Sm = mean(imag(Y));
    Rm(k) = (sqrt(Cm^2 + Sm^2))/sqrt(M);
end
Rm_crit = prctile(Rm,100*(1-a));%% Valor cr�tico por simula��o de Monte Carlo
%% Testando q-sample %%
clc, clear all
M = 180;    % N�mero de janelas
Mmin = 20;  % N�mero m�nimo de janelas para realizar o primeiro teste
Mstep = 2;  % Quantidade de janelas a ser acrescentada para realizar outro teste
MM = Mmin:Mstep:M;  % Quantidades de janelas que ocorrer�o os testes
QT = length(MM);    % Quantidade de testes
a = 0.01;   % N�vel de signific�ncia
nRuns = 10000;
Rm = zeros(nRuns,QT);
t = tic();
L = 1000;
bins = [80 160 240 320 400 480 560 640 720 800 880 960];
q = length(bins);
parfor k = 1:nRuns
    % Modified Rayleigh Test
    YY = randn(L,M) + 1j*randn(L,M);       % Simulando as componenetes espectrais de cada janela na frequ�ncia de interesse.
    for j = 1:QT
        Y = sort(YY(bins,1:MM(j)));%,'ComparisonMethod','abs'); % Colocar em ordem crescente em rela��o ao m�dulo
        for i = 1:size(Y,1)
            Y(i,:) = (Y(i,:)./abs(Y(i,:))).*(1:MM(j)); % Trocando o m�dulo pelo rank
        end
        Cm = sum(real(Y),2);
        Sm = sum(imag(Y),2);
        Rm(k,j) = sum(Cm.^2 + Sm.^2);
    end
end

Rm = 8*Rm/(q*q*(q+1)*(q+1)*M);
plot(mean(Rm,1))

disp(toc(t))

%% Valor cr�tico para testes sequenciais (St�rzebecher, E., Cebulla, M., & Elberling, C. (2005))
clc, clear all
M = 240;    % N�mero de janelas
Mmin = 1;  % N�mero m�nimo de janelas para realizar o primeiro teste
Mstep = 1;  % Quantidade de janelas a ser acrescentada para realizar outro teste
MM = Mmin:Mstep:M;  % Quantidades de janelas que ocorrer�o os testes
QT = length(MM);    % Quantidade de testes
a = 0.01;   % N�vel de signific�ncia
nRuns = 10000;
Rm = zeros(nRuns,QT);
t = tic();
parfor k = 1:nRuns
    % Modified Rayleigh Test
    YY = randn(1,M) + 1j*randn(1,M);       % Simulando as componenetes espectrais de cada janela na frequ�ncia de interesse.
    for j = 1:QT
        Y = sort(YY(1:MM(j)));%,'ComparisonMethod','abs'); % Colocar em ordem crescente em rela��o ao m�dulo
        Y = (Y./abs(Y)).*(1:MM(j)); % Trocando o m�dulo pelo rank
        Cm = mean(real(Y));
        Sm = mean(imag(Y));
        Rm(k,j) = (sqrt(Cm^2 + Sm^2))/(sqrt(MM(j)));
    end
end

disp(toc(t))
%% Apresenta resultados

alpha = 0.01;
figure
[f,xi] = ksdensity(Rm(:,end));
%[f,xi] = ksdensity(max(Rm0(:,:),[],2));

limiar0 = quantile(Rm(:,end), 1-alpha);
f = f/trapz(f);
plot(xi,f,'LineWidth', 2,'color','b');
hold on
plot([limiar0 limiar0],[0,0.05],'b--')

%[f2,xi2] = ksdensity(max(Rm0(:,:),[],2));
%limiar02 = quantile(f2, 1-alpha);
%plot(xi2,f2,'LineWidth', 2);
%plot([limiar02 limiar02],[0,3],'g--')

%[f3,xi3] = ksdensity(max(Rm162,[],2));
[f3,xi3] = ksdensity(max(Rm(:,1:6:91),[],2));
limiar16 = quantile(max(Rm(:,1:6:91),[],2), 1-alpha);
f3 = f3/trapz(f3);
plot(xi3,f3,'LineWidth', 2,'color','y');
plot([limiar16 limiar16],[0,0.05],'y--')

%[f4,xi4] = ksdensity(max(Rm91,[],2));
[f4,xi4] = ksdensity(max(Rm,[],2));
limiar91 = quantile(max(Rm,[],2), 1-alpha);
f4 = f4/trapz(f4);
plot(xi4,f4,'LineWidth', 2,'color','r');
plot([limiar91 limiar91],[0,0.05],'r--')


%limiarBon = quantile(Rm0(:,end),1-alpha/91);
%plot([limiarBon limiarBon],[0,0.05],'b--')

hold off
grid on
xlim([0,2])

%legenda = {'H_{01}',['Limiar_{01} =' num2str(limiar0)],...
%           'H_{0216}',['Limiar_{0216} =' num2str(limiar16)],...
%           'H_{0291}',['Limiar_{0291} =' num2str(limiar91)]};
%           'H_{01Bon}',['Limiar_{01Bon} =' num2str(limiarBon)]};

%legend(legenda,'Location','NorthWest')