%clear all
clc
N = 88; %n�mero de amostras utilizadas no c�lculo da m�dia 
%N = 4032; 

p = 0.05; %n�vel de signific�ncias 
limite_conf = .9; %intervalo de confian�a
Nrun = 10000;
for ii = 1:Nrun 
r(ii) = binornd(N,p);
end
fprintf('\nLimite Interior: %d, alpha1 = %f',quantile(r,1-limite_conf)+1,(quantile(r,1-limite_conf)-1)/N)
fprintf('\nLimite de Amostras -Superior: %d, alpha Aceitavel = %f',quantile(r,limite_conf)-1,(quantile(r,limite_conf)-1)/N)



