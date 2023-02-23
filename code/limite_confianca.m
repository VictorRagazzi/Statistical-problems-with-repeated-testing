%clear all
clc
N = 88; %número de amostras utilizadas no cálculo da média 
%N = 4032; 

p = 0.05; %nível de significãncias 
limite_conf = .9; %intervalo de confiança
Nrun = 10000;
for ii = 1:Nrun 
r(ii) = binornd(N,p);
end
fprintf('\nLimite Interior: %d, alpha1 = %f',quantile(r,1-limite_conf)+1,(quantile(r,1-limite_conf)-1)/N)
fprintf('\nLimite de Amostras -Superior: %d, alpha Aceitavel = %f',quantile(r,limite_conf)-1,(quantile(r,limite_conf)-1)/N)



