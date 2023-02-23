%fun��o da estrat�gia para testes sequenciais 
function [dr,time] = ETS_Pnv(ord,MM,alfa,NDC, VC_not)

%parametros: 
% MM = Mmin:Mstep:Mmax; %parametros aonde vai aplicar o teste
load(['VC_not_Unitario_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa) '.mat'],'VC_not') 

%% 

NDC = ceil(NDC);
if size(MM,2)>1
    MM = MM';
end

if size(ord,2)>1
    ord = ord';
end


%VALOR CR�TICO - MUDA PARA CADA DETECTOR--------
valor_critico = VC_MSC(MM,alfa);  
det = ord(MM)> valor_critico;
ndet = ord(MM) < VC_not(MM);
%-----------------------------------------------------

%avaliar se atende o NDC 
cont = 0; 
dr =0; 
for ii = 1:size(MM,1) 
   
    cont = det(ii)+cont*det(ii); %conta o n�mero de detec��es consecutivas 
%     cont1 = ndet(ii); %+cont*ndet(ii);  %conta o n�mero de n�o detec��es
    if cont == NDC     %NDC =1
        dr = 1;
        time = MM(ii);
        break
    end    
    if ndet(ii) == 1; % conferir sqdo sai 
        dr = 0;
        time = MM (ii);
        break
    end
end

if ii == size (MM,1) 
    time = MM(ii);
end

 