%função da estratégia para testes sequenciais 
function [dr,time] = ETS_2013(ord,MM,alfa,NDC,L)

%parametros: 
% MM = Mmin:Mstep:Mmax; %parametros aonde vai aplicar o teste


%%

NDC = ceil(NDC);
if size(MM,2)>1
    MM = MM';
end

if size(ord,2)>12
    ord = ord';
end


%VALOR CRÍTICO - MUDA PARA CADA DETECTOR--------
%valor_critico = VC_MSC(MM,alfa);  
for pp = 1:length(MM)
    valor_critico = L(pp);
    if pp == 1
       det = ord(MM(pp))>=valor_critico;
    end
    if pp > 1
       det = [det ord(MM(pp))>= valor_critico]; 
    end
 end
%-----------------------------------------------------
%avaliar se atende o NDC 
cont = 0; 
dr =0; 
for ii = 1:size(MM,1) 
   
    cont = det(ii)+cont*det(ii); %conta o número de detecções consecutivas 
    
    if cont == NDC
        dr = 1;
        time = MM(ii);
        break
    end
end

if dr == 0 
    time = MM(ii);
end

 