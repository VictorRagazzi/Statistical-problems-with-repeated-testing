%compara��o do ndc unit�rio e ndc com alfa corrigido.

%%
%carregando variveis
%load('Dados_NDC_corrigido.mat')
%load('Dados_NDC_Unitario.mat')

pd_corrigido=TXD;
pd_ndc_unitario=TXD;

aumentapd=zeros(size(pd_corrigido,1),size(pd_ndc_unitario,2));
aumento=pd_corrigido-pd_ndc_unitario;
for ff=1:3
    for ii=1:size(pd_corrigido,1)
       
        if pd_ndc_unitario(ii,ff)~=0 %alguns protocolos em 30db podem ter txd igual a zero
            aumentapd(ii,ff)=aumento(ii,ff)/pd_ndc_unitario(ii,ff);
       end
    end
end
aumentapd=aumentapd*100;
% por coluna 50db 40db e 30db
%Na m�dia o ndc unit�rio teve uma pd ligeiramente maior que o ndc com
%corre��o
mean(aumentapd)
 
%mediana igual a zero em todas as intensidades garante a equivalencia
median(aumentapd)

%a partir dos valores de minimo e m�ximo pode perceber que tem grandes
%extremos ou seja protocolos em que cada t�cnica pode ser considerada coomo
%otima 
min(aumentapd) 
max(aumentapd) %para intensidades menores ndc corrigido permite um ganho maior de pd 

% novamente 50 40 e 30db 
%distribui��o do aumento da taxa de detec��o
figure %a figura mostra que em grande parte dos protocolos teve equival�ncia ou pequenas varia��es na txd
hist(aumentapd(:,1),100)
title('Distribui��o do aumento da txd 50db')

figure
hist(aumentapd(:,2),100)
title('Distribui��o do aumento da txd 40db')

figure
hist(aumentapd(:,3),100)
title('Distribui��o do aumento da pd 30db')

%% analisando o tempo

%tempo_ndc_unitario=timeMedio;
%tempo_NDC_alfaCorrigdo=timeMedio;

reduztempo=zeros(size(pd_corrigido,1),size(pd_corrigido,2));
%reducao=tempo_NDC_semCorrecao-tempo_NDC_alfaCorrigdo;
for gg=1:3
  reducao(:,gg)=tempo_ndc_unitario(:,gg)-tempo_NDC_alfaCorrigdo(:,gg);
    for hh=1:size(pd_corrigido,1)
       reduztempo(hh,gg)=reducao(hh,gg)/tempo_NDC_alfaCorrigdo(hh,gg);
  end
end

%por ordem 50db 40db e 30db
reduztempo=reduztempo*100;

mean(reduztempo)%para intensidade de 50 db a redu��o no tempo � maior e chega a 12,93% na compara��o do nc unit�rio
%em rela��o ao ndc corrigido, pode associar isso ao fato(hip�tese) da
%altera��o do crit�rio de parada ou seja no ndc unit�rio apenas um teste
%positivo � suficiente para a rejei��o da hip�tese nula, mesmo que o teste
%sej a mais rigroso uma vez que diminui o nivel de signific�ncia de cada
%teste em melhores rela��es sinal ru�do garante a redu��o no tempo 
 
median(reduztempo) 

min(reduztempo)
max(reduztempo) 

%distribui��o da redu��o do tempo
figure % teve mais protocolos com redu��o
hist(reduztempo(:,1),100)
xlabel('Redu��o do tempo (%)')
ylabel('Numero de protocolos')
title('Distribui��o da redu��o do tempo 50db')

figure
hist(reduztempo(:,2),100) %ndc unit�rio mais rapido em todos os protocolos
title('Distribui��o da redu��o do tempo 40db')
xlabel('Redu��o do tempo (%)')
ylabel('Numero de protocolos')

figure
hist(reduztempo(:,3),100)
title('Distribui��o da redu��o do tempo 30db')
xlabel('Redu��o do tempo (%)')
ylabel('Numero de protocolos')
%% Plotando os protocolos que obtem o melhor desempenho txd

X=parametros(:,1:2);
X(:,3)=aumentapd(:,1);
X(:,4)=aumentapd(:,2);
X(:,5)=aumentapd(:,3);


%observa que para protocolos em que mstep � pequeno ou seja em que s�o
%feitos poucos testes o ndc unitario teve melhor desempenho e para
%protocolos com mstep pequeno ndc com corre��o teve melhor desempenho, o
%numero minimo de janelas n�o se mostrou como fator que possa provocar
%diferen�a na probabilidade de detec��o(isso comparando as estrat�gias).


for jj=3:5
    for ii = 1:size(X,1)
    
        if X(ii,jj)<-10  %menor que 10% de aumento
             cor(ii,:) = [0.4660 0.6740 0.1880]; %cor verde representa onde a pd do ndc unitario foi maior
            tamanho(ii,1) = 70;
        elseif X(ii,jj)<10 %entre -10 e 10 fica considerado como parametro estavel cor azul
            cor(ii,:) = [0 0.4470 0.7410];
            tamanho(ii,1) = 10;
        %elseif  X(ii,3)<10
     %   cor(ii,:) = [0.9290 0.6940 0.1250];
     %   tamanho(ii,1) = 55;
        else
            cor(ii,:) = [0.6350 0.0780 0.1840]; %cor vermelha ndc alfa corrigido teve melhor desempenho
            tamanho(ii,1) = 75;
        end
        
    end
 %por ordem de figura 50 40 e 30db   
figure,
s=tamanho;
scatter3(X(:,1),X(:,2),X(:,jj),s,cor)
xlabel('M_{MIN}')
ylabel('M_{STEP}')
zlabel('Aumento PD(%)')    
text(180,200,'NDC unitario teve aumento da PD acima dos 10%','Color','green','FontSize',10)
text(180,210,'limiar "est�vel"','Color','blue','FontSize',10)
text(180,220,'NDC alfa corrigido teve aumento na pd acima dos 10%','Color','red','FontSize',10)
view(0,90)    
title('An�lise protocolo de detec��o')
end
%a vista superior fornece a imagem semelhenate ao exemplo proposto
figure,
s=tamanho;
scatter3(X(:,1),X(:,2),X(:,3),s,cor)
xlabel('M_{MIN}')
ylabel('M_{STEP}')
zlabel('Aumento PD(%)')
%legend({'Line 1','Line 2','Line 3'})
%title('Protocolos de detec��o 50db')
view(0,90)
colorbar


%percebe que para 30db tem mais protocolos com valores acima do limiar (ou seja fora da faixa de -15% a 15%) pode ser percebido tamb�m no histogram em que este est� mais distribuido
%nesse caso mesmos protocolos com mstep pequeno o ndc unit�rio foi melhor e
%os protocolos com poucos testes mantiverem mehlhor o ndc unitario
%j� o ndc com corre��o manteve melhor o desempenho apenas em protocolos com
%v�rios testes
%% %% Plotando os protocolos que obtem o melhor desempenho em rela��o ao tempo

Y=parametros(:,1:2);
Y(:,3)=reduztempo(:,1);
Y(:,4)=reduztempo(:,2);
Y(:,5)=reduztempo(:,3);

for jj=3:5
    for ii = 1:size(Y,1)
    
        if Y(ii,jj)<-10  %menor que 10% de aumento cor verde
             cor(ii,:) = [0.4660 0.6740 0.1880];
            tamanho(ii,1) = 70;
        elseif Y(ii,jj)<10 %entre -10 e 10 fica considerado como parametro est�vel
            cor(ii,:) = [0 0.4470 0.7410];
            tamanho(ii,1) = 10; %cor azul
        %elseif  X(ii,3)<10
     %   cor(ii,:) = [0.9290 0.6940 0.1250];
     %   tamanho(ii,1) = 55;
        else
            cor(ii,:) = [0.6350 0.0780 0.1840]; %cor vermelha
            tamanho(ii,1) = 75;
        end
        
    end
    
figure,
s=tamanho;
scatter3(Y(:,1),Y(:,2),Y(:,jj),s,cor)
xlabel('M_{MIN}')
ylabel('M_{STEP}')
zlabel('Redu��o do tempo(%)')    
view(0,90)    
end
%em 50db varios protocolos do ndc unitario tiverem redu��o do tempo
%em 30db paticamente todos os protocolos mantiverem dentro do limiar
%em nenhum caso o ndc com corre��o foi mais r�pido

%% analise do tempo e txd. Fronteira pareto

TXD=auxpd50db(:,1);
timeM=auxtime50db(:,1);


figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');  %estudar
plot([0 1]*100, [Mmax Mmax],'-.k','linewidth',1) %linha do mmax 240
plot([TXD(1) TXD(1)], [min(timeM) max(timeM)],'-.k','linewidth',1) %linha do teste unico
%Teve altera��o nessa parte

for hh=1:3
    TXD=auxpd50db(:,hh);
    timeM=auxtime50db(:,hh);
    
    %plota os segmentos em pontos
     if  hh==1 
        for ii = 1:size(parametros,1)
            plot(TXD(ii),timeM(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
        end
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
    if hh==1
        plot(TXD(idxs),timeM(idxs),'-or','Markersize',8,'linewidth',1.2) %plot da fronteira
        %legend('NDC Com Alfa Corrigido');
    elseif hh==2
        hold on
        plot(TXD(idxs),timeM(idxs),'-xb','Markersize',8,'linewidth',1.2) %plot da fronteira
        %legend('NDC Sem Corre��o');
    elseif hh==3
        hold on
        plot(TXD(idxs),timeM(idxs),'-xk','Markersize',8,'linewidth',1.2) %plot da fronteira
        %legend('NDC Sem Corre��o');
    else
        hold on
        plot(TXD(idxs),timeM(idxs),'-^g','Markersize',8,'linewidth',1.2) %plot da fronteira
        %ndc unitario
    end


    for ii = 1:size(idxs,1)
     [I] = find((TXD == TXD(idxs(ii))) & (timeM==timeM(idxs(ii))));
     fprintf('\nPD = %f Tempo = %f ',TXD(idxs(ii)),timeM(idxs(ii))); 
     fprintf(' NI = %d ',length(I)); 
     I = I(1); 
        for jj = I  
        
            fprintf(' - Buffer:%d, M_step:%d', parametros(jj,1),parametros(jj,2)); 
            if hh==1
                text(TXD(jj),timeM(idxs(ii))*.995,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ],'Fontsize',8); %plotar o texto dos procolos
            elseif hh==2
                 text(TXD(jj)*.994,timeM(idxs(ii))*1.005,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ],'Fontsize',8); %plotar o texto dos procolos
            else
                text(TXD(jj)*.994,timeM(idxs(ii))*.995,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ],'Fontsize',8); %plotar o texto dos procolos
            end
        end    

    end
end
hold off
xlabel('Taxa de Detec��o(%)','fontsize',12); 
ylabel('Tempo de Exame M�dio (min)','fontsize',12);
%legend({'NDC Com Alfa Corrigido','NDC Sem Corre��o'},'Location','southeast')
%legend('boxoff')
text(max(auxpd50db(:,1)),246,'NDC com ajuste do n�vel de signific�ncia','Color','red','FontSize',10)
text(max(auxpd50db(:,1)),249,'NDC sem corre��o','Color','blue','FontSize',10)
text(max(auxpd50db(:,1)),243,'NDC Unit�rio','Color','green','FontSize',10)%trocar pra auxpd 50db
grid on
fprintf('\n');  
xlim([min(TXD(idxs))*.95,max(TXD(idxs))*1.08]) %alterar par 0.95
ylim([min(timeM(idxs))*.95 Mmax*1.05])
title('An�lise fronteira pareto 50db')



%% salvando variaveis
%auxtimec4(:,3)=timeMedio(:,7);
%auxpdc4(:,3)=TXD(:,7);
