%%
clearvars
linha = [line(),line(),line(),line()];
close all
clc

%%
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');

cor = ['m','r','g','b'];

for metodo = 1:3
    caminho1 = ['timeM_30db_M',num2str(metodo),'.mat'];
    caminho2 = ['pareto_30db_M',num2str(metodo), '.mat'];
    
    load(caminho1,'timeM')
    load(caminho2,'TXD','Mmax','parametros');
    
    % pontos e linhas pretos:
    plot([0 1]*100, [Mmax Mmax],'-.k','linewidth',1)
    hold on
    plot([TXD(end) TXD(end)], [min(timeM) max(timeM)],'-.k','linewidth',1)
   


    for ii = 1:size(parametros,1)
        plot(TXD(ii),timeM(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
    end

    [p, idxs] = paretoFront([TXD,(-timeM)] ); 
    auxL = p(:,1)<0.5; 
    p(auxL,:) = [];
    idxs(auxL,:) = [];

    [~,ind] = sort(p(:,1));
    p = p(ind,:);
    idxs = idxs(ind,:);
    
    % Linha colorida:
    tipo = ['-o' cor(metodo)];
    linha(metodo) = plot(TXD(idxs),timeM(idxs),tipo,'Markersize',8,'linewidth',1.8) ;
    
    for ii = 1:size(idxs,1)

        [I] = find((TXD == TXD(idxs(ii))) & (timeM==timeM(idxs(ii))));
        I = I(1); 

        for jj = I  
            text(TXD(jj),timeM(idxs(ii))*.975,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ]);
        end    
    end
end

%%
set(axes1,'XMinorTick','on');
set(axes1,'YMinorTick','on');
box(axes1,'off');


xlim([min(TXD(idxs))*.95,max(TXD(idxs))*1.05])
ylim([min(timeM(idxs))*.95 Mmax*1.05])

grid on
hold on
title('Curva de Pareto - 30 dB')
legend([linha(1),linha(2),linha(3),linha(4)],'M�todo 1','M�todo 2','M�todo 3','M�todo 4')
xlabel('Taxa de detecção [%]')
ylabel('Tempo de exame [s]')
