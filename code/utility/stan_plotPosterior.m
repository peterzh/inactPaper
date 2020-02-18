function plotPosterior(posterior)

pNames = fieldnames(posterior);
posterior = rmfield(posterior,pNames(contains(pNames,{'log_lik','logOdds','rho','z','corr','log_pRpL'})) );
pNames = fieldnames(posterior);

%remove any variable with more than 2 dims
dim = structfun(@(f) ndims(f), posterior, 'uni', 1);
posterior = rmfield(posterior,pNames(dim>2) );
pNames = fieldnames(posterior);


%Compute means
M = structfun(@(f) mean(f,1), posterior, 'uni', 0);

%Compute quantiles
Q1 = structfun(@(f) quantile(f, [0.2 0.8], 1), posterior, 'uni', 0);
Q2 = structfun(@(f) quantile(f, [0.025 0.975], 1), posterior, 'uni', 0);

figure('color','w');
ha=axes; hold(ha,'on');
ii=1;
ytickLabels = {};
for param = 1:length(pNames)
    
    thisP_means = M.(pNames{param});
    thisP_Q1 = Q1.(pNames{param});
    thisP_Q2 = Q2.(pNames{param});
    
    for vals = 1:numel(thisP_means)
        plot(ha, thisP_means(vals), ii, '.','markersize',20,'Color',[1 1 1]*0.3);
        line(ha, thisP_Q1(:,vals), [1 1]*ii,'Color',[1 1 1]*0.3);
        plot(ha, thisP_Q2(:,vals), ii, '+','markersize',8,'Color',[1 1 1]*0.8);
        
        ytickLabels{ii}=[pNames{param} ' ' num2str(vals)];
        
        ii=ii+1;
    end
end
xlabel(ha,'Value');
set(ha,'ytick',1:(ii-1),'YGrid','on','GridLineStyle',':','YTickLabel',ytickLabels,'TickLabelInterpreter','none')
line(ha,[0 0],ylim(ha),'Color',[1 1 1]*0.8,'linestyle','--')

end