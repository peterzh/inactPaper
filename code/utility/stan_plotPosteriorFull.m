function plotPosteriorFull(p)
pNames = fieldnames(p);
p = rmfield(p, pNames(find(strcmp(pNames,'z')):end));
p = rmfield(p, pNames(find(strcmp(pNames,'log_lik')):end));
pNames = fieldnames(p);


%Expand the param labels if any variable in the posterior
%struct is more than 1D
pNames_expanded = pNames;
for param = 1:length(pNames)
    numVars = size(p.(pNames{param}),2);
    pNames_expanded{param} = repmat(pNames(param),numVars,1);
end
pNames_expanded = cat(1,pNames_expanded{:});

%Plot posterior distribution
figure('color','w');
[S,AX,BigAx,H,HAx] = plotmatrix(struct2array(p),'k.');
set(AX,'box','off');
set(HAx,'box','off');
set(H,'DisplayStyle','bar', 'edgealpha', 0);
delete(AX(find(triu(ones(length(pNames_expanded)),1))));

for param = 1:length(pNames_expanded)
    xlabel(AX(end,param),pNames_expanded{param},'interpreter','none', 'Rotation', 45);
    ylabel(AX(param,1),pNames_expanded{param},'interpreter','none', 'Rotation', 45);
end
end