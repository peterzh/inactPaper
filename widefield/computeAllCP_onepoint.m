

function [cp, dp, cpM, dpM] = computeAllCP_onepoint(r, stimWin, moveWin, nShuf)



% nShuf = 2000;

addRand = true; % to break ties, with the idea that it will result in 5% of neurons
% actually being significant at p=0.05 level. 

allAcr = vertcat(r.acr); allR = vertcat(r.r);

% last dim is [value, p-val, shuffleMean, max shuffle p, min shuffle p, first shuffle val, first shuffle p]
cp = zeros(numel(allAcr), 7); 
dp = zeros(numel(allAcr), 7);
cpM = zeros(numel(allAcr), 7);
dpM = zeros(numel(allAcr), 7);


% for each recording
nidx = 1;
for n = 1:numel(r)
    fprintf(1, '%d/%d\n', n, numel(r));
    tic
    cweA = r(n).cweA; cwtA = r(n).cwtA;
    sp = r(n).sp;
    
    % compute binned spikes anew here
    fprintf(1, '  compute bs\n');
    maxT = cwtA.feedbackTime(end)+5;
    theseST = sp.st;
    clu = sp.clu;
    cids = find(r(n).inclU)-1;
    theseClu = double(clu(theseST<maxT & theseST>0 & ismember(clu,cids)));
    theseST = theseST(theseST<maxT & theseST>0 & ismember(clu,cids));     
    
    
    stimOn = cwtA.stimOn;    
    moveTimeAbs = cwtA.moveTimeAbs;
    inclT = cweA.inclT;
    hasL = cweA.hasLeftMove; hasR = cweA.hasRightMove; hasNG = cweA.hasNoMoveInWin;
    hasL = hasL(inclT); hasR = hasR(inclT); hasNG = hasNG(inclT);
    
    trialChoice = false(size(hasL)); 
    trialChoice(hasR) = true;
    trialChoice = trialChoice(hasL|hasR);
    
    trialDetect = false(size(hasL)); 
    trialDetect(hasL|hasR) = true;
    
    
    
    conL = cweA.contrastLeft(inclT); conR = cweA.contrastRight(inclT);
    [~,~,cLnorm] = unique(conL); [~,~,cRnorm] = unique(conR);
    trialCondition = sub2ind([4 4], cLnorm, cRnorm);
         
    q = crosstab(trialCondition(hasL|hasR), trialChoice);
    nComp = sum(q(:,1).*q(:,2));
    if nComp<10
        excludeThisRec = true;
    else
        excludeThisRec = false;
    end
    
    % for each neuron, get spike counts    
    if ~isempty(stimWin)
        fprintf(1, '  compute stimRates\n');
        [stimRates, ~] = makeAllBA(theseST, theseClu, stimOn(inclT), diff(stimWin), stimWin);
        stimRates = stimRates'; % trials x neurons
    end
    
    if ~isempty(moveWin)
        fprintf(1, '  compute moveRates\n');
        [moveRates, ~] = makeAllBA(theseST, theseClu, moveTimeAbs(inclT), diff(moveWin), moveWin);
        moveRates = moveRates';
    end
    
    if addRand
        % this resolves ties in the spike count randomly
        if ~isempty(stimWin)
            stimRates = stimRates+rand(size(stimRates))/1e6;
        end
        if ~isempty(moveWin)
            moveRates = moveRates+rand(size(moveRates))/1e6;
        end
    end
    
    % create shuffles
    
    tch = trialChoice; tcond = trialCondition(hasL|hasR);
    uCond = unique(tcond);
    for c = 1:numel(uCond)
        inclT = tcond==uCond(c);
        
        chA = tch & inclT;
        nA = sum(chA);
        chB = ~tch & inclT;
        nB = sum(chB);
        q = arrayfun(@(x)randperm(nA+nB,nA), 1:nShuf, 'uni', false);
        shufLabelsC{c} = vertcat(q{:})';
        shufLabelsC{c} = [(1:nA)' shufLabelsC{c}];
    end
    
    
    tch = trialDetect; tcond = trialCondition;
    uCond = unique(tcond);
    for c = 1:numel(uCond)
        inclT = tcond==uCond(c);
        
        chA = tch & inclT;
        nA = sum(chA);
        chB = ~tch & inclT;
        nB = sum(chB);
        q = arrayfun(@(x)randperm(nA+nB,nA), 1:nShuf, 'uni', false);
        shufLabelsD{c} = vertcat(q{:})';
        shufLabelsD{c} = [(1:nA)' shufLabelsD{c}];
    end
    
    % for each time point

    fprintf(1, '  compute cp\n');
    startNidx = nidx;
    
    for nn = 1:numel(cids)
        
        if ~excludeThisRec
            
            if ~isempty(stimWin)
                
                spikeCounts = squeeze(stimRates(:,nn));
                
                % choice prob
                [cc,p,cpSummary] = choiceProbShuf(spikeCounts(hasL|hasR), trialChoice, ...
                    trialCondition(hasL|hasR), shufLabelsC);
                
                cp(nidx,:) = cpSummary;
                
                % detect prob
                [cc,p,cpSummary] = choiceProbShuf(spikeCounts, trialDetect, ...
                    trialCondition, shufLabelsD);
                
                dp(nidx,:) = cpSummary;
                
            end
            
            if ~isempty(moveWin)
                
                spikeCounts = squeeze(moveRates(:,nn));
                
                % choice prob
                [cc,p,cpSummary] = choiceProbShuf(spikeCounts(hasL|hasR), trialChoice, ...
                    trialCondition(hasL|hasR), shufLabelsC);
                
                cpM(nidx,:) = cpSummary;
                
                % detect prob
                [cc,p,cpSummary] = choiceProbShuf(spikeCounts, trialDetect, ...
                    trialCondition, shufLabelsD);
                
                dpM(nidx,:) = cpSummary;
            end
            
        end
        nidx = nidx+1;
    end
    toc
    
end

% correct p-values that are 0 and 1 based on nshuf
for q = [2 4 5 7]
    p1 = cp(:,q); p1(p1==1) = (nShuf-0.5)/nShuf;
    cp(:,q) = p1; 
    
    p1 = cpM(:,q); p1(p1==1) = (nShuf-0.5)/nShuf;
    cpM(:,q) = p1; 

    p1 = dp(:,q); p1(p1==1) = (nShuf-0.5)/nShuf;
    dp(:,q) = p1; 

    p1 = dpM(:,q); p1(p1==1) = (nShuf-0.5)/nShuf;
    dpM(:,q) = p1; 
end