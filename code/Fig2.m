
%% Load and average

decoding_all = cell(height(sessionList),1);
for sess = 1:height(sessionList)
    eRef = sessionList.expRef{sess};
    fprintf('Session %d %s\n',sess,eRef);
    
    activityFile = [ '../data/widefield/eventAlignedActivity/' eRef '.mat'];
    
    a = load(activityFile);
    
    decoding_all{sess} = a.stimAlignedDecoding;
    
end
decoding_all = cat(5,decoding_all{:});

decoding_avg = nanmean(decoding_all,5);
