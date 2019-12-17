%% 52 coord dataset
clear all;
expRefs = readtable('../data/inactivation/sessionList_52Coord.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

%Load standard coordset
load('../data/inactivation/26CoordSet.mat','coordSet');
coordSet = [coordSet; -coordSet(:,1), coordSet(:,2)];

D = table;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadData(expRefs{session});

    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    dd.subjectID = ones(length(dd.response),1)*subjID(session);
    
    %Add extra inactivation labels
    dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
    dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
    dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
    dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
    
    %Remove trials with inactivation on the non-standard
    %coordinate set
    keep = ones(size(dd.response));
    for tr = 1:length(dd.response)
        if dd.laserType(tr)>0
            if any(sum(dd.laserCoord(tr,:) == coordSet,2)==2)
                keep(tr)=1;
            else
                keep(tr)=0;
            end
        end
    end
    dd = structfun(@(f) f(find(keep),:), dd, 'uni', 0);
    
    %If there aren't any laser trials
    if ~any(dd.laserCoord(~isnan(dd.laserCoord(:,1)),1) ~= 0)
        %                         keyboard;
    end
    
    if any(dd.laserType>0)
        D = vertcat(D,struct2table(dd));
        ii=ii+1;
    end
end
D = D(D.repeatNum==1,:);

D.laserIdx = zeros(size(D.response));
%Define laser Idx
for i = 1:size(coordSet,1)
    id = D.laserCoord(:,1) == coordSet(i,1) & D.laserCoord(:,2) == coordSet(i,2);
    D.laserIdx(id) = i;
end

wheel_stimulusOn_timesteps = D.wheel_stimulusOn_timesteps(1,:);
D.wheel_stimulusOn_timesteps = [];
D.wheel_startMove_timesteps=[];
D.wheel_startMove=[];

%Save the data
save('../data/inactivation/Inactivation_52Coord.mat','D','wheel_stimulusOn_timesteps','-v7.3');

%% Pulse inactivation dataaset
clear all;
expRefs = readtable('../data/inactivation/sessionList_Pulse.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = table;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadData(expRefs{session});

    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    dd.subjectID = ones(length(dd.response),1)*subjID(session);
    
    %Add extra inactivation labels
    dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
    dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
    dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
    dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
    
%     %add marker for whether the choice is contra or ipsi to the inactivated side is contra or ipsi to stimulus side
%     dd.ipsiversive = ( contains(string(dd.laserRegion),'Left') & (dd.response==1) ) | ...
%                ( contains(string(dd.laserRegion),'Right') & (dd.response==2) );

    if any(dd.laserType>0)
        D = vertcat(D,struct2table(dd));
        ii=ii+1;
    end
end
D = D(D.repeatNum==1,:);
% 
% %define stimulus conditions where one side is higher than the other
% D.CL_gt_CR = D.stimulus(:,1) > D.stimulus(:,2);
% D.CR_gt_CL = D.stimulus(:,2) > D.stimulus(:,1);
% 
% %prepare laser onset times for plotting later
% D.laserOnset = D.laserOnset + randn(length(D.laserOnset),1)/100000;
% [~,sortIdx]=sort(D.laserOnset);
% D = D(sortIdx,:);

%Remove trials with pre-stimulus wheel movement
keepIdx = zeros(size(D.response));
for sess = 1:max(D.sessionID)
    times = D.wheel_stimulusOn_timesteps(D.sessionID==sess,:);
    times = times(1,:);
    wheel = D.wheel_stimulusOn(D.sessionID==sess,:);
    choice = D.response(D.sessionID==sess);
    
    %Whiten
    wheel_whiten =  wheel/std(wheel(:,end));
    
    %Threshold
    wheel_whiten_window = mean( abs(wheel_whiten(:,-0.150 < times & times < 0.05 & times ~= 0)) , 2);
    keepIdx(D.sessionID==sess) = wheel_whiten_window < 0.05;
end
fprintf('Removing %0.2f%% of trials becuase of pre-stim wheel movements\n',100*mean(keepIdx==0));
D = D(keepIdx==1,:);

wheel_stimulusOn_timesteps = D.wheel_stimulusOn_timesteps(1,:);
D.wheel_stimulusOn_timesteps = [];
D.wheel_startMove_timesteps=[];
D.wheel_startMove=[];

%Save the data 
save('../data/inactivation/Inactivation_Pulse.mat','D','wheel_stimulusOn_timesteps','-v7.3');

%% Higher power inactivation dataset
clear all;
expRefs = readtable('../data/inactivation/sessionList_HigherPower.csv','FileType','text','Delimiter',',');
expRefs = expRefs.expRef;
mouseName = cellfun(@(e)e{3},cellfun(@(s) strsplit(s,'_'), expRefs,'uni',0),'uni',0);
[mice,~,subjID] = unique(mouseName);

D = table;
ii=1;
for session = 1:length(expRefs)
    [dd,meta] = loadData(expRefs{session});
    dd = structfun(@(x)(x(6:(end-14),:)),dd,'uni',0); %trim first 5 trials and last 15
    dd.sessionID = ones(length(dd.response),1)*ii;
    dd.subjectID = ones(length(dd.response),1)*subjID(session);
    
    %Add extra inactivation labels
    dd.laserRegion(dd.laserCoord(:,1)==0.5 & dd.laserCoord(:,2)==-2) = 'RightRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-0.5 & dd.laserCoord(:,2)==-2) = 'LeftRSP';
    dd.laserRegion(dd.laserCoord(:,1)==-3.5 & dd.laserCoord(:,2)==-0.5) = 'LeftS1';
    dd.laserRegion(dd.laserCoord(:,1)==3.5 & dd.laserCoord(:,2)==-0.5) = 'RightS1';
    dd.laserRegion(dd.laserCoord(:,1)==1 & dd.laserCoord(:,2)==-0.5) = 'RightM1';
    dd.laserRegion(dd.laserCoord(:,1)==-1 & dd.laserCoord(:,2)==-0.5) = 'LeftM1';
    
    if any(dd.laserType>0)
        D = vertcat(D,struct2table(dd));
        ii=ii+1;
    end
end
D = D(D.repeatNum==1,:);

coordSet = unique(D.laserCoord(D.laserType==1,:),'rows');
D.laserIdx = zeros(size(D.response));
%Define laser Idx
for i = 1:size(coordSet,1)
    id = D.laserCoord(:,1) == coordSet(i,1) & D.laserCoord(:,2) == coordSet(i,2);
    D.laserIdx(id) = i;
end

%Remove trials with pre-stimulus wheel movement
keepIdx = zeros(size(D.response));
for sess = 1:max(D.sessionID)
    times = D.wheel_stimulusOn_timesteps(D.sessionID==sess,:);
    times = times(1,:);
    wheel = D.wheel_stimulusOn(D.sessionID==sess,:);
    choice = D.response(D.sessionID==sess);
    
    %Whiten
    wheel_whiten =  wheel/std(wheel(:,end));
    
    %Threshold
    wheel_whiten_window = mean( abs(wheel_whiten(:,-0.150 < times & times < 0.05 & times ~= 0)) , 2);
    keepIdx(D.sessionID==sess) = wheel_whiten_window < 0.05;
end
fprintf('Removing %0.2f%% of trials becuase of pre-stim wheel movements\n',100*mean(keepIdx==0));
D = D(keepIdx==1,:);

wheel_stimulusOn_timesteps = D.wheel_stimulusOn_timesteps(1,:);
D.wheel_stimulusOn_timesteps = [];
D.wheel_startMove_timesteps=[];
D.wheel_startMove=[];

%Save the data 
save('../data/inactivation/Inactivation_HigherPower.mat','D','wheel_stimulusOn_timesteps','-v7.3');