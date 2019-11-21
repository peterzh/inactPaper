

function [cweA, cwtA, stimOn, moveTimeAbs, cR, cL, hasR, hasL] = ...
    loadSelectedTrials(mouseName, thisDate)


[~, cweA, cwtA, moveData, lickTimes, passiveStim] = alf.loadCWAlf(...
    mouseName, thisDate, []);

% choose trials/events
winCountMove = [0.1 0.4];
winNoMove = [-0.05 0.4];
timeBinSize = 0.005;

% trying to fix edge effects by rounding each stim event to the nearest
% time bin... :/ 
% cwtA.stimOn = round(cwtA.stimOn/timeBinSize)*timeBinSize;

stimOn = cwtA.stimOn;
cL = cweA.contrastLeft; 
cR = cweA.contrastRight;
choice = cweA.choice;
fb = cweA.feedback; 
inclT = cweA.inclTrials;

ucl = unique(cL);
ucr = unique(cR);

[moveTimeUse, firstMoveType, firstMoveTime, hasNoMoveInWin] = ...
                findMoveTimes(cweA, cwtA, moveData, winCountMove, winNoMove);

hasLeftMove = firstMoveTime<winCountMove(2) & firstMoveTime>winCountMove(1) & firstMoveType==1;
hasRightMove = firstMoveTime<winCountMove(2) & firstMoveTime>winCountMove(1) & firstMoveType==2;
moveTimeAbs = moveTimeUse+cwtA.stimOn;

% also fixing these times with rounding
moveTimeAbs = round(moveTimeAbs/timeBinSize)*timeBinSize;

inclT = inclT & (hasLeftMove | hasRightMove | hasNoMoveInWin);           

% moveTimeUse = moveTimeUse+cwtA.stimOn;

cweA = horzcat(cweA, table(inclT, hasLeftMove, hasRightMove, hasNoMoveInWin, firstMoveType));
% cweA = addvars(cweA, inclT, hasLeftMove, hasRightMove, hasNoMoveInWin, firstMoveType);
% cwtA = addvars(cwtA, moveTimeUse, moveTimeAbs);
cwtA = horzcat(cwtA, table(moveTimeUse, moveTimeAbs));
stimOn = cwtA.stimOn(cweA.inclT);
moveTimeAbs = cwtA.moveTimeAbs(cweA.inclT);
cR = cweA.contrastRight(cweA.inclT);
cL = cweA.contrastLeft(cweA.inclT);
hasR = cweA.hasRightMove(cweA.inclT);
hasL = cweA.hasLeftMove(cweA.inclT);