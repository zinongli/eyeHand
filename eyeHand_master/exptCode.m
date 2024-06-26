clear all
Screen('Preference', 'SkipSyncTests', 1);
subjectID = 'pilo';
subj = subjectID;
eyeDataDir = 'eyedata';
dateTime = clock;                %gets time for seed
rng(sum(100*dateTime) );
expName = 'practice';
session = 01;
redoCalib = 0;
eyeFile = sprintf('%s%s', subjectID, datestr(now, 'mmdd'));

nTrials = 10;



[displayInfo] = startExp(subjectID,datetime,rng);
[displayInfo] = screenVisuals(displayInfo);
if exist(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_tform.mat'],'file') && redoCalib == 0
    load(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_tform.mat']) %load calibration incase of restart
    load(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_calibration.mat'])
    
else
    [tform, calibration,startPhase] = penCalib(displayInfo);
    save(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_tform.mat'],'tform')
    save(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_calibration.mat'],'calibration')
end
mode = 0; % lift = 1, slide = 0
%%
saccadeStartSize = 50; % radius of allowable eye movement in projection pixels
reachStartSize = 20; % radius of allowable eye and hand movement in projection pixels

cursor_size = 5;
pixellength = 0.248;
Affine2d =tform.T(1:2,1:2);
[~,s,~] = svd(Affine2d);
proj2tablet = 1./mean([s(1,1),s(2,2)]);
wait = 1;
patience = 0.5;
topBuff = [0 0 displayInfo.screenXpixels displayInfo.screenAdj/2]; %black bar at top of screen
bottomBuff = [0 displayInfo.screenYpixels-displayInfo.screenAdj/2 displayInfo.screenXpixels displayInfo.screenYpixels]; %black bar at bottom of screen

%% Task Parameters
dists_n = 1;
mmPerProjPx = proj2tablet .* pixellength;
UniRandRadius = 70 .* mmPerProjPx;
edgesize = 50;
sizes_n = 6;
WindowWidth = 1024; % projector window width
yCenter = 384; % projector screen y center
rep = 5;
totalTime = 3;
block_n = 6;

distances = linspace(edgesize,WindowWidth-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,sizes_n*rep) .* mmPerProjPx;

% size granularity TBD, 5:10:55 mm for now
target_sizes = (5:10:55) ./ pixellength ./ proj2tablet;
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';

lifespan = repmat(totalTime,1,block_n);
%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
traXtotalE = [];
traYtotalE = [];
testtimes = zeros(1,10000); % 10 seconds
instruct = 'The experiment will begin.\n Please fixate at the center crosshair and move the cursor there too. \n When a target appears, please look at the target and reach it with your cursor. \n Put the pen on the tablet to start. Press any key to exit.';
%% Screen
% pscreenNumber = max(Screen('Screens'))-1;
% [window,rect] = Screen('OpenWindow', screenNumber);
% [window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]); % for testing
[cx,cy] = RectCenter(displayInfo.windowRect);
Screen('TextSize', displayInfo.window, 24);
Screen('TextColor', displayInfo.window, 255);
Screen('TextFont', displayInfo.window, 'Verdana');
Screen('FillRect', displayInfo.window, 128)
Screen('Flip', displayInfo.window);

%% Initialize eye tracker
[el, exitFlag] = rd_eyeLink('eyestart', displayInfo.window, eyeFile);
if exitFlag
    return
end

%% Calibrate eye tracker
[cal, exitFlag] = rd_eyeLink('calibrate', displayInfo.window, el);
if exitFlag
    return
end
pause(0.5)
%% Present trials
HideCursor;

while true
    DrawFormattedText(displayInfo.window,instruct,'center','center', displayInfo.blackVal);
    Screen('Flip', displayInfo.window);
    [~,~,b] = GetMouse;
    if b(1)
        break
    end
    if KbCheck
        Screen('CloseAll')
        break
    end
end

wait = 1;
flashDur = 0.06;
trialDur = 2;
tSize = 10;
tColor = [0 128 0];
cSize = 10;
cColor = [128 0 0];
block_n = 6;
t_0 = [70,70;-70,70] ./ mmPerProjPx + [cx,cy];
freqE = 3;
freqH = 5;
handPertAxis = "theta";
eyePertAxis = "rho";
p_handMag = deg2rad(8);
p_eyeMag = 0.1;
cursorHomeRadius = 50;
arrowTipRadius = 20;
backGround = 50;
data = [];
traXtotalE = [];
traYtotalE = [];
traXtotalH = [];
traYtotalH = [];
waveInd = 1;
for j = 1:block_n
    seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
    randdists = distances(seeds);
    randdists = randdists(:);
    randsizes = target_sizes(seeds);
    randsizes = randsizes(:);
    params = NaN(length(randdists),16);
    traHx = NaN(length(randdists),round(displayInfo.framerate * (wait+lifespan(j)+patience)));
    traHy = NaN(length(randdists),round(displayInfo.framerate * (wait+lifespan(j)+patience)));
    traEx = NaN(length(randdists),round(displayInfo.framerate * (wait+lifespan(j)+patience)));
    traEy = NaN(length(randdists),round(displayInfo.framerate * (wait+lifespan(j)+patience)));
    trial_n = length(randdists);
    pDuration = length(randdists).*block_n;
    trials = ones(1,trial_n);
    
    i = 0;
    DrawFormattedText(displayInfo.window,['Next Block: ',num2str(j),'/',num2str(block_n)],'center','center',displayInfo.whiteVal); % not sure how to get this centered yet
    Screen('Flip', displayInfo.window);
    pause(2);
    
    while sum(trials) > 0
        i = i+1;
        stage = 0;
        frame = 0;
        if i == trial_n + 1
            randdists = [randdists ; randdists(trials==true,:)];
            randsizes = [randsizes ; randsizes(trials==true,:)];
            params = [params ; params(trials==true,:)];
            traHx = [traHx ; traHx(trials==true,:)];
            traHy = [traHy ; traHy(trials==true,:)];
            wrong_n = sum(trials);
            origin_trial_n = trial_n;
            trial_n = trial_n + wrong_n;
            trials = zeros(1,trial_n);
            trials(1,origin_trial_n+1:end) = 1;
        end
        tar_i = t_0(rem(waveInd,2)+1,:);
        theta = atan2(tar_i(2) - cy, tar_i(1) - cx);
        rho = norm(tar_i - [cx,cy]);
        switch eyePertAxis
            case "theta"
                p_eye_i  = getPertbScalar(freqE,pDuration,waveInd,p_eyeMag);
                [p_eyeCart(1),p_eyeCart(2)] = pol2cart(theta + p_eye_i,rho);
            case "rho"
                p_eye_i  = getPertbScalar(freqE,pDuration,waveInd,p_eyeMag);
                [p_eyeCart(1),p_eyeCart(2)] = pol2cart(theta,rho.*(1+p_eye_i));
        end
        % f_eye = tar_i + p_eyeCart;
        p_eyeCart = p_eyeCart - (tar_i-[cx,cy]);
        
        params(i,1:2) = tar_i;
        params(i,3) = p_eye_i;
        params(i,5:6) = [cx,cy];
        params(i,7:8) = tar_i + p_eyeCart; % f_eye
        %         rd_eyeLink('trialstart',  displayInfo.window, {el, i, cx, cy, start_size});
        % Displays a title at the bottom of the eye tracker display
        Eyelink('Command', 'record_status_message ''Starting trial %d''', i);
        err=Eyelink('CheckRecording');
        if err~=0
            rd_eyeLink('startrecording',  displayInfo.window, el);
        end
        Eyelink('Message', 'TRIAL_START %d', i);
        Eyelink('Message', 'SYNCTIME');		% zero-plot time for EDFVIEW
        
        % staring at fixation
        % present fixation
        Screen('FillRect', displayInfo.window, backGround)
        DrawFormattedText(displayInfo.window,'+','center','center', displayInfo.blackVal);
        timeFix = Screen('Flip', displayInfo.window);
        % direct hand back to fixation
        while true
            [~,~,keyCode] = KbCheck;
            if find(keyCode) == 27 % KbName(27) = 'ESCAPE'
                Screen('CloseAll');
                ShowCursor;
                rd_eyeLink('stoprecording');
                break
            end
            [x,y,buttons] = GetMouse(displayInfo.window2);
            [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
            DrawFormattedText(displayInfo.window,'+','center','center', displayInfo.blackVal);
            % talk to eyelink to find eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            Ex = evt.gx(domEye);
            Ey = evt.gy(domEye);
            % check for blink
            if (~isempty(Ex) || ~isempty(Ey)) && norm([Ex, Ey] - [cx, cy])<saccadeStartSize
                fixation = 1;
            else
                fixation = 0;
            end
            %              % If fixating, colored. If not, black
            if norm(xy - [cx,cy]) >= cursorHomeRadius
                drawArrow2Center(displayInfo.window,cx,cy,xy(1),xy(2),arrowTipRadius,[255 0 0] .* fixation)
            else
                Screen('DrawDots', displayInfo.window, xy, cSize, [255 0 0] .* fixation,[],1);
            end
            if buttons(1) && fixation
                if norm(xy - [cx,cy]) <= saccadeStartSize % if in start area
                    break
                end
            end
            Screen('Flip', displayInfo.window);
        end
        
        
        % start the trial when the eyetracker is recording and the subject is
        % holding fixation
        %         rd_eyeLink('trialstart', displayInfo.window, {el, i, cx, cy, start_size});
        fprintf('\n\nTrial %d \n\n', i)
        % present first stimulus
        %         fixation = rd_eyeLink('fixholdcheck', displayInfo.window, {cx, cy, start_size});
        Eyelink('Message', 'WaitStarts');
        %         timeStim1(i) = Screen('Flip', displayInfo.window);
        for frame = 1:round(wait*displayInfo.framerate)
            %             [x,y,buttons] = GetMouse(displayInfo.window2);
            %         Screen('FillRect', displayInfo.window, [0.5 0.5 0.5])
            % talk to eyelink to find eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            Ex = evt.gx(domEye);
            Ey = evt.gy(domEye);
            % check for blink
            if (~isempty(Ex) || ~isempty(Ey)) && norm([Ex, Ey] - [cx, cy])<saccadeStartSize
                fixation = 1;
            else
                fixation = 0;
            end
            %             DrawFormattedText(displayInfo.window, '+', 'center', 'center', [0 0 0]);
            if ~fixation
                DrawFormattedText(displayInfo.window,'+','center','center', displayInfo.whiteVal);
            else
                DrawFormattedText(displayInfo.window,'+','center','center', displayInfo.blackVal);
                %                 Screen('DrawDots', displayInfo.window, xy, 5, [0 0 0],[],1);
            end
            %             [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
            Screen('Flip',displayInfo.window);
            [~,~,keyCode] = KbCheck;
            if find(keyCode) == 27 % KbName(27) = 'ESCAPE'
                Screen('CloseAll');
                ShowCursor;
                rd_eyeLink('stoprecording');
                break
            end
        end
        %             fixation_buffer_check = 1;
        onset_recorded = 0;
        %             fixation_cache = ones(1,10);
        saccadeOnset = 0;
        Eyelink('Message', 'TargetAppear');
        frame = 1;
        while ~saccadeOnset
            frame = frame + 1;
            DrawFormattedText(displayInfo.window,'+','center','center', displayInfo.blackVal);
            Screen('DrawDots',displayInfo.window, params(i,1:2), tSize,tColor,[],1);
            Screen('Flip',displayInfo.window);
            % talk to eyelink to find eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            Ex = evt.gx(domEye);
            Ey = evt.gy(domEye);
            % check for blink
            if (~isempty(Ex) || ~isempty(Ey))
                if norm([Ex, Ey] - [cx, cy])<saccadeStartSize
                    fixation = 1;
                else
                    saccadeOnset = 1;
                end
            else
                fixation = 0;
            end
            %                 fixation_cache = [fixation_cache , fixation];
            %                 fixation_buffer_check = sum(fixation_cache(end-10:end)) > 0;
            cache = [x,y];
            [x,y,buttons] = GetMouse(displayInfo.window2);
            locdiff = norm(cache - [x,y]);
            [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
            if norm(xy - [cx, cy]) >= reachStartSize && ~onset_recorded
                Eyelink('Message', 'ReachOnset');
                abs_onset_time = GetSecs;
                params(i,9) = abs_onset_time;
                onset_recorded = 1;
                break
            end
        end
        Eyelink('Message', 'SaccadeOnset');
        frameTelo = frame;
        % saccade onsets
        % flash
        for frame = 1:round(flashDur * displayInfo.framerate)
            Screen('FillRect', displayInfo.window, 200);
            Screen('Flip', displayInfo.window);
            cache = [x,y];
            
            % talk to eyelink to find eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            Ex = evt.gx(domEye);
            Ey = evt.gy(domEye);
            if isempty(Ex) || isempty(Ey)
                [Ex , Ey] = deal(NaN);
            end
            [x,y,buttons] = GetMouse(displayInfo.window2);
            traEx(i,frame+frameTelo) = Ex;
            traEy(i,frame+frameTelo) = Ey;
            traHx(i,frame+frameTelo) = x;
            traHy(i,frame+frameTelo) = y;
            locdiff = norm(cache - [x,y]);
            if norm(xy - [cx, cy]) >= reachStartSize && ~onset_recorded
                Eyelink('Message', 'ReachOnset');
                abs_onset_time = GetSecs;
                params(i,9) = abs_onset_time;
                onset_recorded = 1;
            end
        end
        frameTelo = frame;
        Screen('FillRect', displayInfo.window, backGround);
        
        % jump the target
        [x,y,~] = GetMouse(displayInfo.window2); % for speed cache
        for frame = 1: round(displayInfo.framerate * trialDur)
            
            Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
            Screen('Flip', displayInfo.window);
            
            cache = [x,y];
            traHx(i,frame+frameTelo) = x;
            traHy(i,frame+frameTelo) = y;
            % talk to eyelink to find eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            Ex = evt.gx(domEye);
            Ey = evt.gy(domEye);
            if isempty(Ex) || isempty(Ey)
                [Ex , Ey] = deal(NaN);
            end
            traEx(i,frame+frameTelo) = Ex;
            traEy(i,frame+frameTelo) = Ey;
            [x,y,buttons] = GetMouse(displayInfo.window2);
            locdiff = norm(cache - [x,y]);
            [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
            % reach onset
            if norm(xy - [cx, cy]) >= reachStartSize
                theta = atan2(xy(2) - cy, xy(1) - cx);
                rho = norm(xy - [cx,cy]);
                switch handPertAxis
                    case "theta"
                        p_hand_i  = getPertbScalar(freqH,pDuration,waveInd,p_handMag);
                        [p_handCart(1),p_handCart(2)] = pol2cart(theta + p_hand_i,rho);
                    case "rho"
                        p_hand_i  = getPertbScalar(freqH,pDuration,waveInd,p_handMag);
                        [p_handCart(1),p_handCart(2)] = pol2cart(theta,rho .* (1+p_hand_i));
                end
                params(i,4) = p_hand_i;
                p_handCart = p_handCart - (xy-[cx,cy]);
                Screen('DrawDots', displayInfo.window, xy + p_eyeCart + p_handCart, cSize, [128 0 0],[],1);
                Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
                Screen('Flip', displayInfo.window);
                if buttons(1) ==0
                    DrawFormattedText(displayInfo.window,'Stylus Lifted','center','center', displayInfo.blackVal);
                    Screen('Flip', displayInfo.window);
                    trials(i) = 1;
                    pause(1)
                    break
                elseif ~onset_recorded
                    Eyelink('Message', 'ReachOnset');
                    abs_onset_time = GetSecs;
                    params(i,9) = abs_onset_time;
                    onset_recorded = 1;
                end
                
                if locdiff <= speedthreshold/displayInfo.framerate
                    params(i,10:11) = xy; % e_hand
                    params(i,12:13) = xy + p_eyeCart + p_handCart; % f_hand
                    end_t = frame / displayInfo.framerate;
                    rest_of_trial = trialDur - end_t;
                    trials(i) = 0;
                    waveInd = waveInd + 1;
                    Eyelink('Message', 'ReachFinish');
                    %                         Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
                    %                         Screen('Flip', displayInfo.window);
                    %                         pause(delayB4Fhand);
                    Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
                    Screen('DrawDots', displayInfo.window, xy + p_eyeCart + p_handCart, cSize, [128 0 0],[],1);
                    Screen('Flip', displayInfo.window);
                    pause(rest_of_trial);
                    break
                end
            end
            %                 Screen('Flip', displayInfo.window);
        end
    end
    data = [data;params];
    save(['data_eyeHand\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_rawtotal.mat'],'data');
    xTrajTelomereE = NaN(size(traEx,1),size(traXtotalE,2));
    xTrajTelomereE(1:size(traEx,1),1:size(traEx,2)) = traEx;
    traXtotalE = [traXtotalE;xTrajTelomereE];
    yTrajTelomereE = NaN(size(traEy,1),size(traYtotalE,2));
    yTrajTelomereE(1:size(traEy,1),1:size(traEy,2)) = traEy;
    traYtotalE = [traYtotalE;yTrajTelomereE];
    save(['data_eyeHand\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_traEXtotal.mat'],'traXtotalE')
    save(['data_eyeHand\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_traEYtotal.mat'],'traYtotalE')
    
    xTrajTelomereH = NaN(size(traHx,1),size(traXtotalH,2));
    xTrajTelomereH(1:size(traHx,1),1:size(traHx,2)) = traHx;
    traXtotalH = [traXtotalH;xTrajTelomereH];
    yTrajTelomereH = NaN(size(traHy,1),size(traYtotalH,2));
    yTrajTelomereH(1:size(traHy,1),1:size(traHy,2)) = traHy;
    traYtotalH = [traYtotalH;yTrajTelomereH];
    save(['data_eyeHand\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_traHXtotal.mat'],'traXtotalH')
    save(['data_eyeHand\' subj '\' subj '_' expName '_S' num2str(session) '_' date,'_traHYtotal.mat'],'traYtotalH')

    while true
        DrawFormattedText(displayInfo.window,'Block finished. Press any key to proceed to next block.','center','center', displayInfo.whiteVal); % not sure how to get this centered yet
        Screen('Flip', displayInfo.window);
        if KbCheck
            break
        end
    end
end


%% Save the eye data and shut down the eye tracker
if ~exist(eyeDataDir,'dir')
    mkdir(eyeDataDir)
end
rd_eyeLink('eyestop', displayInfo.window, {eyeFile, eyeDataDir});

%% Close screen
Screen('CloseAll')
ShowCursor;

