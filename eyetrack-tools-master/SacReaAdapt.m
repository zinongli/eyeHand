% rd_eyeTrackTester.m
%
% This is designed to be a very simple "experiment" that can be used to
% test and/or demonstrate the use of rd_eyeLink.m.
Screen('Preference', 'SkipSyncTests', 1);
% cd('C:\Users\labadmin\Documents\ZinongGitHub\eyeHand\eyetrack-tools-master');
subjectID = 'pilo';
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
if exist(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_tform.mat']) && redoCalib == 0
    load(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_tform.mat']) %load calibration incase of restart
    load(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_calibration.mat'])
    
else
    [tform, calibration,startPhase] = penCalib(displayInfo);
    save(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_tform.mat'],'tform')
    save(['data_eyeHand\' subjectID '\' subjectID '_' expName '_S' num2str(session) '_' date,'_calibration.mat'],'calibration')
end
mode = 0; % lift = 1, slide = 0
%%
start_size = 20; % radius of allowable eye and hand movement in projection pixels
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
dists_n = 3;
proj2mm = proj2tablet .* pixellength;
UniRandRadius = 70 .* proj2mm;
edgesize = 50;
sizes_n = 6;
WindowWidth = 1024; % projector window width
yCenter = 384; % projector screen y center
rep = 10;
totalTime = 3;
block_n = 6;

distances = linspace(edgesize,WindowWidth-edgesize,dists_n+2)-edgesize;
distances = repmat(distances(2:end-1),1,sizes_n*rep) .* proj2mm;

% size granularity TBD, 5:10:55 mm for now
target_sizes = (5:10:55) ./ pixellength ./ proj2tablet;
target_sizes = repmat(target_sizes,1,dists_n*rep);
target_sizes = target_sizes';
target_sizes = target_sizes(:)';

lifespan = repmat(totalTime,1,block_n);
%% Trial
speedthreshold = 10; % pixel per second, equals to 2.48 mm/s
data = [];
traXtotal = [];
traYtotal = [];
testtimes = zeros(1,10000); % 10 seconds
instruct = 'The experiment will begin. Put the pen on the tablet to start. Press any key to exit.';
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

flashDur = 0.1;
trialDur = 3;
tSize = 10;
tColor = [0 128 0];
cSize = 10;
cColor = [128 0 0];
block_n = 6;
t_0 = [400,600];
p_eye = [10,10];
p_hand = [-10,10];

for j = 1:block_n
    seeds = [randperm(size(distances,2)), randperm(size(distances,2))];
    randdists = distances(seeds);
    randdists = randdists(:);
    randsizes = target_sizes(seeds);
    randsizes = randsizes(:);
    params = NaN(length(randdists),16);
    trax = NaN(length(randdists),round(displayInfo.framerate * (wait+lifespan(j)+patience)));
    tray = NaN(length(randdists),round(displayInfo.framerate * (wait+lifespan(j)+patience)));
    trial_n = length(randdists);
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
            trax = [trax ; trax(trials==true,:)];
            tray = [tray ; tray(trials==true,:)];
            wrong_n = sum(trials);
            origin_trial_n = trial_n;
            trial_n = trial_n + wrong_n;
            trials = zeros(1,trial_n);
            trials(1,origin_trial_n+1:end) = 1;
        end
        params(i,1:2) = t_0;
        params(i,3:4) = p_eye;
        params(i,5:6) = p_hand;
        params(i,7:8) = params(i,1:2) + params(i,3:4); % f_eye
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
        Screen('FillRect', displayInfo.window, 100)
        Screen('DrawDots',displayInfo.window, [cx, cy], start_size,[255 255 255],[],1);
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
            Screen('DrawDots',displayInfo.window, [cx, cy], start_size,[255 255 255],[],1);
            % talk to eyelink to find eye position
            evt = Eyelink('newestfloatsample');
            domEye = find(evt.gx ~= -32768);
            Ex = evt.gx(domEye);
            Ey = evt.gy(domEye);
            % check for blink
            if isempty(Ex) || isempty(Ey)
                fixation = 0;
            else
                % check eye position
                if sqrt((Ex-cx)^2+(Ey-cy)^2)<start_size
                    fixation = 1;
                else
                    fixation = 0;
                end
            end
            Screen('DrawDots', displayInfo.window, xy, cSize, [255 0 0] .* fixation,[],1); % If fixating, colored. If not, black
            if buttons(1)
                if norm(xy - [cx,cy]) <= start_size % if in start area
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
        Eyelink('Message', 'EVENT_T_0');
%         timeStim1(i) = Screen('Flip', displayInfo.window);
            for t = 1:round(wait*displayInfo.framerate)
                %             [x,y,buttons] = GetMouse(displayInfo.window2);
                %         Screen('FillRect', displayInfo.window, [0.5 0.5 0.5])
                % talk to eyelink to find eye position
                evt = Eyelink('newestfloatsample');
                domEye = find(evt.gx ~= -32768);
                Ex = evt.gx(domEye);
                Ey = evt.gy(domEye);
                % check for blink
                if isempty(Ex) || isempty(Ey)
                    fixation = 0;
                else
                    % check eye position
                    if sqrt((Ex-cx)^2+(Ey-cy)^2)<start_size
                        fixation = 1;
                    else
                        fixation = 0;
                    end
                end
                %             DrawFormattedText(displayInfo.window, '+', 'center', 'center', [0 0 0]);
                if ~fixation
                    Screen('DrawDots',displayInfo.window, [cx, cy], start_size,[255 255 255],[],1);
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
            fixation_buffer_check = 1;
            onset_recorded = 0;
            fixation_cache = ones(1,10);
            while fixation_buffer_check
                DrawFormattedText(displayInfo.window,'+','center','center', displayInfo.blackVal);
                Screen('DrawDots',displayInfo.window, [400 600], tSize,[0 128 0],[],1);
                Screen('Flip',displayInfo.window);
                % talk to eyelink to find eye position
                evt = Eyelink('newestfloatsample');
                domEye = find(evt.gx ~= -32768);
                Ex = evt.gx(domEye);
                Ey = evt.gy(domEye);
                % check for blink
                if isempty(Ex) || isempty(Ey)
                    fixation = 0;
                else
                    % check eye position
                    if sqrt((Ex-cx)^2+(Ey-cy)^2)<start_size
                        fixation = 1;
                    else
                        fixation = 0;
                    end
                end
                fixation_cache = [fixation_cache , fixation];
                fixation_buffer_check = sum(fixation_cache(end-10:end)) > 0;
                cache = [x,y];
                [x,y,buttons] = GetMouse(displayInfo.window2);
                locdiff = norm(cache - [x,y]);
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                if norm(xy - [cx, cy]) >= start_size && ~onset_recorded
                    Eyelink('Message', 'Reach_Onset');
                    abs_onset_time = GetSecs;
                    params(i,9) = abs_onset_time;
                    onset_recorded = 1;
                    break
                end
            end
            Eyelink('Message', 'Saccade_Onset');
            % saccade onsets
            % flash
            for frame = 1:round(flashDur * displayInfo.framerate)
                Screen('FillRect', displayInfo.window, 200);
                Screen('Flip', displayInfo.window);
                
                cache = [x,y];
                trax(i,frame) = x;
                tray(i,frame) = y;
                [x,y,buttons] = GetMouse(displayInfo.window2);
                locdiff = norm(cache - [x,y]);
                if norm(xy - [cx, cy]) >= start_size && ~onset_recorded
                    Eyelink('Message', 'Reach_Onset');
                    abs_onset_time = GetSecs;
                    params(i,9) = abs_onset_time;
                    onset_recorded = 1;
                end
            end
            frameTelo = frame;
            Screen('FillRect', displayInfo.window, 100);
            
            % jump the target
            [x,y,~] = GetMouse(displayInfo.window2); % for speed cache
            for frame = 1: round(displayInfo.framerate * trialDur)
                
                Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
                Screen('Flip', displayInfo.window);

                cache = [x,y];
                trax(i,frame+frameTelo) = x;
                tray(i,frame+frameTelo) = y;
                [x,y,buttons] = GetMouse(displayInfo.window2);
                locdiff = norm(cache - [x,y]);
                [xy(1), xy(2)]  = transformPointsForward(tform,x,y);
                % reach onset
                if norm(xy - [cx, cy]) >= start_size
                    Screen('DrawDots', displayInfo.window, xy + params(i,3:4) + params(i,5:6), cSize, [128 0 0],[],1);
                    Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
                    Screen('Flip', displayInfo.window);
                    if buttons(1) ==0
                        DrawFormattedText(displayInfo.window,'Stylus Lifted!','center','center', displayInfo.blackVal);
                        Screen('Flip', displayInfo.window);
                        trials(i) = 1;
                        pause(1)
                        break
                    elseif ~onset_recorded
                        Eyelink('Message', 'Reach_Onset');
                        abs_onset_time = GetSecs;
                        params(i,9) = abs_onset_time;
                        onset_recorded = 1;
                    end
                    
                    if locdiff <= speedthreshold/displayInfo.framerate
                        params(i,10:11) = xy; % e_hand
                        params(i,12:13) = xy + params(i,3:4) + params(i,5:6); % f_hand
                        end_t = frame / displayInfo.framerate;
                        rest_of_trial = trialDur - end_t;
                        trials(i) = 0;
%                         Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
%                         Screen('Flip', displayInfo.window);
%                         pause(delayB4Fhand);
                        Screen('DrawDots',displayInfo.window, params(i,7:8), tSize,[0 128 0],[],1);
                        Screen('DrawDots', displayInfo.window, xy + params(i,3:4) + params(i,5:6), cSize, [128 0 0],[],1);
                        Screen('Flip', displayInfo.window);
                        pause(rest_of_trial);
                        break
                    end
                end
%                 Screen('Flip', displayInfo.window);
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
Showcursor;

