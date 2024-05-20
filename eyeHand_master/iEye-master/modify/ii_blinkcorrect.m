function [ii_data,ii_cfg] = ii_blinkcorrect(ii_data,ii_cfg,chan, pchan, pval, pri, fol,varargin)
%ii_blinkcorrect Cleans data by identifying and removing blinks
%   [ii_data,ii_cfg] =
%   ii_blinkcorrect(ii_data,ii_cfg,chan,pchan,pval,pri,fol) looks for
%   blinks in pchan (values less than pval), removes pri ms of samples before and
%   fol ms of samples after pchan dips below pval, and sets all channels in chan
%   to NaN.
%
%   ii_blinkcorrect(...,'prctile') uses a percentile-based threshold,
%   whereby pval is assumed to be an input to prctile (0-100). Threshold is
%   set at prctile(ii_data.(pchan)(ii_data.(pchan)~=0),pval) - so only
%   detected pupil values used. ~1.5 works well.
%
% ii_data.(chan) at identified blinks will be NaN
% ii_cfg.blinks is 'cursel'-style list of on/offsets of each blink
% ii_cfg.blinkvec is binary mask of blink samples (including pri/fol)
% 
% Inputs (all required):
%   chan:  channel(s) from which blinks should be removed (cell array or
%          string, field(s) of ii_data)
%   pchan: channel name containing pupil data (string)
%   pval:  threshold, in pupil channel units, below which a blink is
%          identified
%   pri:   time before detected blink onset to exclude (ms)
%   fol:   time following detected blink offset to exclude (ms)
%
% Example (run from iEye toolbox home):
% 
% load('exdata1.mat');
% figure; 
% subplot(2,1,1);hold on;
% plot(ii_data.X,'k-','LineWidth',0.5); plot(ii_data.Y,'k-','LineWidth',0.5);
% [ii_data, ii_cfg] = ii_blinkcorrect(ii_data,ii_cfg,{'X','Y'},'Pupil',1500,150,50);  
% plot(ii_data.X,'-','LineWidth',2);plot(ii_data.Y,'-','LineWidth',2);
% xlabel('Samples');
% title('Before (thin black) and after (thick) blink correction');
% subplot(2,1,2);
% plot(ii_data.Pupil);
% xlabel('Samples');
% ylabel('Pupil size');

% Updated TCS 8/11/2017 - no base space variables
% Updated TCS 8/23/2017 - support for percentile-based thresholding
% Updated TCS 11/21/2018 - now save a PupilZ channel as well into ii_data
%
% TODO: implement a pupil velocity/acceleration version of this

% If not passed, get arguments


% TODO: check this...
if nargin == 2 % only data, cfg
    prompt = {'Channel to Correct', 'Pupil Channel', 'Pupil Threshold', '# Prior Samples', '# Following Samples'};
    dlg_title = 'Blink Correction';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines);
    
    chan = answer{1};
    pchan = answer{2};
    pval = str2num(answer{3});
    pri = str2num(answer{4});
    fol = str2num(answer{5});
end


if nargin < 3
    chan = {'X','Y'};
end

if nargin < 4
    pchan = 'Pupil';
end

if nargin < 5 || isempty(pval)
    pval = 1000;
end

if nargin < 6 || isempty(pri)
    pri = 100;
end

if nargin < 7 || isempty(fol)
    fol = 100;
end

if ismember(varargin,'prctile')
    pval_prctile = pval;
    pval = prctile(ii_data.(pchan)(ii_data.(pchan)~=0),pval);
end


% convert pri from milliseconds to samples given sampling rate
pri = round((pri/1000)*ii_cfg.hz);
fol = round((fol/1000)*ii_cfg.hz);


if ~iscell(chan)
    chan = {chan};
end

% first, find blinks
pupil = ii_data.(pchan);

blink_mask = pupil <= pval; 
blink_onsets  = find(diff(blink_mask)== 1);
blink_offsets = find(diff(blink_mask)==-1);

% handle case(s) where blink at beginning/end, so blink_onsets &
% blink_offsets are different lengths. in general, we'd hope that
% blink_onsets(1) is always < blink_offsets(1), etc. this will not be true
% if the first offset (upwards threshold crossing) precedes the first
% actual blink, or if the eye is closed at end of run (no offset)

if blink_onsets(1) > blink_offsets(1)
    
    blink_onsets = [1; blink_onsets];
    
end

if length(blink_offsets) == (length(blink_onsets)-1)
    blink_offsets = [blink_offsets; length(pupil)];
end

blink_window = [blink_onsets-pri blink_offsets+fol];
blink_window(blink_window(:,1)<1,1) = 1;
blink_window(blink_window(:,2)>length(pupil),2) = length(pupil);

blinkvec = zeros(size(pupil));

for tt = 1:size(blink_window,1)
    blinkvec(blink_window(tt,1):blink_window(tt,2)) = 1;
end

% then, replace those with nan, and save an appropriate variable in ii_cfg

for cc = 1:length(chan)
    ii_data.(chan{cc})(blinkvec==1) = nan;
end


ii_cfg.blinks = blink_window;
ii_cfg.blinkvec = blinkvec;

% save some information about how blink correction was performed
if ismember(varargin,'prctile')
    ii_cfg.blink.thresh_mode = 'prctile';
else
    ii_cfg.blink.thresh_mode = 'raw';
end
ii_cfg.blink.thresh = pval;
ii_cfg.blink.mask = blink_mask; % detected blink samples; prior to censoring samples before/after
ii_cfg.blink.window = [pri fol];


ii_cfg.history{end+1} = sprintf('ii_blinkcorrect - corrected channels %s with pupil channel %s, threshold of %d (%d/%d) on %s ', horzcat(chan{:}), pchan, pval, pri,fol, datestr(now,30));


% TODO: create a 'new channel' function that implements a few of the
% bookkeeping things required below...

% save a new channel into ii_data called sprintf('%sZ',pchan) that's just
% Zscore ii_data.Pupil(ii_data.Pupil~=0); nan where Pupil==0
if ismember(varargin,'prctile')
    pchanZ = sprintf('%sZ',pchan);
    ii_data.(pchanZ) = nan(size(ii_data.Pupil));
    ii_data.(pchanZ)(ii_data.(pchan)~=0) = zscore(ii_data.(pchan)(ii_data.(pchan)~=0));
    
    ii_cfg.lchan{1}{end+1} = pchanZ;
    
    % save the threshold as a Z-score too, in case we want to check it..    
    ii_cfg.blink.threshZ = prctile(ii_data.(pchanZ)(~isnan(ii_data.(pchanZ))),pval_prctile);

    % because this channel is used to look at pupil size and look for
    % lingering, unscored blinks, we also nan-out detected blinks
    ii_data.(pchanZ)(ii_cfg.blinkvec==1)=NaN;

end

return

