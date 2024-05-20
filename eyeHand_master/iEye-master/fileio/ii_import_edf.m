function [ii_data,ii_cfg] = ii_import_edf(edf_file, ifg_file, data_file, varargin)
%IMPORT EYELINK EDF FILES
%   This function will import Eyelink EDF files but requires 'edf2asc'
%   command (from Eyelink) be installed in MATLAB's path. A config (*.ifg)
%   file is also required. After import, you can save as binary MAT file
%   for faster file I/O. Config information will be baked into ii_cfg
%   structure. Saved data file (data_file) will contain ii_data and ii_cfg
%   structures, for use by other iEye functions. 
%
% ii_import_edf prompts user for EDF file, IFG file
%
% ii_import_edf(edf_file, ifg_file) saves to edf_file_iEye.mat
%
% ii_import_edf(edf_file, ifg_file, []) prompts user for name/location of
% output file
%
% ii_import_edf(edf_file, ifg_file, data_file) reads
% edf_file with configuration ifg_file and saves data into data_file. if
% any arguments are empty, they are prompted via GUI.
%
% [ii_data, ii_cfg] = ii_import_edf(edf_file, ifg_file, data_file) returns
% ii_data, ii_cfg structures for use by other iEye programs
%
% ii_import_edf(edf_file, ifg_file, data_file, 'oldstyle') saves each
% channel separately, rather than packing into ii_data (for backwards
% compatibility)
%
% Misc info:
% MONOCULAR ONLY AT THE MOMENT
% Ignores corneal reflection flag information at the moment
% Automatically pulls out x,y,pupil sample data
% Will always be saved samples, messages
% i.e. x,y,pupil, message1, message2, message3
% Message variables are CASE SENSITIVE
%
% *IMPORTANT* MAKE SURE EDF2ASC COMMAND IS IN ENV PATH
% Use ii_init.m - but check the pathname within that function for your
% local install

% Updated TCS 8/17/2017 & prior - saves out ii_data, ii_cfg structs instead
% of channel variables. Can use ii_unpackdata to put all channels into base
% workspace (not advised!)

if nargin < 1 || isempty(edf_file)
    [filename, pathname] = uigetfile('*.edf', 'Select EDF file');
    edf_file = fullfile(pathname, filename);
    if filename==0
        fprintf('User selected cancel\n');
        return
    end 
end

if nargin < 2 || isempty(ifg_file)
    [filename_ifg, pathname] = uigetfile('*.ifg', 'Select IFG file');
    ifg_file = fullfile(pathname, filename_ifg);
    
    if filename_ifg==0
        fprintf('User selected cancel\n');
        return
    end
end

if nargin < 3
   data_file = sprintf('%s_iEye.mat',edf_file(1:(end-4)));
end

if isempty(data_file)
    iEye_file = sprintf('%s_iEye.mat',edf_file(1:(end-4)));
    [filename_data, pathname] = uiputfile(iEye_file, 'Create data file');
    data_file = fullfile(pathname, filename_data);    
    if filename_data==0
        fprintf('User selected cancel\n');
        return
    end    
end

% GET CONFIG
[nchan,lchan,schan,cfg] = ii_openifg(ifg_file);
nchan = str2num(nchan);
schan = str2num(schan);
vis = lchan; % mrugank: what's the point of this? 11/01
lchan = textscan(lchan,'%s','delimiter',',');

echan = nchan - 3;

% EXTRACT SAMPLES
% TCS: edf2asc =64 doesn't seem to work on 64bit macs, but old version is
% fine
%if strcmpi(computer(),'MACI64')
%    [status,result] = system(['edf2asc_x86_64 -t -c -s -miss 0 ' edf_file]);
%else
[status,result] = system(['edf2asc -t -c -v -y -s -miss 0 ' edf_file]);
%end
disp(status);
disp(result);
asc_samp_file = strrep(edf_file, '.edf', '.asc');

fid = fopen(asc_samp_file,'r');
%M = textscan(fid,'%f %f %f %f %*s');
sample_data = textscan(fid,'%f %f %f %f %*s %*s %*s %*s %*s');
delete (asc_samp_file);
sample_data = cell2mat(sample_data);

% Added from Nathan *(fixes temporal order of events)
sample_data = sortrows(sample_data, 1);
%s_num = sample_data(:,1);
%x = sample_data(:,2);
%y = sample_data(:,3);
%pupil = sample_data(:,4);


% EXTRACT EVENTS
[status,result] = system(['edf2asc -t -c -v -y -e -miss 0 ' edf_file]);
disp(status);
disp(result);
asc_evnt_file = strrep(edf_file, '.edf', '.asc');

fid = fopen(asc_evnt_file,'r');
event_data = textscan(fid,'%s', 'delimiter','\n');
delete(asc_evnt_file);
event_data = event_data{1};
mline = 1;
Mess = {};

% GET MSG EVENTS
token = strtok(event_data);
for v = 1:length(token)
    val_channame = strcmpi(token(v),'MSG');
    if val_channame == 1
        Mess(mline) = event_data(v);
        mline = mline + 1;
    end
end

Mess = Mess';
% Using strtok to separate text by first space. This produces samp_n which 
% has sample number, varbl which has value_holders and vval which has
% values
[~, remain] = strtok(Mess);
[samp_n, remain] = strtok(remain);
[varbl, vval] = strtok(remain);
samp_n = str2double(samp_n);
vval = str2double(vval);

% Added from Nathan *(fixes temporal order of events)
[samp_n, sort_I] = sort(samp_n);
varbl = varbl(sort_I);
vval = vval(sort_I);

% SEARCH MSG EVENTS FOR VARIABLES XDAT, TarX and TarY
for i = 4:nchan
    mline = 1;
    channame = lchan{1}{i};
    samp_vval = [];
    
    for v = 1:length(varbl)
        val_channame = strcmpi(varbl(v),channame); % TS: no longer case-sensitive!!!
        if val_channame == 1
            samp_vval(mline,:) = [samp_n(v) vval(v)];
            mline = mline + 1;
        end
    end
    
    % GET INDICES & SET VALUES
    li = 1;
    % ci = 1;
    cv = 0;
    sample_data(:,(i+1)) = 0;
    
     % find sample index from samp_vval 1st column in sample_data
    for h = 1:length(samp_vval)
        sample_idx = find(sample_data(:,1)==samp_vval(h,1));
        if ~isempty(sample_idx)
            sample_data((sample_idx:length(sample_data)),(i+1)) = samp_vval(h,2);
            li = sample_idx;
        else
            samp_vval(h,1) = samp_vval(h,1) - 1;
            sample_idx = find(sample_data(:,1)==samp_vval(h,1));
            sample_data((sample_idx:length(sample_data)),(i+1)) = samp_vval(h,2);
            li = sample_idx;
        end
    end
end

% CREATE FILE MATRIX

sample_data(:,1) = [];

% also insert channel data into eyedata, to be saved out into data_file
eyedata = [];
for i = 1:nchan
    channame = lchan{1}{i};
    cvalue = sample_data(:,i);
    eyedata.(lchan{1}{i})=sample_data(:,i);
end

x = sample_data(:,1);


% CREATE II_CFG STRUCT
dt = datestr(now, 30);%'mmmm dd, yyyy HH:MM:SS.FFF AM');

ii_cfg.cursel = [];
ii_cfg.sel = x*0;
ii_cfg.cfg = cfg;
ii_cfg.vis = vis;
ii_cfg.nchan = nchan;
ii_cfg.lchan = lchan;
ii_cfg.hz = schan;
ii_cfg.velocity = [];
ii_cfg.tcursel = [];
ii_cfg.tsel = x*0;
ii_cfg.tindex = 0;
ii_cfg.saccades = [];
ii_cfg.history{1} = ['EDF imported ', dt];
ii_cfg.edf_file = edf_file;
ii_cfg.microsacc =[]; 
ii_data = eyedata;

% SAVE FILE
% below focuses on only the necessary variables - sample_data is probably redundant,
% but some later code may use it so holding onto it for now. both save
% commands are necessary as it's impossible to save struct fields and
% variables simultaneously it seems. This *vastly* cuts down on amount of
% data saved (58-->3 MB for a single run)
% [only need -append for 'oldstyle' saving]

% only save if you give a filename, or you don't and don't get output args
if nargin >= 3 || (nargin<3 && nargout == 0)
    if ~isempty(varargin) && strcmpi(varargin{1},'oldstyle')
        save(data_file,'ii_cfg','edf_file','sample_data','ii_data');
        save(data_file,'-struct','eyedata','-append');
    else % default save state - ii_cfg and ii_data encapsulate all info
        save(data_file,'ii_cfg','ii_data','edf_file');
    end
end

% maybe will help us avoid weird errors?
fclose('all');

end


