function [nchan,lchan,schan,cfg] = ii_openifg(ifg_file)
%Load configuration file
%   This function is used to load the contents of an iEye configuration
%   file (*.ifg).

if nargin ~= 1
    [filename, pathname] = uigetfile('*.ifg', 'Select *.ifg file');
    
    if isequal(filename,0)
        disp('User selected Cancel');
    else
        disp(['User selected', fullfile(pathname, filename)]);
        ifg_file = fullfile(pathname, filename);
    end
end

str = fopen(ifg_file);
cfg = ifg_file;

C = textscan(str,'%s');
nchan = C{1}{1};
lchan = C{1}{2};
schan = C{1}{3};
end

