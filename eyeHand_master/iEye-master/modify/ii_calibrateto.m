function ii_calibrateto(chan, chan2, polynum)
%Calibrate one channel to another
%   This function will calibrate the values of one channel to the values of
%   another. It takes the channel to calibrate as the first argument,
%   followed by the channel to calibrate to, and a polynomial degree.

ii_cfg = evalin('base', 'ii_cfg');
cursel = ii_cfg.cursel;
sel = ii_cfg.sel;

if length(cursel) < 1
    disp('No selections made for calibration');
else
    if nargin ~= 3
        prompt = {'Enter channel to calibrate:', 'Enter channel to calibrate to:', 'Polynomial Degree'};
        dlg_title = 'Calibrate To';
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines);
        
        chan = answer{1};
        chan2 = answer{2};
        polynum = str2num(answer{3});
    end
    
    basevars = evalin('base','who');
    
    if ismember(chan,basevars)
        if ismember(chan2,basevars)
            c1 = evalin('base',chan);
            c2 = evalin('base',chan2);
            
            c1s = c1(sel==1);
            c2s = c2(sel==1);
            
            z = polyfit(c1s,c2s,polynum);
            fit = polyval(z,c1);
            
            disp('Calibration saved');
            assignin('base',chan,fit);
            ii_replot;
            
%             % Plot fit at selections
%             figure;            
%             scatter(fit(sel==1),c2(sel==1));
%             axis equal
%             
%             % Plot overall fit
%             figure;            
%             scatter(fit,c2);
%             axis equal

        dt = datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM');
        ii_cfg = evalin('base','ii_cfg');
        ii_cfg.history{end+1,1} = sprintf('Calibrated %s to %s with polyum %d on %s ', chan, chan2,polynum, dt);
        putvar(ii_cfg);
            
        else
            disp('Channel to calibrate to does not exist in workspace');
        end
    else
        disp('Channel to calibrate does not exist in worksapce');
    end
end
end

