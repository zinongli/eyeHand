% copy column contents:
% 1,2: target location t_0 in wac pixels 
% 3: p_eye scalar (rho or theta)
% 4: p_hand scalar (rho or theta)
% 5,6: start location fixation in wac pixels
% 7,8: f_eye location in wac pixels
% 9: absolute onset time
% 10,11: e_hand location in wac pixels
% 12,13: f_hand location in wac pixels


% 13,14: endpoint x and y in mm
% 15: target size in mm
% 16: actual duration of the reach
% 17: error size in mm
% 18: switch logical array
% 19,20: start position in mm
% 21: actual reach distance
% 22: average speed
% 23: error along the reach direction (vector projection) in mm
% 24,25: relative target position in mm

%%
% Define the file name or path
filename = '0425Sacc1st.csv';

% Read the CSV file into a matrix
eyeData = readmatrix(filename);

set(0, 'DefaultAxesFontSize', 16); % Set global default font size for axes

index = NaN(length(data),1);
for i = 1:length(data)
    index(i) = ~isnan(sum(data(i,10:13)));
    data(i,5:6) = [512,384];
end

valid = data(index==true,:);
validTraXH = traXtotalH(index==true,:);
validTraYH = traYtotalH(index==true,:);
validTraXE = traXtotalE(index==true,:);
validTraYE = traYtotalE(index==true,:);

Affine2d =tform.T(1:2,1:2);
[~,s,~] = svd(Affine2d);
proj2tablet = 1./mean([s(1,1),s(2,2)]);
pixellength = 0.248;
copy = valid;

copyWac(:,1:2) = transformPointsInverse(tform,valid(:,1:2));
copyWac(:,5:6) = transformPointsInverse(tform,valid(:,5:6));
copyWac(:,7:8) = transformPointsInverse(tform,valid(:,7:8));
copyWac(:,10:11) = transformPointsInverse(tform,valid(:,10:11));
copyWac(:,12:13) = transformPointsInverse(tform,valid(:,12:13));

copyProj(:,1:2) = valid(:,1:2);
copyProj(:,5:6) = valid(:,5:6);
copyProj(:,7:8) = valid(:,7:8);
copyProj(:,10:11) = valid(:,10:11);
copyProj(:,12:13) = valid(:,12:13);

%%
[tarThetas, tarRhos] = cart2pol(copyProj(:,1) - copyProj(:,5), copyProj(:,2) - copyProj(:,6));
[eHandThetas, eHandRhos] = cart2pol(copyProj(:,10) - copyProj(:,5), copyProj(:,11) - copyProj(:,6));
[fEyeThetas, fEyeRhos] = cart2pol(copyProj(:,7) - copyProj(:,5), copyProj(:,8) - copyProj(:,6));
[fHandThetas, fHandRhos] = cart2pol(copyProj(:,12) - copyProj(:,5), copyProj(:,13) - copyProj(:,6));
[eEyeThetas, eEyeRhos] = cart2pol(eyeData(:,5) - eyeData(:,3), eyeData(:,6) - eyeData(:,4));

%%
eHandRhosRecent = reshape(eHandRhos - tarRhos,72,5);
fHandRhosRecent = reshape(fHandRhos - tarRhos,72,5);

plot(1:72,mean(eHandRhosRecent,2),'-o')
hold on
plot(1:72,mean(fHandRhosRecent,2),'-o')
hold off

%%
eHandRhosRecent = eHandRhos - tarRhos;
% fHandRhosRecent = fHandRhos - tarRhos;

plot(1:360,eHandRhosRecent,'-o')
hold on
plot(1:360,copy(:,4),'-o')
yline(0,'--')
hold off
xlabel('Trial #')
ylabel('Rho Error from Initial Target T_0 (rad)')
legend('Endpoint','Perturbation')
title('Reach Error vs Reach Perturbation')
%%
figure(1)
subplot(2,2,1)
eHandThetasRecent = eHandThetas - tarThetas;

plot(1:360,rad2deg(eHandThetasRecent),'-o')
hold on
plot(1:360,rad2deg(copy(:,4)),'-o')
yline(0,'--')
hold off
ylim([-12,12])
xlabel('Trial #')
ylabel('Rotational Error (deg)')
legend('Endpoints','Perturbation')
title('Reach Error vs Reach Perturbation')

subplot(2,2,3)
eHandRhosRecent = eHandRhos - tarRhos;

plot(1:360,eHandRhosRecent .* proj2tablet .* pixellength,'-o')
hold on
plot(1:360,copy(:,3) .* proj2tablet .* pixellength,'-o')
yline(0,'--')
hold off
ylim([-28,28])
xlabel('Trial #')
ylabel('Gain Error (mm)')
% legend('Endpoints','Perturbation')
title('Reach Error vs Saccade Perturbation')

subplot(2,2,2)
% Assuming the data is already loaded in a variable named 'data'
data = eHandThetasRecent .* proj2tablet .* pixellength;  % Example data, replace this with your actual data

% Parameters
N = length(data);    % Number of data points
Fs = 360;           % Sampling rate in Hz (replace with your actual sampling rate)

% FFT computation
fft_data = fft(data);

% Amplitude spectrum calculation
amplitudes = abs(fft_data(1:N/2+1));

% Frequency axis preparation
frequencies = (0:N/2) * (Fs / N);

% Plotting the amplitude spectrum
plot(frequencies, amplitudes,'-');
hold on
xline(5,'--r')
hold off
title('Rotational Adaptation Amplitude');
xlabel('Cycle per session');
ylabel('Amplitude');
xlim([-20,180])
legend('','Perturbation')
% Enhance plot visibility
grid on;

subplot(2,2,4)
% Assuming the data is already loaded in a variable named 'data'
data = eHandRhosRecent .* proj2tablet .* pixellength;  % Example data, replace this with your actual data

% Parameters
N = length(data);    % Number of data points
Fs = 360;           % Sampling rate in Hz (replace with your actual sampling rate)

% FFT computation
fft_data = fft(data);

% Amplitude spectrum calculation
amplitudes = abs(fft_data(1:N/2+1));

% Frequency axis preparation
frequencies = (0:N/2) * (Fs / N);

% Plotting the amplitude spectrum
plot(frequencies, amplitudes,'-');
hold on
xline(3,'--r')
hold off
title('Gain Adaptation Amplitude');
xlabel('Cycle per session');
ylabel('Amplitude');
xlim([-20,180])

% Enhance plot visibility
grid on;
%%
figure(2)
subplot(2,2,1)
eEyeThetasRecent = eEyeThetas - tarThetas;

plot(1:360,rad2deg(eEyeThetasRecent),'-o')
hold on
plot(1:360,rad2deg(copy(:,4)),'-o')
yline(0,'--')
hold off
ylim([-12,12])
xlabel('Trial #')
ylabel('Rotational Error (deg)')
legend('Endpoints','Perturbation')
title('Saccade Error vs Reach Perturbation')

subplot(2,2,3)
eEyeRhosRecent = eEyeRhos - tarRhos;

plot(1:360,eEyeRhosRecent .* proj2tablet .* pixellength,'-o')
hold on
plot(1:360,copy(:,3) .* proj2tablet .* pixellength,'-o')
yline(0,'--')
hold off
ylim([-28,28])
xlabel('Trial #')
ylabel('Gain Error (mm)')
% legend('Endpoints','Perturbation')
title('Saccade Error vs Saccade Perturbation')

subplot(2,2,2)
% Assuming the data is already loaded in a variable named 'data'
data = eEyeThetasRecent .* proj2tablet .* pixellength;  % Example data, replace this with your actual data

% Parameters
N = length(data);    % Number of data points
Fs = 360;           % Sampling rate in Hz (replace with your actual sampling rate)

% FFT computation
fft_data = fft(data);

% Amplitude spectrum calculation
amplitudes = abs(fft_data(1:N/2+1));

% Frequency axis preparation
frequencies = (0:N/2) * (Fs / N);

% Plotting the amplitude spectrum
plot(frequencies, amplitudes,'-');
hold on
xline(5,'--r')
hold off
title('Rotational Adaptation Amplitude');
xlabel('Cycle per session');
ylabel('Amplitude');
xlim([-20,180])
legend('','Perturbation')
% Enhance plot visibility
grid on;

subplot(2,2,4)
% Assuming the data is already loaded in a variable named 'data'
data = eEyeRhosRecent .* proj2tablet .* pixellength;  % Example data, replace this with your actual data

% Parameters
N = length(data);    % Number of data points
Fs = 360;           % Sampling rate in Hz (replace with your actual sampling rate)

% FFT computation
fft_data = fft(data);

% Amplitude spectrum calculation
amplitudes = abs(fft_data(1:N/2+1));

% Frequency axis preparation
frequencies = (0:N/2) * (Fs / N);

% Plotting the amplitude spectrum
plot(frequencies, amplitudes,'-');
hold on
xline(3,'--r')
hold off
title('Gain Adaptation Amplitude');
xlabel('Cycle per session');
ylabel('Amplitude');
xlim([-20,180])

% Enhance plot visibility
grid on;
%%
figure('Position', [100, 100, 1400, 300]); % Create a figure window 600x400 pixels, positioned 100 pixels from the left and bottom of the screen
subplot(1,2,1)
eEyeThetasRecent = eHandThetas - tarThetas;

plot(1:360,rad2deg(eEyeThetasRecent),'-o')
hold on
plot(1:360,rad2deg(copy(:,4)),'-o')
yline(0,'--')
hold off
ylim([-12,12])
xlim([0,380])
xlabel('Trial #')
ylabel('Rotational Error (deg)')
legend('Endpoints','Perturbation')
title('Reach Error vs Reach Perturbation')

subplot(1,2,2)
% Assuming the data is already loaded in a variable named 'data'
data = eEyeThetasRecent .* proj2tablet .* pixellength;  % Example data, replace this with your actual data

% Parameters
N = length(data);    % Number of data points
Fs = 360;           % Sampling rate in Hz (replace with your actual sampling rate)

% FFT computation
fft_data = fft(data);

% Amplitude spectrum calculation
amplitudes = abs(fft_data(1:N/2+1));

% Frequency axis preparation
frequencies = (0:N/2) * (Fs / N);

% Plotting the amplitude spectrum
plot(frequencies, amplitudes,'-');
hold on
xline(5,'--r')
hold off
title('Rotational Adaptation Frequency Amplitude');
xlabel('Cycle per session');
ylabel('Amplitude');
xlim([-20,180])
legend('','Perturbation')
% Enhance plot visibility
grid on;

%%
figure('Position', [100, 100, 1400, 400]); % Create a figure window 600x400 pixels, positioned 100 pixels from the left and bottom of the screen

subplot(1,2,1)
eEyeRhosRecent = eHandRhos - tarRhos;

plot(1:360,eEyeRhosRecent .* proj2tablet .* pixellength,'-o')
hold on
plot(1:360,copy(:,3) .* proj2tablet .* pixellength,'-o')
yline(0,'--')
hold off
ylim([-28,28])
xlabel('Trial #')
ylabel('Gain Error (mm)')
legend('Endpoints','Perturbation')
title('Reach Error vs Saccade Perturbation')


subplot(1,2,2)
% Assuming the data is already loaded in a variable named 'data'
data = eEyeRhosRecent .* proj2tablet .* pixellength;  % Example data, replace this with your actual data

% Parameters
N = length(data);    % Number of data points
Fs = 360;           % Sampling rate in Hz (replace with your actual sampling rate)

% FFT computation
fft_data = fft(data);

% Amplitude spectrum calculation
amplitudes = abs(fft_data(1:N/2+1));

% Frequency axis preparation
frequencies = (0:N/2) * (Fs / N);

% Plotting the amplitude spectrum
plot(frequencies, amplitudes,'-');
hold on
xline(3,'--r')
hold off
title('Gain Adaptation Frequency Amplitude');
xlabel('Cycle per session');
ylabel('Amplitude');
xlim([-20,180])
legend('','Perturbation')
% Enhance plot visibility
grid on;
%%
eEyeThetasRecent = reshape(eHandThetas - tarThetas,120,3);

plot(1:120,mean(eEyeThetasRecent,2),'-o')
hold on
plot(1:120,copy(1:120,3),'-o')
hold off
%%
for i = 1:360
plot(validTraXE(i,:),validTraYE(i,:),'-o')
hold on
% plot(copyPro(i,5),copyPro(i,6),'*k','MarkerSize',15)
plot(copyProj(i,1),copyProj(i,2),'*k','MarkerSize',15)
plot(copyProj(i,7),copyProj(i,8),'*r','MarkerSize',15)
plot(copyProj(i,12),copyProj(i,13),'*b','MarkerSize',15)
hold off
pause(2)
end
%%
% copy(:,10) = sqrt(sum((copy(:,1:2) - copy(:,8:9)).^2,2)) .* pixellength;
% copy(:,[11,12]) = [copy(:,1)*pixellength (1080 - copy(:,2))*pixellength];
% copy(:,[13,14]) = [copy(:,6)*pixellength (1080 - copy(:,7))*pixellength]; % 1080 = tablet pixel height
% copy(:,15) = valid(:,10) .* pixellength .* proj2tablet ./ 2; % proj2tablet = projetor size to tablet size (physical size), /2 is diameter vs radius
% copy(:,16) = copy(:,5) - copy(:,4);
% copy(:,17) = sqrt( (copy(:,13)-copy(:,11)).^2 + (copy(:,14)-copy(:,12)).^2 );
% copy(:,27) = 1:length(copy);
% copy(:,19:20) = (copy(:,6:7) - copy(:,8:9)) .* pixellength;
% copy(:,21) = sqrt(sum((copy(:,6:7) - copy(:,8:9)).^2,2)) .* pixellength;
% copy(:,22) = copy(:,21) ./ copy(:,16);
% copy(:,24:25) = (copy(:,1:2) - copy(:,8:9)) .* pixellength;% relative target coordinate
% copy(:,23) = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)) - 1) .*copy(:,10);
% copy(:,26) = valid(:,11);
% copy(:,28) = copy(:,26) ~= 0;
% endPoints = (copy(:,6:7) - copy(:,8:9)) .* pixellength;% relative endpoint coordinate
% projScale = (abs(dot(copy(:,19:20),copy(:,24:25),2) ./ dot(copy(:,24:25),copy(:,24:25),2)));
% rejections = endPoints - projScale.* copy(:,24:25);
% tooRight = NaN(1,length(copy));
% for i = 1:length(copy)
%     tooRight(i) = sum(cross([endPoints(i,:),0],[rejections(i,:),0])<0);
% end
% tooRight(tooRight==0) = -1;
% rejLength = sqrt(rejections(:,1).^2 + rejections(:,2).^2);
% copy(:,29) = rejLength .* tooRight';
% 
% % sigmoid fitting max speed from trajs
% reCenteredTrajX = NaN(360,270);
% reCenteredTrajY = NaN(360,270);
