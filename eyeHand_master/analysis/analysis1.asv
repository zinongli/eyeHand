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

copyPro(:,1:2) = valid(:,1:2);
copyPro(:,5:6) = valid(:,5:6);
copyPro(:,7:8) = valid(:,7:8);
copyPro(:,10:11) = valid(:,10:11);
copyPro(:,12:13) = valid(:,12:13);

%%
[tarRhos, tarThetas] = cart2pol(copyPro(:,1) - copyPro(:,5), copyPro(:,2) - copyPro(:,6));
[eHandRhos, eHandThetas] = cart2pol(copyPro(:,10) - copyPro(:,5), copyPro(:,11) - copyPro(:,6));
[fEyeRhos, fEyeThetas] = cart2pol(copyPro(:,7) - copyPro(:,5), copyPro(:,8) - copyPro(:,6));
[fHandRhos, fHandThetas] = cart2pol(copyPro(:,12) - copyPro(:,5), copyPro(:,13) - copyPro(:,6));

%%
eHandRhosRecent = reshape(eHandRhos - tarRhos,72,5);
fHandRhosRecent = reshape(fHandRhos - tarRhos,72,5);

plot(1:72,mean(eHandRhosRecent,2),'-o')
hold on
plot(1:72,mean(fHandRhosRecent,2),'-o')
hold off

%%
eHandRhosRecent = eHandRhos - tarRhos;
fHandRhosRecent = fHandRhos - fEyeRhos;

plot(1:360,eHandRhosRecent,'-o')
hold on
plot(1:360,fHandRhosRecent,'-o')
hold off

%%
for i = 1:360
plot(validTraXE(i,:),validTraYE(i,:),'-o')
hold on
% plot(copyPro(i,5),copyPro(i,6),'*k','MarkerSize',15)
plot(copyPro(i,1),copyPro(i,2),'*k','MarkerSize',15)
plot(copyPro(i,7),copyPro(i,8),'*r','MarkerSize',15)
plot(copyPro(i,12),copyPro(i,13),'*b','MarkerSize',15)
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
