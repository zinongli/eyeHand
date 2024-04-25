function drawArrow2Center(window,cx,cy,mx,my,r,color)
angle = atan2(my-cy,mx-cx);
arrowLength = 30;
arrowWingLength = 10;
arrowWidth = 2;

arrowTipX = cx + r .* cos(angle);
arrowTipY = cy + r .* sin(angle);

arrowBaseX = cx + (r+arrowLength) .* cos(angle);
arrowBaseY = cy + (r+arrowLength) .* sin(angle);

arrowBackLeftX = arrowTipX + arrowWingLength * cos(angle - pi./4);
arrowBackLeftY = arrowTipY + arrowWingLength * sin(angle - pi./4);

arrowBackRightX = arrowTipX + arrowWingLength * cos(angle + pi./4);
arrowBackRightY = arrowTipY + arrowWingLength * sin(angle + pi./4);

Screen('DrawLine', window, color,arrowBaseX,arrowBaseY,arrowTipX,arrowTipY,arrowWidth);
Screen('DrawLine', window, color,arrowTipX,arrowTipY,arrowBackLeftX,arrowBackLeftY,arrowWidth);
Screen('DrawLine', window, color,arrowTipX,arrowTipY,arrowBackRightX,arrowBackRightY,arrowWidth);

end