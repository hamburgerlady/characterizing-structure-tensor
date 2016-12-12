function [J_xx_sm,J_xy_sm,J_yy_sm,vf11,vf12]=str_Jsm(J_xx,J_xy,J_yy,sigma)
% function [J_xx_sm,J_xx_sm,J_xx_sm]=str_Jsm(J_xx,J_xy,J_yy,sigma)


fs = ceil(6*sigma);
[xx,yy]=meshgrid(-fs:fs);

gk=exp(-(xx.^2+yy.^2)/2/sigma^2)/(2*pi*sigma^2);
gk = gk/sum(gk(:));
vf1 = sum(gk(:).^2);
vf2 = sum(sum(gk(:,1:end-2).*gk(:,3:end)))*2/4; % cov = 1/4 var
vf11 = vf1+vf2;
vf12 = vf1;
J_xx_sm = conv2(J_xx,gk,'same');
J_xy_sm = conv2(J_xy,gk,'same');
J_yy_sm = conv2(J_yy,gk,'same');
