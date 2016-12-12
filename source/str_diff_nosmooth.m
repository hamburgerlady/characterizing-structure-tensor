function [imx,imy,vf]=str_diff_nosmooth(im)
% function [imx,imy]=str_diff_nosmooth(im)


gkx=[1 0 -1];
gky=gkx';



vf = sum(gkx(:).^2);

imx = conv2(im,gkx,'same');
imy = conv2(im,gky,'same');