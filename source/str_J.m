function [J_xx,J_xy,J_yy]=str_J(imx,imy)
% function [J_xx,J_yy,J_xy]=str_J(imx,imy)

J_xx = imx.*imx;
J_xy = imx.*imy;
J_yy = imy.*imy;
