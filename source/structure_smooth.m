function imout = structure_smooth(im,s1)

im = double(im);
ksfactor=1.5;
sigmamin=0.4;
sigmamax=5;


sigma2 = 1.7;
xxyysmcovf = 1/1.6;

[imx,imy,vf]=str_diff_nosmooth(im);
vdx = vf*s1^2;



[J_xx,J_xy,J_yy]=str_J(imx,imy);
[J_xx_sm,J_xy_sm,J_yy_sm,vfJ11,vfJ12]=str_Jsm(J_xx,J_xy,J_yy,sigma2);

j11 = J_xx_sm;
j12 = J_xy_sm;
j22 = J_yy_sm;

jdisc = abs(j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;
lille = 1e-9;
lam2 = j11/2 + j22/2 + jdisc;
lam1 = j11/2 + j22/2 - jdisc;

E1_1 = (j11/2 + j22/2 - jdisc)./(lille+j12) - j22./(lille+j12);
E1_2 = 1;

sum1 = sqrt(E1_1.^2+E1_2.^2);

E1_1 = E1_1./sum1;
E1_2 = E1_2./sum1;

A = -E1_1;
B = -E1_2;


lam1e_m = 0.5*(vdx+vdx-sqrt(2/pi)*sqrt(2*vdx^2*vfJ12)-sqrt(2*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11)/pi));
vcovfac1 = 1/4;
lam1e_v = vcovfac1*(2*vdx^2*vfJ11+2*vdx^2*vfJ11+(pi-2)/pi*2*vdx^2*vfJ12+(pi-2)/pi*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11));



lambdamin = lam1e_m;
dfac = sqrt(lam1e_v)/lambdamin;
d = dfac/lambdamin;

%keyboard
imout = structure_smooth_sum(im,lam1,lam2,A,B,sigmamin,sigmamax,d,lambdamin,ksfactor);




