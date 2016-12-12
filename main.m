addpath('source');


%% plotfigs

disp('Creating Fig. 1 ...');

N=100000;
k1 = rand*10;
k2 = rand*10;
th1 = rand*10;
th2 = rand*10;

m1 = k1*th1;
v1 = k1*th1^2;
m2 = k2*th2;
v2 = k2*th2^2;

mp = m1*m2;
vp = m1^2*v2+m2^2*v1+v1*v2;

kp = mp^2/vp;
thp = vp/mp;

ms = m1+m2;
vs = v1+v2;

ks = ms^2/vs;
ths = vs/ms;


X1 = gamrnd(k1,th1,N,1);
X2 = gamrnd(k2,th2,N,1);


[hh1,bb1]=hist(X1,N/500);
hh1 = hh1/sum(hh1)/mean(abs(diff(bb1)));
hh10 = gampdf(bb1,k1,th1);

figure(1);
clf
subplot(2,2,1);
bar(bb1,hh1,'c');
hold on
ll=plot(bb1,hh10,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})

[hh2,bb2]=hist(X2,N/500);
hh2 = hh2/sum(hh2)/mean(abs(diff(bb2)));
hh20 = gampdf(bb2,k2,th2);

subplot(2,2,2);
bar(bb2,hh2,'c');
hold on
ll=plot(bb2,hh20,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})

X = X1+X2;
[hhs,bbs]=hist(X,N/500);
hhs = hhs/sum(hhs)/mean(abs(diff(bbs)));
hhs0 = gampdf(bbs,ks,ths);

subplot(2,2,3);
bar(bbs,hhs,'c');
hold on
ll=plot(bbs,hhs0,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})


X = X1.*X2;
[hhp,bbp]=hist(X,N/500);
hhp = hhp/sum(hhp)/mean(abs(diff(bbp)));
hhp0 = gampdf(bbp,kp,thp);

subplot(2,2,4);
bar(bbp,hhp,'c');
hold on
ll=plot(bbp,hhp0,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})

set(gcf,'NextPlot','add');
axes;
h = title('Fig. 1. Top row: The distributions of two input gamma ... ');
set(gca,'Visible','off');
set(h,'Visible','on');




%%

disp('Creating Fig. 2. running 5000 iterations...');
nn = 5000;
boffa = zeros(6,nn);
boffa2 = zeros(6,nn);

fprintf('Iteration ');
for iii = 1:nn,
if mod(iii,100)==0,
    fprintf('%d ',iii);
end

N=100000;
k1 = rand*10;
k2 = rand*10;
th1 = rand*10;
th2 = rand*10;

m1 = k1*th1;
v1 = k1*th1^2;
m2 = k2*th2;
v2 = k2*th2^2;

mp = m1*m2;
vp = m1^2*v2+m2^2*v1+v1*v2;

kp = mp^2/vp;
thp = vp/mp;

ms = m1+m2;
vs = v1+v2;
ks = ms^2/vs;
ths = vs/ms;

X1 = gamrnd(k1,th1,N,1);
X2 = gamrnd(k2,th2,N,1);

X = X1+X2;
[hhs,bbs]=hist(X,N/500);
hhs = hhs/sum(hhs)/mean(abs(diff(bbs)));
hhs0 = gampdf(bbs,ks,ths);

X = X1.*X2;
[hhp,bbp]=hist(X,N/500);
hhp = hhp/sum(hhp)/mean(abs(diff(bbp)));
hhp0 = gampdf(bbp,kp,thp);
boffa(:,iii)=[norm(hhs-hhs0) norm(hhp-hhp0) k1 k2 th1 th2]';
boffa2(:,iii)=[KLDiv(hhs,hhs0) KLDiv(hhp,hhp0) k1 k2 th1 th2]';

end
fprintf('\n');

[h_err_sum,b_err_sum] = hist(log10(boffa(1,:)),80);
[h_err_prod,b_err_prod] = hist(log10(boffa(2,:)),80);
h_err_sum=h_err_sum/sum(h_err_sum)*100;
h_err_prod=h_err_prod/sum(h_err_prod)*100;

figure(2);
subplot(1,2,1);
bar(b_err_sum,h_err_sum,'c');
xlabel('Fig. 2. Histogram of the logarithm of the Kullback-Leibler divergence between data histograms and...');
subplot(1,2,2);
bar(b_err_prod,h_err_prod,'c');

%%
disp('Creating Fig. 3 ...');

N = 100000;
x = 0:0.2:50;
k = 2.5;
theta1 = 2;
theta2 = 3;

p1 = gampdf(x,k,theta1);
p2 = gampdf(x,k,theta2);

X1 = gamrnd(k,theta1,N,1);
X2 = gamrnd(k,theta2,N,1);

X = X1+X2;

%[hh,bb]=hist(X,x);
%hh = hh/sum(hh)/mean(abs(diff(bb)));

ks=k*(theta1+theta2)^2/(theta1^2+theta2^2);
ths=(theta1^2+theta2^2)/(theta1+theta2);
ps = gampdf(x,ks,ths);

figure(3);
clf
hold on


b = (theta2-theta1)/(theta1*theta2);
N1 = (b*x+2*k).*whittakerM(0, k+1/2, b*x);
D = (2*k+1)*gamma(k+1/2);
N2 = 2*(k+1).*whittakerM(1, k+1/2, b*x);
CC = x.^k*b^(-k+1).*exp(-theta1*x/(theta1*theta2))/(gamma(k)*theta1^k*theta2^k);
CC2 = b^(k-2)*sqrt(pi)*exp(-(1/2)*b*x).*(b*x).^(-k).*x.^(k-2)*2^(-2*k);
pt = CC.*(N1+N2)/D.*CC2;
pt(1)=0;

ll=plot(x,pt,'r');
set(ll,'LineWidth',2);
ll2 = plot(x,ps,'b--');
set(ll2,'LineWidth',2);

legend({'True distribution PDF','Approximate Gamma distribution PDF'})
xlabel('Fig. 3. The figure shows the model distribution and the exact theoretical distribution ...');

%%
disp('Creating Fig. 4 ...');

N = 1000;
s1 = 2;
sigma2 = 1.7;
%xxyysmcovf = 1/1.6;
im = randn(N)*s1;

[imx,imy,vf]=str_diff_nosmooth(im);
vdx = vf*s1^2;


% disp('imx:');
% disp([mean(imx(:)) 0])
% disp([var(imx(:)) vdx])

[J_xx,J_xy,J_yy]=str_J(imx,imy);

figure(4)
clf
subplot(1,4,1);
hold off
[hh,bb]=hist(J_xx(:),N);
hh = hh/sum(hh)/mean(abs(diff(bb)));
bar(bb,hh,'c');
hold on
hh2 = gampdf(bb,1/2,2*vdx);
ll=plot(bb,hh2,'r');
set(ll,'LineWidth',2);
% disp('J11:');
% disp([mean(J_xx(:)) vdx])
% disp([var(J_xx(:)) 2*vdx^2])
mb = max(bb)*0.2;
mh = max(hh)*1.1;
legend({'Histogram','Estimated PDF'})
axis([0 mb 0 mh])
title('J_{11}')

subplot(1,4,2);
hold off
[hh,bb]=hist(J_xy(:),N/5);
hh = hh/sum(hh)/mean(abs(diff(bb)));
bar(bb,hh,'c');
hold on
hh2 = besselk(0,abs(bb)/vdx)/pi/vdx;
ll=plot(bb,hh2,'r');
set(ll,'LineWidth',2);
% disp('J12:');
% disp([mean(J_xy(:)) 0])
% disp([var(J_xy(:)) vdx^2])
mb=0.7*min([max(bb) -min(bb)]);
mh = max(hh)*1.1;
legend({'Histogram','Estimated PDF'})
axis([-mb mb 0 mh])
title('J_{12}')

%[J_xx_sm,J_xy_sm,J_yy_sm,vfJ11,vfJ12]=str_Jsm(J_xx,J_xy,J_yy,sigma2);
[J_xx_sm,J_xy_sm,~,vfJ11,vfJ12]=str_Jsm(J_xx,J_xy,J_yy,sigma2);


subplot(1,4,3);
hold off
[hh,bb]=hist(J_xx_sm(:),N/5);
hh = hh/sum(hh)/mean(abs(diff(bb)));
bar(bb,hh,'c');
hold on
hh2 = gampdf(bb,1/2/vfJ11,2*vdx*vfJ11);
mh = max(hh2)*1.1;
ll=plot(bb,hh2,'r');
set(ll,'LineWidth',2);
mb = max(bb);
legend({'Histogram','Estimated PDF'})
axis([0 mb 0 mh])
title('J_{11} smoothed')

% disp('J_11_sm:');
% disp([mean(J_xx_sm(:)) vdx])
% disp([var(J_xx_sm(:)) 2*vdx^2*vfJ11])


subplot(1,4,4);
hold off
[hh,bb]=hist(J_xy_sm(:),N/5);

hh = hh/sum(hh)/mean(abs(diff(bb)));
bar(bb,hh,'c');
hold on
hhn = normpdf(bb,0,sqrt(2*vdx^2*vfJ12/1.1)); % varför dela med 1.1??
ll=plot(bb,hhn,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})
title('J_{12} smoothed')


% disp('J_12_sm:');
% disp([mean(J_xy_sm(:)) 0])
% disp([var(J_xy_sm(:)) 2*vdx^2*vfJ12])


%j11 = J_xx_sm;
%j12 = J_xy_sm;
%j22 = J_yy_sm;
%lam1 = j11/2 + j22/2 + (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;
%lam2 = j11/2 + j22/2 - (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;

%lam1e = j11/2 + j22/2 + abs(j12)/2+abs(j11-j22)/2;
%lam2e = j11/2 + j22/2 - abs(j12)/2-abs(j11-j22)/2;



%lam1e_m = 0.5*(vdx+vdx+sqrt(2/pi)*sqrt(2*vdx^2*vfJ12)+sqrt(2*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11)/pi));
%vcovfac1 = 1/1.9;
%lam1e_v = vcovfac1*(2*vdx^2*vfJ11+2*vdx^2*vfJ11+(pi-2)/pi*2*vdx^2*vfJ12+(pi-2)/pi*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11));

%lam2e_m = 0.5*(vdx+vdx-sqrt(2/pi)*sqrt(2*vdx^2*vfJ12)-sqrt(2*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11)/pi));
%vcovfac2 = 1/4;
%lam2e_v = vcovfac2*(2*vdx^2*vfJ11+2*vdx^2*vfJ11+(pi-2)/pi*2*vdx^2*vfJ12+(pi-2)/pi*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11));



% disp('Eigenvalues...')
% 
% disp([mean(lam1(:)) mean(lam1e(:)) lam1e_m]);
% disp([var(lam1(:)) var(lam1e(:)) lam1e_v]);
% 
% disp([mean(lam2(:)) mean(lam2e(:)) lam2e_m]);
% disp([var(lam2(:)) var(lam2e(:)) lam2e_v]);

subplot(1,4,2);
xlabel('Fig. 4. The figure shows histograms of data  and model distributions. From left to right: the unsmoothed j_{11} and j_{12}, the smoothed j_{11} and j_{12}.');


%%
disp('Creating Fig. 5 ...');
im0 = double(imread('cameraman.tif'))/10;
s1 = 4;
im = im0 + randn(size(im0))*s1;

figure(5);
subplot(3,2,1);
imagesc(im0)
colormap(gray(256))
axis off
axis equal
subplot(3,2,2);
imagesc(im);
colormap(gray(256))
axis off
axis equal

nn = 200;
sigma2 = 1.7;
xxyysmcovf = 1/1.6;
[imx,imy,vf]=str_diff_nosmooth(im);
[im0x,im0y,vf0]=str_diff_nosmooth(im0);
grady = sqrt(imx.^2+imy.^2);
gradybunch = 4;

vdx = vf*s1^2;
[J_xx,J_xy,J_yy]=str_J(imx,imy);
[J0_xx,J0_xy,J0_yy]=str_J(im0x,im0y);
[J_xx_sm,J_xy_sm,J_yy_sm,vfJ11,vfJ12]=str_Jsm(J_xx,J_xy,J_yy,sigma2);
[J0_xx_sm,J0_xy_sm,J0_yy_sm,vfJ110,vfJ120]=str_Jsm(J0_xx,J0_xy,J0_yy,sigma2);


subplot(3,2,3);
hold off
%[hhxx,bbxx]=hist(J_xx_sm(:)-J0_xx_sm(:),nn);
[hhxx,bbxx]=hist(J_xx_sm(grady<gradybunch),nn);

hhxx = hhxx/sum(hhxx)/mean(abs(diff(bbxx)));
bar(bbxx,hhxx,'c');
hold on
hhxxe = gampdf(bbxx,1/2/vfJ11,2*vdx*vfJ11);
mh = max(hhxxe)*1.1;
ll=plot(bbxx,hhxxe,'r');
set(ll,'LineWidth',2);
%mb = max(bbxx)/2;
legend({'Histogram','Estimated PDF'})
%axis([0 mb 0 mh])
axis([0 100 0 mh])
title('J_{11} smoothed')

subplot(3,2,4);
hold off
%[hhxy,bbxy]=hist(J_xy_sm(:)-J0_xy_sm(:),nn);
[hhxy,bbxy]=hist(J_xy_sm(grady<gradybunch),nn);

hhxy = hhxy/sum(hhxy)/mean(abs(diff(bbxy)));
bar(bbxy,hhxy,'c');
hold on

hhn = normpdf(bbxy,0,sqrt(2*vdx^2*vfJ12/1.1)); % varför dela med 1.1??
ll=plot(bbxy,hhn,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})
title('J_{12} smoothed')
axis([-30 30 0 0.065])


j011 = J0_xx_sm;
j012 = J0_xy_sm;
j022 = J0_yy_sm;
lam01 = j011/2 + j022/2 + (j011.^2 - 2*j011.*j022 + 4*j012.^2 + j022.^2).^(1/2)/2;
lam02 = j011/2 + j022/2 - (j011.^2 - 2*j011.*j022 + 4*j012.^2 + j022.^2).^(1/2)/2;

lam01e = j011/2 + j022/2 + abs(j012)/2+abs(j011-j022)/2;
lam02e = j011/2 + j022/2 - abs(j012)/2-abs(j011-j022)/2;


j11 = J_xx_sm;
j12 = J_xy_sm;
j22 = J_yy_sm;
lam1 = j11/2 + j22/2 + (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;
lam2 = j11/2 + j22/2 - (j11.^2 - 2*j11.*j22 + 4*j12.^2 + j22.^2).^(1/2)/2;

lam1e = j11/2 + j22/2 + abs(j12)/2+abs(j11-j22)/2;
lam2e = j11/2 + j22/2 - abs(j12)/2-abs(j11-j22)/2;


lam1e_m = 0.5*(vdx+vdx+sqrt(2/pi)*sqrt(2*vdx^2*vfJ12)+sqrt(2*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11)/pi));
vcovfac1 = 1/1.9;
lam1e_v = vcovfac1*(2*vdx^2*vfJ11+2*vdx^2*vfJ11+(pi-2)/pi*2*vdx^2*vfJ12+(pi-2)/pi*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11));

lam2e_m = 0.5*(vdx+vdx-sqrt(2/pi)*sqrt(2*vdx^2*vfJ12)-sqrt(2*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11)/pi));
vcovfac2 = 1/4;
lam2e_v = vcovfac2*(2*vdx^2*vfJ11+2*vdx^2*vfJ11+(pi-2)/pi*2*vdx^2*vfJ12+(pi-2)/pi*(2*(1-xxyysmcovf)*2*vdx^2*vfJ11));


% disp('Eigenvalues...')
% 
% disp([mean(lam1(:)) mean(lam1e(:)) lam1e_m]);
% disp([var(lam1(:)) var(lam1e(:)) lam1e_v]);
% 
% disp([mean(lam2(:)) mean(lam2e(:)) lam2e_m]);
% disp([var(lam2(:)) var(lam2e(:)) lam2e_v]);


subplot(3,2,5);
hold off

%[hhl1,bbl1]=hist(lam1(:)-lam01(:),nn);
[hhl1,bbl1]=hist(lam1(grady<gradybunch),nn);
hhl1 = hhl1/sum(hhl1)/mean(abs(diff(bbl1)));
bar(bbl1,hhl1,'c');
hold on
hhl1e = gampdf(bbl1,lam1e_m^2/lam1e_v,lam1e_v/lam1e_m);
ll=plot(bbl1,hhl1e,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})
%mb = max(bbl1);
mh = 1.3*max(hhl1(:));
%axis([0 mb 0 mh])
axis([0 100 0 mh])
title('Eigenvalue 1')


subplot(3,2,6);
hold off
%[hhl2,bbl2]=hist(lam2(:)-lam02(:),nn);
[hhl2,bbl2]=hist(lam2(grady<gradybunch),nn);
hhl2 = hhl2/sum(hhl2)/mean(abs(diff(bbl2)));
bar(bbl2,hhl2,'c');
hold on
hhl2e = gampdf(bbl2,lam2e_m^2/lam2e_v,lam2e_v/lam2e_m);
ll=plot(bbl2,hhl2e,'r');
set(ll,'LineWidth',2);
legend({'Histogram','Estimated PDF'})
mb = max(bbl2);
mh = 1.3*max(hhl2(:));
%axis([0 mb 0 mh])
axis([0 65 0 mh])

title('Eigenvalue 2')

subplot(3,2,5);
xlabel('Fig. 5. Top: The original image and the input image with added...');

%%
disp('Creating Fig. 6 ...');
figure(6);
clf
subplot(1,2,1);
xx = 0:0.1:100;
ll = plot(xx,min(4,(4-0.5)*exp(-(0.5/10)*(xx-10))+0.5));
set(ll,'Color',[0.05 0.5 0]);
set(ll,'LineWidth',2);
hold on
axis([0 105 0 4.5])
plot([16.7 16.7 0],[0 3 3],'k');
plot([48.9 48.9 0],[0 1 1],'k')
plot([0 103],[0.5 0.5],'k--');
plot([10 10],[0 4],'k--');
xlabel('Fig. 6. Depiction of the smoothing kernel construction...'); 
subplot(1,2,2);
drawkernel;

%%
disp('Creating Fig. 7 using 6 images ...');


imnames = {'data/DSC03946.jpg','data/DSC03956.jpg','data/DSC04025.jpg','data/DSC02634.jpg','data/IMG_1084.jpg','data/IMG_1265.jpg'};
nn = length(imnames);
imins = cell(1,nn);
imouts = cell(1,nn);
stds =[1.5 1 1 1.5 2 1.5];
fprintf('Processing image nr ');
for imnr = 1:nn,        
    fprintf('%d ',imnr);
    im = imread(imnames{imnr});
    imins{imnr}=im;
    im3_sm = double(im);
    for iii = 1:3,
        im_sm = structure_smooth(double(im(:,:,iii)),stds(imnr));
        im3_sm(:,:,iii)=im_sm;
    end
    imout1 = dtm_rgb(im3_sm,256,1,2000);
%        imout2 = dtm_rgb(im3_sm,256,0,2000);
    imouts{imnr}=imout1;
end
fprintf('\n');
figure(7);
imshow([imins{1} imins{2} imins{3};imouts{1} imouts{2} imouts{3};...
    imins{4} imins{5} imins{6};imouts{4} imouts{5} imouts{6}]);


%%
disp('Creating Fig. 8 ...');
figure(8);
clf
imshow([imins{1}(501:1200,1001:2000,:)*7 imouts{1}(501:1200,1001:2000,:)]);



