%%

steppy = 0.05;
[xx,yy]=meshgrid(-12:steppy:12);

s1 = 2;
s2 = 4;
th = 30/360*2*pi;

xr = cos(th)*xx+sin(th)*yy;
yr = -sin(th)*xx+cos(th)*yy;
gk = exp(-0.5*(xr.^2/s1^2+yr.^2/s2^2));

imshow(gk);
hold on

[M,N]=size(gk);

colormap(parula)


%l1=plot3([0 12*cos(th)],[0 12*sin(th)],[0.03 0.03],'g');
l1=quiver(N/2, M/2,N/2*cos(th),M/2*sin(th),'k');
set(l1,'LineWidth',2)
%l1=plot3([0 12*sin(th)],[0 -12*cos(th)],[0.03 0.03],'g');
l1=quiver(N/2, M/2, N/2*sin(th),-M/2*cos(th),'k');
set(l1,'LineWidth',2)
%axis equal
colorbar
t1 = text(N/2+N/2*cos(th)-0.3/steppy,M/2+M/2*sin(th)-0.3/steppy,'v_2');
set(t1,'FontSize',20);
set(t1,'Color','k');

t1 = text(N/2+N/2*sin(th),M/2-M/2*cos(th)+0.8/steppy,'v_1');
set(t1,'FontSize',20);
set(t1,'Color','k');


v1 = [cos(th) sin(th)];
v2 = [sin(th) -cos(th)];

x1 = [N M]/2+v1*s1/steppy+v2*0.5/steppy;
x2 = [N M]/2+v1*s1/steppy-v2*0.5/steppy;

l1 = plot([x1(1) x2(1)],[x1(2) x2(2)],'k');
set(l1,'LineWidth',3)
t1 = text(x1(1)+0.3/steppy,x1(2)-0.3/steppy,'\sigma_2');
set(t1,'FontSize',20);
set(t1,'Color','k');

x1 = [N M]/2+v2*s2/steppy+v1*0.5/steppy;
x2 = [N M]/2+v2*s2/steppy-v1*0.5/steppy;

l1 = plot([x1(1) x2(1)],[x1(2) x2(2)],'k');
set(l1,'LineWidth',3)
t1 = text(x1(1)+0.3/steppy,x1(2),'\sigma_1');
set(t1,'FontSize',20);
set(t1,'Color','k');
