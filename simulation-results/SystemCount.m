clear

load('Gain(0,90)2deg.mat');
theta=Gain(:,3);
Gain0=Gain(:,1);
Gain90=Gain(:,2);

Kpdb=15; %db
Kp=db2pow(Kpdb);
Kndb=3.2; %db
Kn=db2pow(Kndb);
p=3;
wavlen=5e-2; %m
scanang=18; %+-deg
beamwidthx=3; %deg
beamwidthy=4; %deg
tmax=-17; %db УБЛ
k=2*pi/wavlen;

syms x;
delta=max(double(solve(-13-13*x-22*x^2==tmax)));

alph=0.5*log(0.5)/log(cosd(scanang));
thetadif=acosd(db2pow(tmax)^(1/(2*alph)));
d=wavlen/(sind(thetadif)+sind(scanang));


N1=ceil((1+0.636*delta^2)*51*wavlen/(beamwidthy*d));
xn=([1:N1]-(N1+1)/2)*d;
ampn=1+delta*cos(2*pi*xn/(N1*d));

t0=-scanang:scanang:scanang;

%%%%%%
phn=[[1:N1]-1].*d.*k.*sind(0);
In=ampn.*exp(1i.*phn);
F0=sum(In.*exp(1i*k*xn*sind(0)));
F0=abs(F0);
Fant0=F0*max(db2pow(Gain90));
%%%%%%

Gain0n=db2pow(Gain90)/max(db2pow(Gain90));

for t0ind=1:size(t0,2)
    phn=[[1:N1]-1].*d.*k.*sind(t0(t0ind));
    In=ampn.*exp(1i.*phn);
    Th=-180:1:180;
    for ind=1:length(Th);
     F(ind)=sum(In.*exp(1i*k*xn*sind(Th(ind))));
    end

    Fa(t0ind,:)=abs(F)/F0;
    Fant(t0ind,:)=Fa(t0ind,:).*(Gain0n');
    Fant1(t0ind,:)=abs(F).*db2pow(Gain90');

    pks=findpeaks(Fa(t0ind,:));
    t=2*pow2db(max(pks(pks<0.6))/max(pks));
    lgo(t0ind)=strcat("$\theta_0 =$ ",num2str(t0(t0ind)),"$^\circ$");
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);clf;hold on;
plot(Th,Fant1)
line([scanang scanang],[0 Fant0],'Color','red','LineStyle','--')
line([-scanang -scanang],[0 Fant0],'Color','red','LineStyle','--')
axis([-90 90 0 Fant0])
title('$F_1(\theta)$ when $\varphi=90^\circ$','interpreter','latex')
plot(Th,F0*db2pow(Gain90),'k:')
legend([lgo,"limits of $\theta_{scan}$","", "$f(\theta)\cdot D_{array}$"],'Interpreter','latex');

yl = ylabel('$E,times$','Interpreter','latex')
yl.Rotation=0;yl.Position=[-82 145 1];

xl=xlabel('$\theta$,deg','Interpreter','latex');
xl.Position = [90 -5 1];