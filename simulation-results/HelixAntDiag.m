load('Gain(0,90)2deg.mat');
theta=Gain(:,3);
Gain0=Gain(:,1);
Gain90=Gain(:,2);
figure(1)
ax=polaraxes;
polarplot(theta,(Gain0));
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
cirdeg=-pi/3:pi/180:pi/3;
hold on
polarplot(theta,(Gain90));
rlim([ceil(abs(min([Gain0; Gain90]))/5)*5*sign(min([Gain0; Gain90])) ceil(abs(max([Gain0; Gain90]))/5)*5*sign(max([Gain0; Gain90]))])
r=repmat(max(Gain0)-3,1,length(cirdeg));
polarplot(cirdeg,r)
ind=find(Gain0==max(Gain0));
polarplot(theta(ind),(Gain0(ind)),'r*')
title('Helix far filed')
pllg=["Phi = 0 deg","Phi = 90 deg"];
lgo=legend(pllg);
lgo.NumColumns=2;



figure(2)
[Th,Ph]=meshgrid(theta);
G=sqrt(db2pow(Gain0)*db2pow(Gain90)');
[x,y,z]=sph2cart(Th,Ph,G);
mesh(x,y,z)