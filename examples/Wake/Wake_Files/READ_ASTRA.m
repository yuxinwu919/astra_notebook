P0=load('test.ini');
P1=load('wake.0210.001');
NP=size(P0,1);
P0(2:NP,3)=P0(1,3)+P0(2:NP,3);
P0(2:NP,6)=P0(1,6)+P0(2:NP,6);
P1(2:NP,3)=P1(1,3)+P1(2:NP,3);
P1(2:NP,6)=P1(1,6)+P1(2:NP,6);

figure(1);
Z1=P1(:,3)-mean(P1(:,3));
za=min(Z1);
zb=max(Z1);
zr=za:0.01*(zb-za):zb;
hist(Z1,zr);

figure(2);
plot(Z1,P1(:,6)-P0(:,6),'.');

T=load('test.dat');
figure(3); plot(T(:,1),T(:,2));
figure(4); plot(T(:,1),T(:,3));
