clear 
% epochPos=load('../data/hkkt/res_spp.txt');
epochPos=load('../data/albh/res_spp.txt');

sitRef = [-2405145.476;5385196.812;2420034.840];
[B_0, L_0, H_0] = XYZ2BLH(sitRef);


t = epochPos(:,1) / 3600;
d_XYZ = epochPos(:,5:7);
sigma = epochPos(:,8);
satCount = epochPos(:,9);

brdcValid = 1:round(22.5*3600/30);

dENU = ECEF2ENU(d_XYZ', B_0, L_0);

subplot(1,2,1)
% plot(t(brdcValid), d_XYZ(brdcValid,:),'*');
% legend('dX','dY','dZ');

plot(t(brdcValid), dENU(:,brdcValid)','*');
legend('dE','dN','dU');


subplot(1,2,2)
plot(t(brdcValid), sigma(brdcValid), 'o');
legend("sigma");



dENU = dENU';
rms = sqrt(sum((d_XYZ(brdcValid,:).^2))/(length(brdcValid)-1))
pos_err = norm(rms)

rms = sqrt(sum((dENU(brdcValid,:).^2))/(length(brdcValid)-1))

pos_err = norm(rms)