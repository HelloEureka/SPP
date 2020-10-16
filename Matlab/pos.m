epochPos=load('../data/ecef.txt');
%lot(epochPos(:,2:4),'*')
epochPos_v = epochPos - [0, -2341333.112, -3539049.520, 4745791.248];

plot(epochPos_v(:,2:4),'*--');
legend('X','Y','Z');
ylim([-10,10]);

