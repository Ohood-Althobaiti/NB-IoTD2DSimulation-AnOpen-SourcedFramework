%plot power consumption profile
%function plotPuser(usrr,mtrxenergy)
usrr=50;
figure()
plot(mtrxenergy(:,usrr))
title('power consumption for user 50')
xlabel('Number of rounds r')
ylabel('Energy disipated J')

