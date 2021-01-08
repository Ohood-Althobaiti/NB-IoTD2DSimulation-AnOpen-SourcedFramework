% Plot Base station latency energy performance 

%for r=1:size(mtrxbtrylftma6,3)
r=3;% required round
L1=mtrxcc1v1(i,:,r);
L2=mtrxcc1v2(i,:,r);
mtrxcc2v1(i,:,r);
mtrxcc2v2(i,:,r);


%L2=mtrxbtrylftmb6(:,:,r);
range_f1=range_f1_BS6;
%-----------------------
length(range_f1)
figure()

for i=1:length(range_f1)
plot(range_c2,L1(i,:),range_c2,L2(i,:))
hold on

end
xlabel('c_2: Number of repetitions for class')
ylabel('L: Battery lifetime [days]')


str = {strcat('class 1, f1= ' , num2str(range_f1(1)))};  
str2 = {strcat('class 2, f1= ' , num2str(range_f1(1)))};  
str3 = {strcat('class 1, f1= ' , num2str(range_f1(2)))};  
str4 = {strcat('class 2, f1= ' , num2str(range_f1(2)))};  

legend(str{:},str2{:},str3{:},str4{:})
%end
title('Mutual impact among coexisting classes in a cell')
grid on
%xticks([0 5 10 15 20])
%yticks([0 500 1000 1500])