%--------------------------------------------------------------------------
% Date: April, 2020
%
%  A. Azari et al. 
%             "On the Latency-Energy Performance of NB-IoT Systems
%              in Providing Wide-Area IoT Connectivity"
%
% 
%--------------------------------------------------------------------------
function [mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2]=LntcyEnrgyPrfrmnc(N,u,tau)

%% Table I: Parameters for performance analysis. 
% Indices 1 and 2 refer to coverage class 1 and 2 respectively.
% -----------------
% Traffic
% -----------------
%N = 20000;      % nodes, number of devices
S = 0.5*24;     % packet generation frequency, 1 packet per 2 hours [per day]
pu = 0.8;       % probability of uplink (res. DL) service request
% moments of uplink and downlink packet lengths la, ma
la = 500;       % [bits]
ma = 5000;      % [bits]
% average length of control and RA signaling

CF = 0.01;      % length of communication frame [sec]
lambda = 1/CF;  % frequency of arrival of BS-initiated control data [per sec]
f1 = 0.5;       % fraction of devices belongs to each coverage class
%f2 = 1-f1;      % 0.5;
% -----------------
% Coverage
% -----------------
c1 = 1;         % number of repetitions         
c2 = 2;       
Dsyn1 = 0.33;   % synchronization delay [sec]    
Dsyn2 = 0.66;       
R1 = 5000;      % uplink data rate [bps]
R2 = 5000;
dR1 = 15000;    % downlink data rate [bps]
dR2 = 15000;
% -----------------
% RRM
% -----------------
b = 0.2;   % fraction of time in which reference signals are scheduled
Tth = 2;   % maximum waiting for receiving RAR message [sec]
M = 16;    % number of RA resources M1; M2 16 preambles
% time interval between two scheduling of NPRACH t design parameter
% time interval between two scheduling of NPCCCH d design parameter
% -----------------
% Other
% -----------------
E0 = 1000; % Device’s battery capacity [joule]
% Device’s power consumption
Pt = 0.2;   % device’s power consumption in transmission        [W]
PI = 0.01;  % device’s power consumption in idle                [W]
Pl = 0.1;   % device’s power consumption in listening           [W]
Pc = 0.01;  % device’s power consumption in electronic circuits [W]
Pt1 = Pt; 
Pt2 = Pt; 
eps = 1;      % power amplifier efficiency, no value provided
Nrmax = 1;    % maximum number of allowed attempts, no value provided
C = 2;        % number of classes
%% Figure 5
range_t = 0.065;
range_d = 0.01;
tau_ex5 = 0.002;
range_f1 = [0.9 0.95 0.98 0.995];
range_c2 = 1 : 2 : 19;
L1 = zeros(4,length(range_c2));
L2 = zeros(4,length(range_c2));
for i = 1 : length(range_f1)
    f1_ex5 = range_f1(i);
    for k = 1 : length(range_c2)
        c2_ex5 = range_c2(k);
        [sLT1,sLT2,~,~,~,~,~,~,~,~,~] = ...
                    simulator(N,S,pu,la,ma,u,tau_ex5,lambda,f1_ex5,c1,c2_ex5,C,...
                              Dsyn1,Dsyn2,R1,R2,dR1,dR2,b,Tth,M,...
                              E0,PI,Pl,Pc,Pt1,Pt2,eps,Nrmax,range_t,range_d,1);
        L1(i,k) = sLT1;
        L2(i,k) = sLT2;
    end
end
figure(6)
plot(range_c2,L1(1,:),'b-+',range_c2,L2(1,:),'r-o',...
     range_c2,L1(2,:),'y-*',range_c2,L2(2,:),'m',...
     range_c2,L1(3,:),'g-x',range_c2,L2(3,:),'c-s',...
     range_c2,L1(4,:),'r-d',range_c2,L2(4,:),'b')
xlabel('c_2: Number of repeatitions for class 2')
ylabel('L: Battery lifetime [days]')
legend('Class 1, f1 = 0.9',  'Class 2, f1 = 0.9',   ...
       'Class 1, f1 = 0.95', 'Class 2, f1 = 0.95',  ...
       'Class 1, f1 = 0.98', 'Class 2, f1 = 0.98',  ...
       'Class 1, f1 = 0.995','Class 2, f1 = 0.995', ...
       'Location', 'southwest')
title('Mutual impact of two coexisting classes')
grid on
xticks([0 5 10 15 20])
yticks([0 500 1000 1500])

%% Figure 8
fctrngt=.025;
fctrngd=.025;


%range_t = [0.05:0.005:0.1,0.1:0.005:1];
%range_d = [0.001:0.001:0.01,0.01:0.001:0.09,0.1:0.001:1];
range_t = [0.05:fctrngt:1];
range_d = [0.001:fctrngt:1];

[sLT1,sLT2,sDu1,sDu2,sDd1,sDd2,mEu1,mEu2,mEd1,mEd2] = ...
                simulator(N,S,pu,la,ma,u,tau,lambda,f1,c1,c2,C,...
                          Dsyn1,Dsyn2,R1,R2,dR1,dR2,b,Tth,M,...
                          E0,PI,Pl,Pc,Pt1,Pt2,eps,Nrmax,range_t,range_d,0);
               
[mtrxtt,mtrxdt] = meshgrid(range_t,range_d);
assignin('base','mtrxtt',mtrxtt)
assignin('base','mtrxdt',mtrxdt)

figure(7)
s1=mesh(mtrxdt,mtrxtt,sLT1);
%s1.FaceColor = [0 0 1];
%s1.FaceAlpha = 0.6;
view(3)
grid
xlabel('td(sec)')
ylabel('t(sec)')
zlabel('Battery lifetime(days)')
hold on
s2=mesh(mtrxdt,mtrxtt,sLT2);
%s2.FaceColor = [1 0 0];
%s2.FaceAlpha = 0.6;
assignin('base','mtrxdt',mtrxdt)
assignin('base','mtrxtt',mtrxtt)
assignin('base','sLT1',sLT1)
assignin('base','sLT2',sLT2)
title('Battery lifetime L versus t and td')
txt1 = {'Class 2'};
text(4,0.5,600,txt1)
txt2 = {'Class 1'};
text(4,0.5,1240,txt2)
v=[-3 -4 1];%view point reference
[caz,cel] = view(v);%change view point
set(gca,'xscale','log')
set(gca,'yscale','log')
%set(gca,'zscale','log')
hold off
m = max(max(sLT1));
[id,it]=find(sLT1==m);
fprintf('Fig.81: Maximum lifetime for class 1 is @ d=%.3f and t=%.2f. \n',range_d(id),range_t(it))
m = max(max(sLT2));
[id,it]=find(sLT2==m);
fprintf('Fig.81: Maximum lifetime for class 2 is @ d=%.3f and t=%.2f. \n',range_d(id),range_t(it))

figure(8)
mesh(mtrxdt,mtrxtt,sDu1)
%view([40,20])
grid
xlabel('td(sec)')
ylabel('t(sec)')
zlabel('Uplink latency(sec)')
hold on
mesh(mtrxdt,mtrxtt,sDu2)
hold off
assignin('base','sDu1',sDu1)
assignin('base','sDu2',sDu2)
title('Uplink latency Du versus t and td')
txt1 = {'Class 2'};
text(4,0.5,0.85,txt1)
txt2 = {'Class 1'};
text(4,0.5,0.35,txt2)
set(gca,'xscale','log')
set(gca,'yscale','log')
%set(gca,'zscale','log')

m = min(min(sDu1));
[id,it]=find(sDu1==m);
fprintf('Fig.82: Minimum uplink latency for class 1 is @ d=%.3f and t=%.2f. \n',range_d(id),range_t(it))
m = min(min(sDu2));
[id,it]=find(sDu2==m);
fprintf('Fig.82: Minimum uplink latency for class 2 is @ d=%.3f and t=%.2f. \n',range_d(id),range_t(it))

figure(9)
mesh(mtrxdt,mtrxtt,mEu1)
%view([40,20])
grid
xlabel('td(sec)')
ylabel('t(sec)')
zlabel('Uplink Energy(J)')
hold on
mesh(mtrxdt,mtrxtt,mEu2)
hold off
assignin('base','mEu1',mEu1)
assignin('base','mEu2',mEu2)

title('Energy consumption in a Uplink service (J)')
txt1 = {'Class 2'};
text(4,0.5,1.7,txt1)
txt2 = {'Class 1'};
text(4,0.5,0.01,txt2)
set(gca,'xscale','log')
set(gca,'yscale','log')
%set(gca,'zscale','log')


figure(10)
mesh(mtrxdt,mtrxtt,sDd1)
assignin('base','sDd1',sDd1)
%view([170,20])
grid
xlabel('td(sec)')
ylabel('t(sec)')
zlabel('Downlink latency(sec)')
hold on
mesh(mtrxdt,mtrxtt,sDd2)
assignin('base','sDd2',sDd2)
hold off
title('Downlink latency Dd versus t and td')
txt1 = {'Class 2'};
text(4,0.5,1.9,txt1)
txt2 = {'Class 1'};
text(4,0.5,0.9,txt2)
v=[-3 -4 2.5];%view point reference
[caz,cel] = view(v);%change view point
set(gca,'xscale','log')
set(gca,'yscale','log')
%set(gca,'zscale','log')
%zlim([1 5])
%zticks([2 3 4 5])



m = min(min(sDd1));
[id,it]=find(sDd1==m);
fprintf('Fig.83: Minimum downlink latency for class 1 is @ d=%.3f and t=%.2f. \n',range_d(id),range_t(it))
m = min(min(sDd2));
[id,it]=find(sDd2==m);
fprintf('Fig.83: Minimum downlink latency for class 2 is @ d=%.3f and t=%.2f. \n',range_d(id),range_t(it))



figure(11)
mesh(mtrxdt,mtrxtt,mEd1)
%view([40,20])
grid
xlabel('td(sec)')
ylabel('t(sec)')
zlabel('Downlink Energy (J)')
hold on
mesh(mtrxdt,mtrxtt,mEd2)
hold off
assignin('base','mEd1',mEd1)
assignin('base','mEd2',mEd2)

title('Energy consumption in a Downlink service (J)')
txt1 = {'Class 2'};
text(4,0.5,0.14,txt1)
txt2 = {'Class 1'};
text(4,0.5,0.09,txt2)
%v=[-3 -4 2.5];%view point reference
%[caz,cel] = view(v);%change view point
set(gca,'xscale','log')
set(gca,'yscale','log')
%set(gca,'zscale','log')




end

function [sLT1,sLT2,sDu1,sDu2,sDd1,sDd2,sE1,mEu1,mEu2,mEd1,mEd2] = ...
         simulator(N,S,pu,la,ma,u,tau,lambda,f1,c1,c2,C,Dsyn1,Dsyn2,R1,R2,dR1,dR2,b,Tth,M,...
                   E0,PI,Pl,Pc,Pt1,Pt2,eps,Nrmax,range_t,range_d,ctrl)
    
    f2 = 1 - f1;
    
    % The arrival rates of uplink/downlink service requests per day, [2] (1)
    Gu = N*S*pu/(24*3600); 
    Gd = N*S*(1-pu)/(24*3600);

    sLT1 = nan(length(range_d),length(range_t));
    sLT2 = nan(length(range_d),length(range_t));
    sDu1 = nan(length(range_d),length(range_t));
    sDu2 = nan(length(range_d),length(range_t));
    sDd1 = nan(length(range_d),length(range_t));
    sDd2 = nan(length(range_d),length(range_t));
    sE1  = nan(length(range_d),length(range_t));

    % Calculate cumulative distribution function values for different n
    i = 1;
    mQ_tot = zeros(length(range_d)*length(range_t),1);
    for td = range_d 
        for t = range_t
            mQ_tot(i) = (Gu+Gd)*max(td,t)+lambda*td; 
            i = i + 1;
        end
    end
    K_max = 2*floor(max(mQ_tot));
    Fcdf_tot = ones(K_max-1,1);
    %disp('CDF calculation');
    if ctrl
        for n = 1 : K_max - 1
           try
               Fcdf_tot(n) = Fcdf1(n,u,Tth,f1,f2,c1,c2);
           catch
               Fcdf_tot(n) = 1;
           end
           %n
        end
    %else
        %for n = 1 : K_max - 1
        %    Fcdf_tot(n) = Fcdf2(n,u,Tth,f1,f2,c1,c2);
        %    %n
        %end
    end
    
    ii = 0;
    id = 0;
    for td = range_d 
        id = id + 1;
        it = 0; 
        for t = range_t
            it = it + 1;
            ii = ii + 1; 

            % [2] (6)
            mQ = (Gu+Gd)*max(td,t)+lambda*td;  % (f1 + f2 = 1) 

            % [2] (7)
            bDt = (f1*c1 + f2*c2) * u;

            % [2] (8)
            Drar1 = 0.5*td + 0.5*mQ*bDt + c1*u;
            Drar2 = 0.5*td + 0.5*mQ*bDt + c2*u;

            Dra1 = 0.5*t + c1*tau;
            Dra2 = 0.5*t + c2*tau;

            % [2] (9)
            Prach1 = 0;        
            Prach2 = 0;
            N1 = f1*(Gu+Gd)*td;%t; 
            N2 = f2*(Gu+Gd)*td;%t;
            for k = 2 : N
                Prach1 = Prach1 + ((N1^k)*exp(-N1)/factorial(k)) * ((M-1)/M)^(k-1);
                Prach2 = Prach2 + ((N2^k)*exp(-N2)/factorial(k)) * ((M-1)/M)^(k-1);
            end
            
            % [2] (10-11) 
            X = 0;
            for K = 2 : 2*floor(mQ)
                for k = 1 : K-1
                    Fc1 = 1 - Fcdf_tot(K-k);   % Fcdf(K-k,u,Tth,f1,f2,c1,c2)
                    if K-k-1 > 0
                        Fc2 = Fcdf_tot(K-k-1); % Fcdf(K-k-1,u,Tth,f1,f2,c1,c2)
                    else
                        Fc2 = 1;
                    end
                    X = X + k/K * ((mQ^K)*exp(-mQ))/factorial(K) * Fc1 * Fc2; 
                end
            end
            Pjrar = 1 - X;  

            Pj1 = Prach1 * Pjrar; 
            Pj2 = Prach2 * Pjrar; 
            % [2] (5)
            Pjl1 = 0;
            Pjl2 = 0;
            for l = 1 : Nrmax
                Pjl1 = Pjl1 + ((1 - Pj1)^(l-1))*Pj1*l;
                Pjl2 = Pjl2 + ((1 - Pj2)^(l-1))*Pj2*l;
            end
            
            Drr1 =  Pjl1 * (Dra1 + Drar1);
            Drr2 =  Pjl2 * (Dra2 + Drar2);

            % [2] (12)
            w = 1 - (c1 + c2)*tau/t;    

            % [2] (13)
            sa = f1*c1*la/(R1*w) + f2*c2*la/(R2*w); 
            lb = la^2;
            sb = f1*c1^2*lb/(R1^2*w^2) + f2*c2^2*lb/(R2^2*w^2);

            G = Gu*t/C;    % f1 + f2 = 1
            ro = C*G*sa/t; 
            if(ro > 0.99 || ro < 0 || isnan(ro)|| isinf(ro))
                continue     
            end   
            
            % [2] (14)
            Dtx1 = 0.5*ro*sb/(sa*(1-ro)) + 0.5*G*sa/(1-ro) + c1*la/(R1*w);
            Dtx2 = 0.5*ro*sb/(sa*(1-ro)) + 0.5*G*sa/(1-ro) + c2*la/(R2*w);

            % [2] (15)
            y = 1 - b - mQ*bDt/td;

            % [2] (16)
            ha = f1*c1*ma/(dR1*y) + f2*c2*ma/(dR2*y);  
            mb = ma^2;
            hb = f1*c1^2*mb/(dR1^2*y^2) + f2*c2^2*mb/(dR2^2*y^2);  

            % [2] (17)
            G = Gd*t/C;    % f1 + f2 = 1
            eta = C*G*ha/t;      
            if(eta > 0.99 || eta < 0 || isnan(eta)|| isinf(eta))
                continue
            end                     
            Drx1 = 0.5*eta*hb/(ha*(1-eta)) + 0.5*G*ha/(1-eta) + c1*ma/(dR1*y);
            Drx2 = 0.5*eta*hb/(ha*(1-eta)) + 0.5*G*ha/(1-eta) + c2*ma/(dR2*y);
            %%
            % [2] (2)
            Du1 = Dsyn1 + Drr1 + Dtx1; 
            Du2 = Dsyn2 + Drr2 + Dtx2;
            Dd1 = Dsyn1 + Drr1 + Drx1; 
            Dd2 = Dsyn2 + Drr2 + Drx2;
            
            
            

            % [2] (18)
            Esyn1 = Pl*Dsyn1;
            Esyn2 = Pl*Dsyn2;

            % [2] (19)
            Erar1 = Pl*Drar1;
            Erar2 = Pl*Drar2;

            % [2] (21)
            Era1 = (Dra1 - c1*tau)*PI + c1*tau*(Pc + eps*Pt1);
            Era2 = (Dra2 - c2*tau)*PI + c2*tau*(Pc + eps*Pt2);

            % [2] (20)
            Pjl1 = 0;
            Pjl2 = 0;
            for l = 1 : Nrmax
                Pjl1 = Pjl1 + ((1 - Pj1)^(l-1))*Pj1;
                Pjl2 = Pjl2 + ((1 - Pj2)^(l-1))*Pj2;
            end
            
            Err1 = Pjl1*(Era1 + Erar1); 
            Err2 = Pjl2*(Era2 + Erar2);

            % [2] (22)
            Etx1 = (Dtx1 - c1*la/(R1*w))*PI + (Pc + eps*Pt1)*c1*la/(R1*w);
            Etx2 = (Dtx2 - c2*la/(R2*w))*PI + (Pc + eps*Pt2)*c2*la/(R2*w);

            % [2] (23)
            Erx1 = (Drx1 - c1*ma/(dR1*y))*PI + Pl*c1*ma/(dR1*y);
            Erx2 = (Drx2 - c2*ma/(dR2*y))*PI + Pl*c2*ma/(dR2*y);
            
            % Eq (24)
                % Eq (3)
                Eu1 = Esyn1 + Err1 + Etx1;
                Eu2 = Esyn2 + Err2 + Etx2;
                Ed1 = Esyn1 + Err1 + Erx1;
                Ed2 = Esyn2 + Err2 + Erx2;
                % Eq (4)
                LT1 = E0/(S*pu*Eu1 + S*(1-pu)*Ed1);
                LT2 = E0/(S*pu*Eu2 + S*(1-pu)*Ed2);

            sLT1(id,it) = LT1;
            sLT2(id,it) = LT2; 

            sDu1(id,it) = Du1;  
            sDu2(id,it) = Du2; 

            sDd1(id,it) = Dd1;
            sDd2(id,it) = Dd2;
            
            mEu1(id,it) = Eu1;  
            mEu2(id,it) = Eu2;
            
            mEd1(id,it) = Ed1;  
            mEd2(id,it) = Ed2;
            

            sE1(id,it) = S*Eu1 + S*Ed1;
        end
    end
end
% [2] (10)
function y = Fcdf1(n,u,x,f1,f2,c1,c2)
   if n == 1
      y = 0;
      if c1*u < x
          y = y + f1;
      end
      if c2*u < x
          y = y + f2;
      end
   else
      y = f1*Fcdf1(n-1,u,x-c1*u,f1,f2,c1,c2) + f2*Fcdf1(n-1,u,x-c2*u,f1,f2,c1,c2);
   end
 end

% If n < 500 and parameters in the Table I, this is alternative for Fcdf
function y = Fcdf2(n,u,x,f1,f2,c1,c2)
    A = [];
    for i = 0 : 1 : n
        A = [A, nchoosek(n,i)];
    end
    if x > n*c1*u
        F1 = f1*ones(1,n+1);
    else
        F1 = zeros(1,n+1);
    end
    if x > n*c2*u
        F2 = f2*ones(1,n+1);
    else
        F2 = zeros(1,n+1);
    end  
    y = (F1.^(n:-1:0)) .* (F2.^(0:1:n)) * A'; 
end