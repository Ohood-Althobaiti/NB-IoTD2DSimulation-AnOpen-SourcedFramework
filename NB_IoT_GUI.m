function varargout = NB_IoT_GUI(varargin)
% NB_IOT_GUI MATLAB code for NB_IoT_GUI.fig
%      NB_IOT_GUI, by itself, creates a new NB_IOT_GUI or raises the existing
%      singleton*.
%
%      H = NB_IOT_GUI returns the handle to a new NB_IOT_GUI or the handle to
%      the existing singleton*.
%
%      NB_IOT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NB_IOT_GUI.M with the given input arguments.
%
%      NB_IOT_GUI('Property','Value',...) creates a new NB_IOT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NB_IoT_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NB_IoT_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NB_IoT_GUI

% Last Modified by GUIDE v2.5 30-Jul-2020 11:54:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NB_IoT_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NB_IoT_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% simulation NB IoT  internet of things


% References: [1]  G. Piro, L. Griego,F. Capozzi, P Camarda. Simulation LTE Cellular Systems:
%                  an open source framework.(G. Piro, “LTE-Sim - the LTE simulator,” [Online] Available: http://telematics.poliba.it/LTE-Sim)

%             [2]  A. Azari, C Stevanovic, P. Popoovski, C. Cavdar. On the Latency-Energy 
%                  Performance of NB-IoT Systems in Providing Wide-Area IoT Connectivity

%             [3]  F. Aderohunmu,  Energy Management Techniques in Wireless Sensor Networks: Protocol
%                  Design and Evaluation

%             [4]  T. Singh. Analysis of Low Energy Adaptive Clustering Hierarchy (LEACH) protocol. 
%                  Department of Computer Science and Engineering National Institute of Technology Rourkela. 
%  
%             [5]  M. Renato, Network planning Model for NB-IoT. Faculdade de engenharia, universidade do porto.

%close all;

% --- Executes just before NB_IoT_GUI is made visible.
function NB_IoT_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NB_IoT_GUI (see VARARGIN)

% Choose default command line output for NB_IoT_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NB_IoT_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NB_IoT_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

set(handles.pushbutton1,'Enable','off');  
set(handles.radiobutton3,'Enable','off');  
set(handles.radiobutton4,'Enable','off');
set(handles.radiobutton5,'Enable','off');
BS1=0;
BS2=0;
BS3=0;
BS4=0;
BS5=0;
BS6=0;
BS7=0;
BS8=0;
BS9=0;
BS10=0;
BS11=0;
BS12=0;
BS13=0;
BS14=0;
BS15=0;
BS16=0;
BS17=0;
BS18=0;
BS19=0;
BS20=0;
BS21=0;
BS22=0;
BS23=0;
BS24=0;



assignin('base','BS1',BS1);
assignin('base','BS2',BS2);
assignin('base','BS3',BS3);
assignin('base','BS4',BS4);
assignin('base','BS5',BS5);
assignin('base','BS6',BS6);
assignin('base','BS7',BS7);
assignin('base','BS8',BS8);
assignin('base','BS9',BS9);
assignin('base','BS10',BS10);
assignin('base','BS11',BS11);
assignin('base','BS12',BS12);
assignin('base','BS13',BS13);
assignin('base','BS14',BS14);
assignin('base','BS15',BS15);
assignin('base','BS16',BS16);
assignin('base','BS17',BS17);
assignin('base','BS18',BS18);
assignin('base','BS19',BS19);
assignin('base','BS20',BS20);
assignin('base','BS21',BS21);
assignin('base','BS22',BS22);
assignin('base','BS23',BS23);
assignin('base','BS24',BS24);

stp=0;
assignin('base','stp',stp);
psckbx=0;

n=str2num(get(handles.edit1,'String'));
leach=get(handles.checkbox2,'Value');
%n=5000;%Number of Nodes in the field
u = 0.002;      % [sec]
tau = 0.01;     % [sec]



%% eNB characteristics
%Established power classes for NB IoT PTXjdB: 14, 20, 23 dBm

%% Celll network topology
% cell dimensions in kilometer
Rcell=2;% Cell radius [km] [1](p5)
Rcell2=Rcell*2;
Rcell3=Rcell*3;
Rcell4=Rcell*4;
Rcell5=Rcell*5;

h=.00;


% The network topology is compsed by a set of cells and network nodes
% (including eNBs, one or more MME/GW, and UEs)

mtrxlxyz=[-Rcell2   Rcell2      h; % 1
             0      Rcell2      h; % 2
           Rcell2   Rcell2      h; % 3
           Rcell4   Rcell2      h; % 4
          -Rcell    Rcell       h; % 5
           Rcell    Rcell       h; % 6
           Rcell3   Rcell       h; % 7
           Rcell5   Rcell       h; % 8
          -Rcell2     0         h; % 9
              0       0         h; % 10
           Rcell2     0         h; % 11
           Rcell4     0         h; % 12
          -Rcell    -Rcell      h; % 13
           Rcell    -Rcell      h; % 14
           Rcell3   -Rcell      h; % 15
           Rcell5   -Rcell      h; % 16
          -Rcell2  -Rcell2      h; % 17
              0    -Rcell2      h; % 18
          Rcell2   -Rcell2      h; %19
          Rcell4   -Rcell2      h; % 20
         -Rcell   -Rcell3       h; % 21     
          Rcell   -Rcell3       h; % 22    
          Rcell3   -Rcell3      h; % 23   
          Rcell5   -Rcell3      h; % 24   
          ];

% Center the grid in the plot and determine the environment
xmin=min(mtrxlxyz(:,1));
xmax=max(mtrxlxyz(:,1));

ymin=min(mtrxlxyz(:,2));
ymax=max(mtrxlxyz(:,2));
 
ypls=((abs(ymax)+abs(ymin))/2)+ymin;
xpls=((abs(xmax)+abs(xmin))/2)+xmin;

mtrxlxyz(:,1)=mtrxlxyz(:,1)-xpls;
mtrxlxyz(:,2)=mtrxlxyz(:,2)-ypls;

xmin2=round(min(mtrxlxyz(:,1)));
xmax2=round(max(mtrxlxyz(:,1)));

ymin2=round(min(mtrxlxyz(:,2)));
ymax2=round(max(mtrxlxyz(:,2)));

% Environment dimensions based in Rcell

lx=xmax2*2+xmax/4; ly=ymax2*2+ymax/4; lz=.03;% dimensions x,y,z (km)

%---------------------------------------------------------------

X1=1; Y1=2;

eNB=size(mtrxlxyz,1);

ele_nod=[% Cell diagram nodes
1 5
5 6
6 3
3 4
4 8
9 5
6 11
11 12
8 12
9 13
13 14
11 14
12 16
17 13
14 19 
19 20 
16 20
17 21
21 22
19 22
20 24
];

num_nod=size(mtrxlxyz,1);%
num_nod2=num_nod-1;

num_ele=size(ele_nod,1);
LaG=ele_nod(:,[1,2]);

nbar=size(LaG,1);%Number de EFs (LaG)

fctr=2;
XT=mtrxlxyz(1); YT=mtrxlxyz(2);%eNB position
Nx=lx*fctr; Ny=ly*fctr;
%
gx=-lx/2:lx/Nx:lx/2;
gy=-ly/2:ly/Ny:ly/2;
[rx,ry]=meshgrid(gx,gy);
eNB=size(mtrxlxyz,1);
sz1=size(rx,1);
sz2=size(rx,2);

for k=1:eNB
lxi=mtrxlxyz(k,1); lyi=mtrxlxyz(k,2); lzi=mtrxlxyz(k,3); 
%Position
d_eNB(:,:,k)= sqrt(((abs(rx-lxi)).^2)+((abs(ry-lyi)).^2));%+((lzi)^2)); [km]
%distance between the eNB and the UE in km
end

%alphai=-(log(lamdai)/thaoi)%% alpha i

%% propagation loss model
%MijdB % Mutipath loss

rdbttn1=get(handles.radiobutton3,'Value');
rdbttn2=get(handles.radiobutton4,'Value');
rdbttn3=get(handles.radiobutton5,'Value');
if rdbttn1==1
CellScenario='Urban';% Cell scenario= 'Urban' 'Rural' 'Micro'
end
if rdbttn2==1
CellScenario='Rural';% Cell scenario= 'Urban' 'Rural' 'Micro'
end
if rdbttn3==1
CellScenario='Micro';% Cell scenario= 'Urban' 'Rural' 'Micro'
end


switch CellScenario% path loss model between BS and CH depends of the cell scenario [1](table III)[1], [25, 5, 27, 28] %@@@@@@@
    case 'Urban'% Macro cell - Urban and Suburban areas
     LidB=103.8+37.6*log10(d_eNB); % path loss %@@@@@@@
    case 'Rural'% Macro cell rural area
     LidB=96+34.1*log10(d_eNB);  %@@@@@@@
    case 'Micro'% Micro dcell
     LidB=24+45*log10(d_eNB+16); %@@@@@@@
end








x=1;


%%-----------------------------------------------------------------
Gtx = 21.10; % anntena gain
Grx = 19.18; % anntena gain

TidB=10; % penetration loss (dB) [1] (15)
% SijdB % shadowing
% eNB transmission power and reception power in the i-th subchanel
% respectively
%  PRXijdB=PTXjdB-MijdB-ttlLidB-TidB-SijdB% [1] (7)
PTXjdB=14;% Transmission power
PRXijdB=abs(PTXjdB*d_eNB-LidB-TidB);% [1] (7) % Reception power
F=2.5;
No=-174;% dBm [1] (p10) Thermal Noise Power
Bj=180000;% 180 kHz [2] (p1)
I=.50;% Interference
% the physical layer computes for each sub channel the SINR for the 
% received signal considering the received power, the noise and the interference 

SINRij=PRXijdB/(F*No*Bj+I);%[1] (8)
Nsfigr=1;% Noise figure [watts]
NsSpctrldnsty=No*Nsfigr;
SCS=15000;%15 kHz SubCarrier-spacing
vTones=[12 6 3 2 1];%[5] (3.10) %@@@@@@@
Tones=vTones(3);
f=12/Tones;%[5] (3.10) 
BW=180000/f;% Effective bandwith SCS*Tones;%
TxSpctrlDnsty=PRXijdB/BW;%[5] (3.10)


SNR=TxSpctrlDnsty/NsSpctrldnsty; % [5] %@@@@@@@

effctvNse=10*log10(BW)-No;

vMCL=[144 155 164];% [dB] Maximum coupling loss


%SNRtarget= PTXjdB - effctvNse - vMCL; % [5] (3.7)

assignin('base','LidB',LidB);
assignin('base','PRXijdB',PRXijdB);
assignin('base','SINRij',SINRij);
assignin('base','SNR',SNR);
assignin('base','effctvNse',effctvNse);
assignin('base','BW',BW);
 

ttlLidB=zeros(size(d_eNB,1),size(d_eNB,2));
ttlLidBtargt=zeros(size(d_eNB,1),size(d_eNB,2));
for k=1:eNB
ttlLidB=ttlLidB+LidB(:,:,k)/eNB;
end


figure(1)
hold off
meshc(rx,ry,ttlLidB);
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Pathloss dB');

switch CellScenario
    case 'Urban'
        title({'Average Pathloss [dB] Urban Cell scenario'});
    case 'Rural'
        title({'Average Pathloss [dB] Rural Cell scenario'});
    case 'Micro'  
        title({'Average Pathloss [dB] Micro Cell scenario'});
end
%zticks([vMCL(1) min(min(ttlLidB))  vMCL(2) vMCL(3) max(max(ttlLidB))])%@@@@@@@@@@@@




figure(2)
hold off
for k=1:eNB
hold on
subplot(6,4,k); 
meshc(rx,ry,(LidB(:,:,k)));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Path loss dB');
% zticks([100 min(min(LidB(:,:,k))) vMCL(1) vMCL(2) vMCL(3) max(max(LidB(:,:,k)))])
switch k
    case 1
      title({'BS21'})  
    case 2
      title({'BS22'})        
    case 3
      title({'BS23'})  
    case 4
      title({'BS24'}) 
    case 5
      title({'BS17'})  
    case 6
      title({'BS18'})        
    case 7
      title({'BS19'})  
    case 8
      title({'BS20'})
    case 9
      title({'BS13'})  
    case 10
      title({'BS14'})        
    case 11
      title({'BS15'})  
    case 12
      title({'BS16'}) 
    case 13
      title({'BS9'})  
    case 14
      title({'BS10'})        
    case 15
      title({'BS11'})  
    case 16
      title({'BS12'})
    case 17
      title({'BS5'})
    case 18
      title({'BS6'})  
    case 19
      title({'BS7'})        
    case 20
      title({'BS8'})  
    case 21
      title({'BS1'}) 
    case 22
      title({'BS2'})  
    case 23
      title({'BS3'})        
    case 24
      title({'BS4'})  
end

end

switch CellScenario
    case 'Urban'
        suptitle({'Path loss [dB] Urban Cell scenario'});
    case 'Rural'
        suptitle({'Path loss [dB] Rural Cell scenario'});
    case 'Micro'  
        suptitle({'Path loss [dB] Micro Cell scenario'});
end


figure(3)
hold off
meshc(rx,ry,(LidB(:,:,15)));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Pathloss [dB]');
minvle=min(min(LidB(:,:,15)));
mxvlue=max(max(LidB(:,:,15)));
%zticks([100 min(min(LidB(:,:,21))) vMCL(1) vMCL(2) vMCL(3) max(max(LidB(:,:,21)))])
switch CellScenario
    case 'Urban'
        title({'Pathloss BS15[dB] Urban Cell scenario'});
    case 'Rural'
        title({'Pathloss BS15 [dB] Rural Cell scenario'});
    case 'Micro'  
        title({'Pathloss BS15 [dB] Micro Cell scenario'});
end


figure(4)
hold off
for k=1:eNB
hold on
subplot(6,4,k); 
meshc(rx,ry,(SNR(:,:,k)));
zticks([-12:4:-2])
zticklabels({'-10','-4'})
xlabel('X (km)');
ylabel('Y (km)');
zlabel('SNR dB');


switch k
    case 1
      title({'BS21'})  
    case 2
      title({'BS22'})        
    case 3
      title({'BS23'})  
    case 4
      title({'BS24'}) 
    case 5
      title({'BS17'})  
    case 6
      title({'BS18'})        
    case 7
      title({'BS19'})  
    case 8
      title({'BS20'})
    case 9
      title({'BS13'})  
    case 10
      title({'BS14'})        
    case 11
      title({'BS15'})  
    case 12
      title({'BS16'}) 
    case 13
      title({'BS9'})  
    case 14
      title({'BS10'})        
    case 15
      title({'BS11'})  
    case 16
      title({'BS12'})
    case 17
      title({'BS5'})
    case 18
      title({'BS6'})  
    case 19
      title({'BS7'})        
    case 20
      title({'BS8'})  
    case 21
      title({'BS1'}) 
    case 22
      title({'BS2'})  
    case 23
      title({'BS3'})        
    case 24
      title({'BS4'})  
end
end
switch CellScenario
    case 'Urban'
        suptitle({'SNR [dB] Urban Cell scenario'});
    case 'Rural'
        suptitle({'SNR [dB] Rural Cell scenario'});
    case 'Micro'  
        suptitle({'SNR [dB] Micro Cell scenario'});
end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
figure(300)
hold off
meshc(rx,ry,(SNR(:,:,15)));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('SNR [dB]');
minvle=min(min(SNR(:,:,15)));
mxvlue=max(max(SNR(:,:,15)));
%zticks([100 min(min(LidB(:,:,21))) vMCL(1) vMCL(2) vMCL(3) max(max(LidB(:,:,21)))])
switch CellScenario
    case 'Urban'
        title({'SNR BS15 [dB] Urban Cell scenario'});
    case 'Rural'
        title({'SNR BS15 [dB] Rural Cell scenario'});
    case 'Micro'  
        title({'SNR BS15 [dB] Micro Cell scenario'});
end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


figure(5)
hold off
for k=1:eNB
hold on
subplot(6,4,k); 
meshc(rx,ry,(SINRij(:,:,k)));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('SINRij dB');
%zticks([-10 -4])
%zticklabels({'-10','-4'})

switch k
    case 1
      title({'BS21'})  
    case 2
      title({'BS22'})        
    case 3
      title({'BS23'})  
    case 4
      title({'BS24'}) 
    case 5
      title({'BS17'})  
    case 6
      title({'BS18'})        
    case 7
      title({'BS19'})  
    case 8
      title({'BS20'})
    case 9
      title({'BS13'})  
    case 10
      title({'BS14'})        
    case 11
      title({'BS15'})  
    case 12
      title({'BS16'}) 
    case 13
      title({'BS9'})  
    case 14
      title({'BS10'})        
    case 15
      title({'BS11'})  
    case 16
      title({'BS12'})
    case 17
      title({'BS5'})
    case 18
      title({'BS6'})  
    case 19
      title({'BS7'})        
    case 20
      title({'BS8'})  
    case 21
      title({'BS1'}) 
    case 22
      title({'BS2'})  
    case 23
      title({'BS3'})        
    case 24
      title({'BS4'})  
end
end
switch CellScenario
    case 'Urban'
        suptitle({'SINRij [dB] Urban Cell scenario'});
    case 'Rural'
        suptitle({'SINRij [dB] Rural Cell scenario'});
    case 'Micro'  
        suptitle({'SINRij [dB] Micro Cell scenario'});
end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
figure(301)
hold off
meshc(rx,ry,(SINRij(:,:,15)));
xlabel('X (km)');
ylabel('Y (km)');
zlabel('SINR [dB]');
minvle=min(min(SINRij(:,:,15)));
mxvlue=max(max(SINRij(:,:,15)));
%zticks([100 min(min(LidB(:,:,21))) vMCL(1) vMCL(2) vMCL(3) max(max(LidB(:,:,21)))])
switch CellScenario
    case 'Urban'
        title({'SINR BS15[dB] Urban Cell scenario'});
    case 'Rural'
        title({'SINR BS15 [dB] Rural Cell scenario'});
    case 'Micro'  
        title({'SINR BS15 [dB] Micro Cell scenario'});
end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

hold off
axes(handles.axes1);
title({'Narrow Band Internet of things','Attocell Network topology'});
plot(mtrxlxyz(:,X1),mtrxlxyz(:,Y1),'cs','MarkerFaceColor','b','MarkerSize',8);
text(mtrxlxyz(:,X1)+.32,mtrxlxyz(:,Y1)+.0,num2str((1:num_nod)'));% eNB
for e=1:num_ele
    line(mtrxlxyz(LaG(e,:),X1),mtrxlxyz(LaG(e,:),Y1),'Color',[.18 .8 .8]); 
end
xlim([-lx/2 lx/2])
ylim([-ly/2 ly/2])
grid on;
set(gca,'Ydir','reverse')
hold off
[mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2]=LntcyEnrgyPrfrmnc(n,u,tau);
%%---------------------------------------------------------------------
cntrsyst=[round(size(PRXijdB,2)/2) round(size(PRXijdB,1)/2)];%determine the center of the power matrix 
%Energy model
Eo=1000;%Initial Energy; Normal nodes: 5000 J, Advanced nodes: 10000 J
m=0.1;%Percentage of nodes than are advanced
a=1;%\alpha
centr=round(size(ttlLidB)/2);%center of network in path loss matrix
mxx2=centr(1);%
mxy2=centr(2);%
mxx1=max(mtrxlxyz(:,1));
mxy1=max(mtrxlxyz(:,2));
xfctr=(mxx2/mxx1);
yfctr=(mxy2/mxy1);
mtrxcc1v1=[];
mtrxcc1v2=[];
mtrxcc2v1=[];
mtrxcc2v2=[];
assignin('base','mtrxcc1v1',mtrxcc1v1);
assignin('base','mtrxcc1v2',mtrxcc1v2);
assignin('base','mtrxcc2v1',mtrxcc2v1);
assignin('base','mtrxcc2v2',mtrxcc2v2); 
%Creation of the random Sensor Network
axes(handles.axes1);
for i=1:1:n
    S(i).xd=(randi([xmin2*10,xmax2*10],1,1))/10;
    %S(i).xd=rand(1,1)*xmax;
    XR(i)=S(i).xd;
    S(i).yd=(randi([ymin2*10,ymax2*10],1,1))/10;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads (CH) only nodes
    S(i).type='N'; 
    temp_rnd0=i;  
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end   
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
        plot(S(i).xd,S(i).yd,'+');
        hold on;
    end      
    % pathloss
    [dif,pntxd]= min(abs(ry(:,1)-(S(i).yd)));%point x in the matrix
    [dif,pntyd]= min(abs(rx(1,:)-(S(i).xd)));%point y in the matrix
    S(i).pntxd=pntxd;
    S(i).pntyd=pntyd;    
  %  if  > Pthlsstargt   
    %--------------------------
    
    for j=1:eNB 
    distnctoBS(i,j)=d_eNB(pntxd,pntyd,j); % Distance from each node to the base station (BS)
    end
    
    dstncclsrBS=min(distnctoBS(i,:));%
    [dif,ClsrBS]= min(abs(distnctoBS(i,:)-dstncclsrBS));% 
    
    
    pthlsspnt=ttlLidB(pntxd,pntyd);
    
 
   
   % if pthlsspnt >= vMCL(1)
        if pthlsspnt <= vMCL(3)
           cvrgclss=2;      
        elseif pthlsspnt <= vMCL(2)         
           cvrgclss=1;  
        end
   % end
    
    
    vclsrBS(i,1)=ClsrBS;
    vclsrBS(i,2)=dstncclsrBS;% closer BS to the point (xd,yd)
    vclsrBS(i,3)=ttlLidB(pntxd,pntyd);% average pathloss from closer BS
    vclsrBS(i,4)=cvrgclss;% Coverage class 
    
    
    S(i).disttoBS=dstncclsrBS;
    S(i).BS=ClsrBS; 

end



%x and y Coordinates of the Sink
sink.x=0.5*xmax;
sink.y=0.5*ymax;

rmax=300;%maximum number of rounds

%variables %@@@@@@@
Ave_CH=0;%@@@@@@@
sum=0;%@@@@@@@
count_ch=0; %@@@@@@@
Throughput=0; %@@@@@@@
mth=0.2;%@@@@@@@@@@@
bth=1.5;%@@@@@@@
ath=2;
%@@@@@@

%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;

do=sqrt(Efs/Emp);%Computation of do


%Optimal Election Probability of a node
%to become cluster head
p=0.1;



       
%First Iteration

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;
Clust_throughput=[];
countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;


mtrxCH=zeros(1,1,n);

r=1;
countCHsth=0;
table1=zeros(6,4);
while r<=rmax

pause(0.1);
ckbx=get(handles.chckbx,'Value');
stp = evalin('base','stp'); 
% Number of rounds
    pause(0.1);
  

switch ckbx
    case 1%pause
    pause
    case 0
        if stp==0
           r=r+1;

        else
           r=rmax+1; 
           set(handles.pushbutton1,'Enable','on');
           set(handles.radiobutton3,'Enable','on');  
           set(handles.radiobutton4,'Enable','on');
           set(handles.radiobutton5,'Enable','on');
           quit cancel
        end
end
    

  %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

hold off;


axes(handles.axes1);

title({'Narrow Band Internet of things','Attocell Network topology'});
hold on
plot(mtrxlxyz(:,X1),mtrxlxyz(:,Y1),'cs','MarkerFaceColor','b','MarkerSize',10);

text(mtrxlxyz(:,X1)+.32,mtrxlxyz(:,Y1)+.0,num2str((1:num_nod)'));% eNB

for e=1:num_ele
    line(mtrxlxyz(LaG(e,:),X1),mtrxlxyz(LaG(e,:),Y1),'Color',[.18 .8 .8]); 
end

xlim([-lx/2 lx/2])
ylim([-ly/2 ly/2])
grid on;
set(gca,'Ydir','reverse')

hold off

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

%figure(2)
for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end

title({'Narrow Band Internet of things','Attocell Network topology'});
hold on
plot(mtrxlxyz(:,X1),mtrxlxyz(:,Y1),'bs','MarkerFaceColor','c','MarkerSize',10);
text(mtrxlxyz(:,X1)+.32,mtrxlxyz(:,Y1)+.0,num2str((1:num_nod)'));% eNB

for e=1:num_ele
    line(mtrxlxyz(LaG(e,:),X1),mtrxlxyz(LaG(e,:),Y1),'Color',[.1 .18 .8]); 
end
xlim([-lx/2 lx/2])
ylim([-ly/2 ly/2])
grid on;
set(gca,'Ydir','reverse')

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;


vBSn1=[];
vBSn2=[];
vBSn3=[];
vBSn4=[];
vBSn5=[];
vBSn6=[];
vBSn7=[];
vBSn8=[];
vBSn9=[];
vBSn10=[];
vBSn11=[];
vBSn12=[];
vBSn13=[];
vBSn14=[];
vBSn15=[];
vBSn16=[];
vBSn17=[];
vBSn18=[];
vBSn19=[];
vBSn20=[];
vBSn21=[];
vBSn22=[];
vBSn23=[];
vBSn24=[];



%create a matrix 
mxvtr=max(vclsrBS(:,2));
MtrxBSn=zeros();


%@@@@@@@@@@@@@@@@@@@@
sum=sum+cluster;
countCHsth=countCHsth+1

Ave_CH=(sum*0.1)/(1+(mth*ath)+(bth*x));
Throughput=Ave_CH*4*PACKETS_TO_BS(r);

STATISTICS(r).ave_clustHd=Ave_CH;
ave_ch(r)=Ave_CH;
STATISTICS(r).throughput=Throughput;
Clust_throughput(r)=Throughput;
if countCHsth==10

Ave_CH=0;
sum=0;
countCHsth=0;
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;
%@@@@@@@@@@@@@@@@@@@@

for i=1:1:n
   if(S(i).E>0)
       switch leach
           case 0%without using LEACH
               temp_rand=0; 
               trshldvle=1;% threshold value 
           case 1
               temp_rand=rand;
               trshldvle=(p/(1-p*mod(r,round(1/p))));% threshold value 
       end
  if leach ==1
  if ( (S(i).G)<=0)
       %Election of Cluster Heads
       
         if(temp_rand<= trshldvle)
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            Cls(cluster).xd=S(i).xd;
            Cls(cluster).yd=S(i).yd;
            %figure(2)
            hold on
            plot(S(i).xd,S(i).yd,'k*');
            d_Clstr=S(i).disttoBS;
            Cls(cluster).distance=d_Clstr;
            Cls(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
           
            vclstrlist(cluster,r)=i;
           % pause
            cluster=cluster+1;
            
            
            %% Queue for each base station BS

    
    
BSn=vclsrBS(i,1);

%----------------------

switch BSn
    case 1
     [vBSn1,S,table1,tocqueBS1]=queueBS(vBSn1,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 2
     [vBSn2,S,table1,tocqueBS2]=queueBS(vBSn2,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 3
     [vBSn3,S,table1,tocqueBS3]=queueBS(vBSn3,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 4
     [vBSn4,S,table1,tocqueBS4]=queueBS(vBSn4,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 5
     [vBSn5,S,table1,tocqueBS5]=queueBS(vBSn5,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 6
     [vBSn6,S,table1,tocqueBS6]=queueBS(vBSn6,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 7
     [vBSn7,S,table1,tocqueBS7]=queueBS(vBSn7,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 8
     [vBSn8,S,table1,tocqueBS8]=queueBS(vBSn8,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 9
     [vBSn9,S,table1,tocqueBS9]=queueBS(vBSn9,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 10
     [vBSn10,S,table1,tocqueBS10]=queueBS(vBSn10,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 11
     [vBSn11,S,table1,tocqueBS11]=queueBS(vBSn11,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 12
     [vBSn12,S,table1,tocqueBS12]=queueBS(vBSn12,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 13
     [vBSn13,S,table1,tocqueBS13]=queueBS(vBSn13,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 14
     [vBSn14,S,table1,tocqueBS14]=queueBS(vBSn14,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 15
     [vBSn15,S,table1,tocqueBS15]=queueBS(vBSn15,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 16
     [vBSn16,S,table1,tocqueBS16]=queueBS(vBSn16,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 17
     [vBSn17,S,table1,tocqueBS17]=queueBS(vBSn17,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 18
     [vBSn18,S,table1,tocqueBS18]=queueBS(vBSn18,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 19
     [vBSn19,S,table1,tocqueBS19]=queueBS(vBSn19,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 20
     [vBSn20,S,table1,tocqueBS20]=queueBS(vBSn20,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 21
     [vBSn21,S,table1,tocqueBS21]=queueBS(vBSn21,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 22
     [vBSn22,S,table1,tocqueBS22]=queueBS(vBSn22,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 23
     [vBSn23,S,table1,tocqueBS23]=queueBS(vBSn23,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 24
     [vBSn24,S,table1,tocqueBS24]=queueBS(vBSn24,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
end    
       %  Mtrxtoc(r,:)=[tocqueBS1 tocqueBS2 tocqueBS3 tocqueBS4 tocqueBS5 tocqueBS6 tocqueBS7 tocqueBS8 tocqueBS9 tocqueBS10 tocqueBS11 tocqueBS13 tocqueBS13 tocqueBS14 tocqueBS15 tocqueBS16 tocqueBS17tocqueBS18 tocqueBS19 tocqueBS20 tocqueBS21 tocqueBS22 tocqueBS23 tocqueBS24 ];    
   
        end  
    
  end
  else
      %----------------
      
          %Election of Cluster Heads
       
         if(temp_rand<= trshldvle)
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=round(1/p)-1;
            Cls(cluster).xd=S(i).xd;
            Cls(cluster).yd=S(i).yd;
            %figure(2)
            hold on
            plot(S(i).xd,S(i).yd,'k*');
            d_Clstr=S(i).disttoBS;
            Cls(cluster).distance=d_Clstr;
            Cls(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
           
            vclstrlist(cluster,r)=i;
           % pause
            cluster=cluster+1;
            
            
            %% Queue for each base station BS

    
    
BSn=vclsrBS(i,1);

%----------------------

switch BSn
    case 1
     [vBSn1,S,table1]=queueBS(vBSn1,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 2
     [vBSn2,S,table1]=queueBS(vBSn2,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 3
     [vBSn3,S,table1]=queueBS(vBSn3,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 4
     [vBSn4,S,table1]=queueBS(vBSn4,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 5
     [vBSn5,S,table1]=queueBS(vBSn5,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 6
     [vBSn6,S,table1]=queueBS(vBSn6,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 7
     [vBSn7,S,table1]=queueBS(vBSn7,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 8
     [vBSn8,S,table1]=queueBS(vBSn8,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 9
     [vBSn9,S,table1]=queueBS(vBSn9,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 10
     [vBSn10,S,table1]=queueBS(vBSn10,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 11
     [vBSn11,S,table1]=queueBS(vBSn11,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 12
     [vBSn12,S,table1]=queueBS(vBSn12,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 13
     [vBSn13,S,table1]=queueBS(vBSn13,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 14
     [vBSn14,S,table1]=queueBS(vBSn14,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 15
     [vBSn15,S,table1]=queueBS(vBSn15,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 16
     [vBSn16,S,table1]=queueBS(vBSn16,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 17
     [vBSn17,S,table1]=queueBS(vBSn17,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 18
     [vBSn18,S,table1]=queueBS(vBSn18,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 19
     [vBSn19,S,table1]=queueBS(vBSn19,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 20
     [vBSn20,S,table1]=queueBS(vBSn20,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 21
     [vBSn21,S,table1]=queueBS(vBSn21,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 22
     [vBSn22,S,table1]=queueBS(vBSn22,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 23
     [vBSn23,S,table1]=queueBS(vBSn23,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
    case 24
     [vBSn24,S,table1]=queueBS(vBSn24,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1);
end    
            
         
        end     
      
      
  end
   end 
  

end


STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;


vclstrlist(vclstrlist==0) = NaN;
if leach==1
for cluster=1:size(vclstrlist,1)
 if isnan(vclstrlist(cluster,r)) 
 else
   vcantCHr(r)=cluster;% total cluster head CH per round r
 end
 end
end

if leach==1
%Election of Associated Cluster Head CH for Normal Nodes
for i=1:1:n
     vclsrCHtoNi(i,1,r)=i;% Node id
   if ( S(i).type=='N' && S(i).E>0 )
     
      
       for CH=1:vcantCHr(r)
           % vdCHtoNN distance from each node to each cluster head CH each round r
           %vdCHtoNN(i,chckbx,r)=sqrt((S(i).xd-S(vclstrlist(chckbx)).xd)^2 + (S(i).yd-S(vclstrlist(chckbx)).yd)^2 );   
           
            vdCHtoNN(i,CH,r)= sqrt(((abs(S(i).xd-S(vclstrlist(CH,r)).xd)).^2)+((abs(S(i).yd-S(vclstrlist(CH,r)).yd)).^2));%+((lzi)^2));    
       end
           vdCHtoNN(vdCHtoNN==0) = NaN;
           clsrCHtoNi=min(vdCHtoNN(i,:,r));% 
           assignin('base','vdCHtoNN',vdCHtoNN);
           
           [dif,indclsrCHtoNi]= min(abs(vdCHtoNN(i,:,r)-clsrCHtoNi));% Index closest CH to node i
           
           
    % pathloss between users
    switch CellScenario% path loss model depends of the cell scenario [3], [29,30, 31, 32] %@@@@@@@
       case 'Urban'% Macro cell - Urban and Suburban areas
        LidBeNB=Grx + Gtx + PTXjdB + 2.7*10*log10(3.5/(4*pi*vdCHtoNN(i,indclsrCHtoNi,r))); % path loss %@@@@@@@
       case 'Rural'% Macro cell rural area
        LidBeNB=Grx + Gtx + PTXjdB + 2*10*log10(2/(4*pi*vdCHtoNN(i,indclsrCHtoNi,r)));  %@@@@@@@
       case 'Micro'% Micro dcell
        LidBeNB=Grx + Gtx + PTXjdB + 4*10*log10(6/(4*pi*vdCHtoNN(i,indclsrCHtoNi,r))); %@@@@@@@
    end

    %Coverage class
    if LidBeNB <= vMCL(2)
       cvrgclsseNB=1;      
    else        
       cvrgclsseNB=2;  
    end
           
           % Create matrix for cluster head CH assignament           
           vclsrCHtoNi(i,2,r)=vdCHtoNN(i,indclsrCHtoNi,r);% Distance closest CH to node i
           vclsrCHtoNi(i,3,r)=vclstrlist(indclsrCHtoNi,r);% Id closest CH to node i     
           vclsrCHtoNi(i,4,r)=LidBeNB;% Pathloss  
           vclsrCHtoNi(i,5,r)=cvrgclsseNB;% Coverage class  
           
           
           vclsrCHtoNi(vclsrCHtoNi==0) = NaN;           
           assignin('base','vclsrCHtoNi',vclsrCHtoNi);
           
           
           
   if(cluster-1>=1)
       min_dis=clsrCHtoNi;
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-Cls(c).xd)^2 + (S(i).yd-Cls(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       


     
       S(i).min_dis_cluster=min_dis_cluster;
           
   end
   end
 
if r>=1
if r <=rmax
mtrxenergy(r,i)=S(i).E;
end
end  
 
end

[mtrxCH,S]=queueCH(vclsrCHtoNi,u,tau,n,r,mtrxCH,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S);


assignin('base','mtrxCH',mtrxCH)

mtrxenergy(mtrxenergy==0) = NaN; 
assignin('base','mtrxenergy',mtrxenergy);
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;


BS1=evalin('base','BS1');
BS2=evalin('base','BS2');
BS3=evalin('base','BS3');
BS4=evalin('base','BS4');            
BS5=evalin('base','BS5');
BS6=evalin('base','BS6');
BS7=evalin('base','BS7');
BS8=evalin('base','BS8');             
BS9=evalin('base','BS9');
BS10=evalin('base','BS10');
BS11=evalin('base','BS11');
BS12=evalin('base','BS12');            
BS13=evalin('base','BS13');
BS14=evalin('base','BS14');
BS15=evalin('base','BS15');
BS16=evalin('base','BS16');           % 
BS17=evalin('base','BS17');
BS18=evalin('base','BS18');
BS19=evalin('base','BS19');
BS20=evalin('base','BS20');            
BS21=evalin('base','BS21');
BS22=evalin('base','BS22');
BS23=evalin('base','BS23');
BS24=evalin('base','BS24');           
            


assignin('base','vBSn1',vBSn1);
assignin('base','vBSn1',vBSn2);
assignin('base','vBSn1',vBSn3);
assignin('base','vBSn1',vBSn3);
assignin('base','vBSn1',vBSn4);
assignin('base','vBSn1',vBSn6);
assignin('base','vBSn1',vBSn7);
assignin('base','vBSn1',vBSn8);
assignin('base','vBSn1',vBSn9);
assignin('base','vBSn1',vBSn10);
assignin('base','vBSn1',vBSn11);
assignin('base','vBSn1',vBSn12);
assignin('base','vBSn1',vBSn13);
assignin('base','vBSn1',vBSn14);
assignin('base','vBSn1',vBSn15);
assignin('base','vBSn1',vBSn16);
assignin('base','vBSn1',vBSn17);
assignin('base','vBSn1',vBSn18);
assignin('base','vBSn1',vBSn19);
assignin('base','vBSn1',vBSn20);
assignin('base','vBSn1',vBSn21);
assignin('base','vBSn1',vBSn22);
assignin('base','vBSn1',vBSn23);
assignin('base','vBSn1',vBSn24);

clc
fprintf('     BS       d(km)   pathloss     class     t(s)     td(s)    eNB   EdsptUJ   EdsptDJ')
vBSn1
vBSn2
vBSn3
vBSn4
vBSn5
vBSn6
vBSn7
vBSn8
vBSn9
vBSn10
vBSn11
vBSn12
vBSn13
vBSn14
vBSn15
vBSn16
vBSn17
vBSn18
vBSn19
vBSn20
vBSn21
vBSn22
vBSn23
vBSn24

if size(vBSn24)>0
figure(125)
hold on
vBSn24
size(vBSn24)
vtrtoplt=vBSn24(:,10)
length(vtrtoplt)

vtrtopltttl=0;
for it=1:length(vtrtoplt)
vtrtopltttl=vtrtopltttl+vtrtoplt(it)
end

    

plot(r,vtrtopltttl,'*r')
title('time at the BS#24')
xlabel('Number of rounds r')
ylabel('time')
end
%% 
assignin('base','mtrxcc1v1',mtrxcc1v1);
assignin('base','mtrxcc1v2',mtrxcc1v2);
assignin('base','mtrxcc2v1',mtrxcc2v1);
assignin('base','mtrxcc2v2',mtrxcc2v2);  
assignin('base','S',S);
Clust_throughput
%@@@@@@@@@@@@@@@@@plot Throughput
if size(Clust_throughput)>0
figure(122)
Clust_throughput
plot(Clust_throughput)
title('data rate from CHs to BS ')
xlabel('Number of rounds r')
ylabel('Throughput in Kbits')
hold off
end
figure(123)
hold on
plot(r,PACKETS_TO_BS(r),'*b')
title('Number of packets to BS')
xlabel('Number of rounds r')
ylabel('Number of packets from CHs to BS')
%@@@@@@@@@@@@@@@@@

nodesperround(r)=n-dead-dead_a-dead_n;
figure(124)
hold on
plot(r,nodesperround(r),'.r')
title('Number of Live Nodes per round')
xlabel('Number of rounds r')
ylabel('Number of Live Nodes per round')

%Send data to GUI  
set(handles.text2,'String',num2str(countCHs));
if r < rmax
set(handles.text2,'String',num2str(r));
end
pause(0.1);
end

% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in chckbx.
function chckbx_Callback(hObject, eventdata, handles)

% hObject    handle to chckbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chckbx


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
psckbx=get(handles.chckbx,'Value');
set(handles.text3,'String',num2str(psckbx));


% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
stp=1;
assignin('base','stp',stp)
% plot path loss



% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
BS21=evalin('base','BS21');
switch BS21
    case 0
        BS21=1;
        set(handles.pushbutton36,'BackgroundColor','blue');
    case 1
        BS21=0;
        set(handles.pushbutton36,'BackgroundColor','cyan');  
end
 assignin('base','BS21',BS21)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton37.
function pushbutton37_Callback(hObject, eventdata, handles)
BS22=evalin('base','BS22');
switch BS22
    case 0
        BS22=1;
        set(handles.pushbutton37,'BackgroundColor','blue');
    case 1
        BS22=0;
        set(handles.pushbutton37,'BackgroundColor','cyan');  
end
 assignin('base','BS22',BS22)
% hObject    handle to pushbutton37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton38.
function pushbutton38_Callback(hObject, eventdata, handles)
BS23=evalin('base','BS23');
switch BS23
    case 0
        BS23=1;
        set(handles.pushbutton38,'BackgroundColor','blue');
    case 1
        BS23=0;
        set(handles.pushbutton38,'BackgroundColor','cyan');  
end
 assignin('base','BS23',BS23)
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton39.
function pushbutton39_Callback(hObject, eventdata, handles)
BS24=evalin('base','BS24');
switch BS24
    case 0
        BS24=1;
        set(handles.pushbutton39,'BackgroundColor','blue');
    case 1
        BS24=0;
        set(handles.pushbutton39,'BackgroundColor','cyan');  
end
 assignin('base','BS24',BS24)
% hObject    handle to pushbutton39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton40.
function pushbutton40_Callback(hObject, eventdata, handles)
BS17=evalin('base','BS17');
switch BS17
    case 0
        BS17=1;
        set(handles.pushbutton40,'BackgroundColor','blue');
    case 1
        BS17=0;
        set(handles.pushbutton40,'BackgroundColor','cyan');  
end
 assignin('base','BS17',BS17)
% hObject    handle to pushbu
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
BS18=evalin('base','BS18');
switch BS18
    case 0
        BS18=1;
        set(handles.pushbutton41,'BackgroundColor','blue');
    case 1
        BS18=0;
        set(handles.pushbutton41,'BackgroundColor','cyan');  
end
 assignin('base','BS18',BS18)
% hObject    handle to pushbu
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
BS19=evalin('base','BS19');
switch BS19
    case 0
        BS19=1;
        set(handles.pushbutton42,'BackgroundColor','blue');
    case 1
        BS19=0;
        set(handles.pushbutton42,'BackgroundColor','cyan');  
end
 assignin('base','BS19',BS19)
% hObject    handle to pushbu
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
BS20=evalin('base','BS20');
switch BS20
    case 0
        BS20=1;
        set(handles.pushbutton43,'BackgroundColor','blue');
    case 1
        BS20=0;
        set(handles.pushbutton43,'BackgroundColor','cyan');  
end
 assignin('base','BS20',BS20)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
BS13=evalin('base','BS13');
switch BS13
    case 0
        BS13=1;
        set(handles.pushbutton44,'BackgroundColor','blue');
    case 1
        BS13=0;
        set(handles.pushbutton44,'BackgroundColor','cyan');  
end
 assignin('base','BS13',BS13)
% hObject    handle to pushbutton44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
BS14=evalin('base','BS14');
switch BS14
    case 0
        BS14=1;
        set(handles.pushbutton45,'BackgroundColor','blue');
    case 1
        BS14=0;
        set(handles.pushbutton45,'BackgroundColor','cyan');  
end
 assignin('base','BS14',BS14)
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
BS15=evalin('base','BS15');
switch BS15
    case 0
        BS15=1;
        set(handles.pushbutton46,'BackgroundColor','blue');
    case 1
        BS15=0;
        set(handles.pushbutton46,'BackgroundColor','cyan');  
end
 assignin('base','BS15',BS15)
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
BS16=evalin('base','BS16');
switch BS16
    case 0
        BS16=1;
        set(handles.pushbutton47,'BackgroundColor','blue');
    case 1
        BS16=0;
        set(handles.pushbutton47,'BackgroundColor','cyan');  
end
 assignin('base','BS16',BS16)
% hObject    handle to pushbu
% hObject    handle to pushbutton47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
BS9=evalin('base','BS9');
switch BS9
    case 0
        BS9=1;
        set(handles.pushbutton48,'BackgroundColor','blue');
    case 1
        BS9=0;
        set(handles.pushbutton48,'BackgroundColor','cyan');  
end
 assignin('base','BS9',BS9)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
BS10=evalin('base','BS10');
switch BS10
    case 0
        BS10=1;
        set(handles.pushbutton49,'BackgroundColor','blue');
    case 1
        BS10=0;
        set(handles.pushbutton49,'BackgroundColor','cyan');  
end
 assignin('base','BS10',BS10)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, ~, handles)
BS11=evalin('base','BS11');
switch BS11
    case 0
        BS11=1;
        set(handles.pushbutton50,'BackgroundColor','blue');
    case 1
        BS11=0;
        set(handles.pushbutton50,'BackgroundColor','cyan');  
end
 assignin('base','BS11',BS11)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton51.
function pushbutton51_Callback(hObject, eventdata, handles)
BS12=evalin('base','BS12');
switch BS12
    case 0
        BS12=1;
        set(handles.pushbutton51,'BackgroundColor','blue');
    case 1
        BS12=0;
        set(handles.pushbutton51,'BackgroundColor','cyan');  
end
 assignin('base','BS12',BS12)
% hObject    handle to pushbutton51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton84.
function pushbutton84_Callback(hObject, eventdata, handles)
BS5=evalin('base','BS5');
switch BS5
    case 0
        BS5=1;
        set(handles.pushbutton84,'BackgroundColor','blue');
    case 1
        BS5=0;
        set(handles.pushbutton84,'BackgroundColor','cyan');  
end
 assignin('base','BS5',BS5)

% hObject    handle to pushbutton84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton85.
function pushbutton85_Callback(hObject, eventdata, handles)
BS6=evalin('base','BS6');
switch BS6
    case 0
        BS6=1;
        set(handles.pushbutton85,'BackgroundColor','blue');
    case 1
        BS6=0;
        set(handles.pushbutton85,'BackgroundColor','cyan');  
end
 assignin('base','BS6',BS6)
% hObject    handle to pushbutton85 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton86.
function pushbutton86_Callback(hObject, eventdata, handles)
BS7=evalin('base','BS7');
switch BS7
    case 0
        BS7=1;
        set(handles.pushbutton86,'BackgroundColor','blue');
    case 1
        BS7=0;
        set(handles.pushbutton86,'BackgroundColor','cyan');  
end
 assignin('base','BS7',BS7)
% hObject    handle to pushbutton86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton87.
function pushbutton87_Callback(hObject, eventdata, handles)
BS8=evalin('base','BS8');
switch BS8
    case 0
        BS8=1;
        set(handles.pushbutton87,'BackgroundColor','blue');
    case 1
        BS8=0;
        set(handles.pushbutton87,'BackgroundColor','cyan');  
end
 assignin('base','BS8',BS8)
% hObject    handle to pushbutton87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton88.
function pushbutton88_Callback(hObject, eventdata, handles)
BS1=evalin('base','BS1');
switch BS1
    case 0
        BS1=1;
        set(handles.pushbutton88,'BackgroundColor','blue');
    case 1
        BS1=0;
        set(handles.pushbutton88,'BackgroundColor','cyan');  
end
 assignin('base','BS1',BS1)
% hObject    handle to pushbutton88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton89.
function pushbutton89_Callback(hObject, eventdata, handles)
BS2=evalin('base','BS2');
switch BS2
    case 0
        BS2=1;
        set(handles.pushbutton89,'BackgroundColor','blue');
    case 1
        BS2=0;
        set(handles.pushbutton89,'BackgroundColor','cyan');  
end
 assignin('base','BS2',BS2)
% hObject    handle to pushbutton89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton90.
function pushbutton90_Callback(hObject, eventdata, handles)
BS3=evalin('base','BS3');
switch BS3
    case 0
        BS3=1;
        set(handles.pushbutton90,'BackgroundColor','blue');
    case 1
        BS3=0;
        set(handles.pushbutton90,'BackgroundColor','cyan');  
end
 assignin('base','BS3',BS3)
% hObject    handle to pushbutton90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton91.
function pushbutton91_Callback(hObject, eventdata, handles)
BS4=evalin('base','BS4');
switch BS4
    case 0
        BS4=1;
        set(handles.pushbutton91,'BackgroundColor','blue');
    case 1
        BS4=0;
        set(handles.pushbutton91,'BackgroundColor','cyan');  
end
 assignin('base','BS4',BS4)

% hObject    handle to pushbutton91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function pushbutton36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
