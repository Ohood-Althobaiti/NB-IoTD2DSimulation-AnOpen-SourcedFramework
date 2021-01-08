        
function [vBSnx,S,table1,tocqueBS]=queueBS(vBSnx,vclsrBS,i,u,tau,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S,table1)
tic;
BNpo=size(vBSnx,1)+1;
  for auxq=1:size(vclsrBS,2)
    vBSnx(BNpo,auxq)=vclsrBS(i,auxq);
    t=BNpo*u;
    td=BNpo*tau;
    vBSnx(BNpo,auxq+1)=t;%  
    vBSnx(BNpo,auxq+2)=td;%
    [diftt,pnttt]= min(abs(mtrxtt(1,:)-t));%point in the matrix  
    [difdt,pntdt]= min(abs(mtrxdt(:,1)-td));%point in the matrix 
    vBSnx(BNpo,auxq+3)=pnttt;%
    vBSnx(BNpo,auxq+4)=pntdt;%  
    class=vBSnx(BNpo,4);

      switch class
        case 1
            EdsptUJ=mEu1(pntdt,pnttt)*BNpo;
            EdsptDJ=mEd1(pntdt,pnttt)*BNpo;
        otherwise
            EdsptUJ=mEu2(pntdt,pnttt)*BNpo;
            EdsptDJ=mEd2(pntdt,pnttt)*BNpo;       
      end   
    vBSnx(BNpo,auxq+5)=EdsptUJ;% Devices' energy consumpmtion for an Uplink request
    vBSnx(BNpo,auxq+6)=EdsptDJ;% Devices' energy consumpmtion for a Downlink request 
    ttlEnrgy=EdsptUJ+EdsptDJ;
    vBSnx(BNpo,auxq+7)=ttlEnrgy;% Devices' energy consumpmtion for a Uplink/Downlink request
    
    S(i).E=S(i).E-ttlEnrgy;
    vBSnx(BNpo,auxq+10)=t;
    vBSnx(BNpo,auxq+11)=td;
  end
vBSnx(BNpo,auxq+8)=i;
tocqueBS=toc;
vBSnx(BNpo,auxq+9)=tocqueBS;

end