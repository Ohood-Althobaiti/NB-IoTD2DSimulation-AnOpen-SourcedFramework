
function [mtrxCH,S]=queueCH(vclsrCHtoNi,u,tau,n,r,mtrxCH,mtrxtt,mtrxdt,mEu1,mEu2,mEd1,mEd2,S)
%u=0.01;
%tau=0.002;
%n=300;
%r=2;
%mtrxCH=zeros(1,1,n);
aux2=1;
aux3=1;
for i=1:n
cluster=vclsrCHtoNi(i,1,r);
d=vclsrCHtoNi(i,2,r);
eNB=vclsrCHtoNi(i,3,r);%associated NCH       
class=vclsrCHtoNi(i,5,r);
aux1=1;

    if isnan(vclsrCHtoNi(i,2,r))% it device is a CH   
        mtrxCH(:,:,cluster)=NaN; 
       else%find the CH and add to it's queue       
       dvcsinqn=size(mtrxCH(:,:,eNB),1);    
       mtrxCH(mtrxCH==0) = NaN;           
       while aux1<=dvcsinqn         
         if isnan(mtrxCH(aux1,1,eNB))%there are not a device assignated
            mtrxCH(aux3,1,eNB)=i;
            mtrxCH(aux3,2,eNB)=d;
            mtrxCH(aux3,3,eNB)=vclsrCHtoNi(i,4,r);
            mtrxCH(aux3,4,eNB)=class;
            t=aux1*u;
            td=aux1*tau;
            [diftt,pnttt]= min(abs(mtrxtt(1,:)-t));%point in the matrix  
            [difdt,pntdt]= min(abs(mtrxdt(:,1)-td));%point in the matrix     
            switch class
              case 1
                  if aux1==1
                   EdsptUJ=mEu1(pntdt,pnttt);
                   EdsptDJ=mEd1(pntdt,pnttt);
                  else 
                   EdsptUJ=mEu1(pntdt,pnttt);%+sum(mtrxCH(:,4,eNB));
                   EdsptDJ=mEd1(pntdt,pnttt);%+sum(mtrxCH(:,5,eNB));   
                  end
                otherwise
                  if aux1==1
                   EdsptUJ=mEu2(pntdt,pnttt);
                   EdsptDJ=mEd2(pntdt,pnttt); 
                  else
                   EdsptUJ=mEu2(pntdt,pnttt);%+sum(mtrxCH(:,4,eNB));
                   EdsptDJ=mEd2(pntdt,pnttt);%+sum(mtrxCH(:,5,eNB));  
                      
                  end
            end  
              mtrxCH(aux3,4,eNB)=class; %@@@@@@@@@@@
              mtrxCH(aux3,5,eNB)=EdsptUJ;%@@@@@@@@@@
              mtrxCH(aux3,6,eNB)=EdsptDJ;%@@@@@@@@@@@
              ttlEnrgy=EdsptUJ+EdsptDJ;
              mtrxCH(aux3,7,eNB)=ttlEnrgy; %@@@@@@@@@@@@@
              S(i).E=S(i).E-ttlEnrgy;
              aux3=1;
              aux1=dvcsinqn+1; 
          else          
             aux1=aux1+1;
             aux2=aux2+1;
             aux3=aux3+1;
         end              
       end     
       if isnan(d)% it device is a CH     
        mtrxCH(:,:,cluster)=NaN; 
       end      
    end
end
end