function [cUp cD]=CAnalyses()
global c k b y_0 y_End 
c=0;
hc=0.01;
N_y=1000;
N_stop=10000;


dy=(y_End(1)-y_0(1))/N_y;
i=0;
Y=zeros(1,N_y+1);
dY=zeros(1,N_y+1);
while i<N_stop
    j=1;
    for y=y_0(1):dy:y_End(1)
        Y(j)=k*y+b+c*(y-y_0(1))*(y_End(1)-y);
        dY(j)=(k+c*(y_0(1)+y_End(1)-2*y))*Y(j);
        j=j+1;
    end
%     if (max(Y)>=0)||(max(dY)>=1)
    if (max(Y)>=0)||(min(dY-Y)<=0.1)
        break
    end
    c=c+hc;
    i=i+1;
end
cUp=c-hc;
    
c=0;
i=0;
while i<N_stop
    j=1;
    for y=y_0(1):dy:y_End(1)
        Y(j)=k*y+b+c*(y-y_0(1))*(y_End(1)-y);
        dY(j)=(k+c*(y_0(1)+y_End(1)-2*y))*Y(j);
        j=j+1;
    end
%     if (max(Y)>=0)||(max(dY)>=1)
    if (max(Y)>=0)||(min(dY-Y)<=0.1)
        break
    end
    c=c-hc;
    i=i+1;
end
cD=c+hc;