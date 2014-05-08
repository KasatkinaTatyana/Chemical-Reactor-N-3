function dy=Traject(t,y)
global cUp cD k b y_0 y_End lambda
dy=zeros(1,1);
yD=k*y+b+cD*(y-y_0(1))*(y_End(1)-y);
yUp=k*y+b+cUp*(y-y_0(1))*(y_End(1)-y);
dy=(lambda/yD+(1-lambda)/yUp)^(-1);