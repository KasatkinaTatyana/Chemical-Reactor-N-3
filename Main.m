close all
clear all
clc
global x_0 x_End
global y_0 y_End
global t_0
global R E1 E2 A10 A20
global k b alpha betta
%--------------------------------------------------------------------------
%-------------------------Константы реактора-------------------------------
R=8.3;
E1=2.09e4;
E2=4.18e4;
A10=1.1;
A20=172.2; 
%----------------------Коэффициенты управления-----------------------------
a0=4.3145;
a1=-0.1099;
b0=1.4962;
b1=0.0515;
g1=41.8;
g2=83.6;
%--------------------------------------------------------------------------
t_0=0;
h_t=0.05;       %шаг по времени t
c_st=[1 3 3];
%--------------------------------------------------------------------------
global N    %число разбиений
N=1000;

x_0=[0.993 0.007 298];
x_End=[0.35 0.55 333];
%--------------------------------------------------------------------------
%Граничные значения канонических переменных
y_0=YTrans(x_0);
y_End=YTrans(x_End);
%--------------------------------------------------------------------------
k=(y_End(2)-y_0(2))/(y_End(1)-y_0(1));
b=(y_0(2)*y_End(1)-y_End(2)*y_0(1))/(y_End(1)-y_0(1));
alpha = (y_0(3)-k*y_0(2))/y_0(2)/(y_0(1)-y_End(1));
betta = (y_0(2)*y_End(3)+y_0(3)*y_End(2)-2*k*y_0(2)*y_End(2))/y_0(2)/y_End(2)/(y_End(1)-y_0(1))^2;
%--------------------------------------------------------------------------
% [cUp cD]=CAnalyses();
% %----------------------------Test CAnalyses--------------------------------
% N_y=1000;
% dy=(y_End(1)-y_0(1))/N_y;
% YUp=zeros(1,N_y+1);
% dYUp=zeros(1,N_y+1);
% YD=zeros(1,N_y+1);
% dYD=zeros(1,N_y+1);
% j=1;
% for y=y_0(1):dy:y_End(1)
%     YUp(j)=k*y+b+cUp*(y-y_0(1))*(y_End(1)-y);
%     dYUp(j)=(k+cUp*(y_0(1)+y_End(1)-2*y))*YUp(j);
%     YD(j)=k*y+b+cD*(y-y_0(1))*(y_End(1)-y);
%     dYD(j)=(k+cD*(y_0(1)+y_End(1)-2*y))*YD(j);
%     j=j+1;
% end
% disp('max(YUp)=');
% disp(max(YUp));
% disp('max(YD)=');
% disp(max(YD));
% disp('min(dYUp-YUp)=');
% disp(min(dYUp-YUp));
% disp('min(dYD-YD)=');
% disp(min(dYD-YD));


dy=(y_End(1)-y_0(1))/N;

tau_End=quad(@yCube,y_0(1),y_End(1));
dtau=tau_End/N;

%% Определение t_End
t_End=t_0;
for y=y_0(1):dy:y_End(1)
    y0=y;
    y1=k*y0+b+alpha*(y0-y_0(1))*(y0-y_End(1))+betta*(y0-y_0(1))^2*(y0-y_End(1));

    y2=(k+alpha*(y0-y_0(1))+alpha*(y0-y_End(1))+2*betta*(y0-y_0(1))*(y0-y_End(1))+betta*(y0-y_0(1))^2)*y1;

    y3=2*(alpha+3*betta*y0-2*betta*y_0(1)-betta*y_End(1))*y1^2+...
    (k+alpha*(y0-y_0(1))+alpha*(y0-y_End(1))+2*betta*(y0-y_0(1))*(y0-y_End(1))+betta*(y0-y_0(1))^2)^2*y1;

    V_y=[y0 y1 y2];
    
    T=InvStateDiff(V_y);
    t_End=t_End-dy/(k1(T)*y^2);
end
t_End
%% Тест правильности закона управления
% i=1;
% y_t=y_0;
% 
% for tau=0:dtau:tau_End
%     U=Control(tau);
%     u=U(1,1);
% %     y_t(i+1,1)=y_t(i,1)+dtau*U(5,1);
% %     y_t(i+1,2)=y_t(i,2)+dtau*U(6,1);
% %     y_t(i+1,3)=y_t(i,3)+dtau*(fKan(U(4:6,1))+gKan(U(4:6,1))*u);
%     y_t(i+1,1)=y_t(i,1)+dtau*y_t(i,2);
%     y_t(i+1,2)=y_t(i,2)+dtau*y_t(i,3);
%     y_t(i+1,3)=y_t(i,3)+dtau*(fKan(y_t(i,:))+gKan(y_t(i,:))*u);
%     i=i+1;  
% end
% 
% disp(y_t(end,:));
% disp(U(4:6,1));
%% Тест моделирования системы неканонического вида, дифференцирование по tau
% x_t=x_0;
% i=1;
% for tau=0:dtau:tau_End
%     U=Control(tau);
%     u=U(1,1);
%     x_t(i+1,1)=x_t(i,1)+dtau*(-x_t(i,2));
%     x_t(i+1,2)=x_t(i,2)+dtau*(x_t(i,2)-k2(x_t(i,3))/k1(x_t(i,3)) * x_t(i,2)^2/x_t(i,1)^2 );
%     x_t(i+1,3)=x_t(i,3)+dtau*u;
%     i=i+1;
% end
% 
% disp(x_t(end,:));
% disp(U(4:6,1));
%% Моделирование исходной системы в реальном времени t
x_r=x_0; tau=t_0;
i=1;
for t=t_0:h_t:t_End
    y_r=YTrans(x_r(i,:));
    
    U=Control(tau);
    v=U(1,1);
    f_tau=U(2,1);
    g_tau=U(3,1);
    y_tau=U(4:6,1);
    
    f_r=fKan(y_r);
    g_r=gKan(y_r);
    
    v=(-f_r + U(7,1) - 0.3*(y_r(3) - y_tau(3)) - 0.03*(y_r(2) - y_tau(2)) - 0.001*(y_r(1) - y_tau(1)) )/g_r;
    
    T=InvStateDiff(y_tau);
    h_tau = h_t / s([y_tau(1) -y_tau(2) T]);
    
    x_r(i+1,1)=x_r(i,1)+(-k1(x_r(i,3))*x_r(i,1)^2)*h_t;
    x_r(i+1,2)=x_r(i,2)+( k1(x_r(i,3))*x_r(i,1)^2 - k2(x_r(i,3))*x_r(i,2) )*h_t;
    x_r(i+1,3)=x_r(i,3)+ v / s(x_r(i,:)) * h_t;
    i=i+1;
    tau=tau+h_tau;
end

disp(x_r(end,:));
disp(U(4:6,1));
