function F=Euler(y_curr, t_curr, t_new, h_t)
% y_curr - ������� �������� y (y0)
% t_curr - ������� �������� tau
% t_new - ����� �������� tau, � ������� ����� �������� y
% h_t - ��� ��������������
global y_0 y_End k b alpha betta
F=y_curr;

if (t_new > t_curr)
    y=y_curr;
    for dt=t_curr : h_t : t_new
        y=y+h_t*(k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1)));
    end
    F=y;
end
