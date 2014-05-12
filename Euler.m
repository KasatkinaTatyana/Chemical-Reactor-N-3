function F=Euler(y_curr, t_curr, t_new, h_t)
% y_curr - текущее значение y (y0)
% t_curr - текущее значение tau
% t_new - новое значение tau, в котором будет вычислен y
% h_t - шаг интегрирования
global y_0 y_End k b alpha betta
F=y_curr;

if (t_new > t_curr)
    y=y_curr;
    for dt=t_curr : h_t : t_new
        y=y+h_t*(k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1)));
    end
    F=y;
end
