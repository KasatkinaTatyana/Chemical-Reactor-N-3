function F=Euler(t,h_t)
global y_0 y_End k b alpha betta
F=y_0(1);

if (t>0)
    y=y_0(1);
    for dt=0:h_t:t
        y=y+h_t*(k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1)));
        F=y;
    end
end
