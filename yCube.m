function F=yCube(y)
global k b alpha betta y_0 y_End
%  F=1./(k*y+b+c*(y-y_0(1)).*(y_End(1)-y));
F=1./(k*y+b+alpha*(y-y_0(1)).*(y-y_End(1))+betta*((y-y_0(1)).^2).*(y-y_End(1)));