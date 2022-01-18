%% Create a funtion for the solver
% Caitlin Sifuentes
function [unew,vnew,hnew] = LW_CS(dx,dy,dt,g,u,v,h,Kx,Ky)
uh = u.*h;
vh = v.*h;
h_onehalf_x = (h(2:end,:)+h(1:end-1,:))/2 - (dt/dx)*(.5)*(uh(2:end,:)-uh(1:end-1,:));
h_onehalf_y = (h(:,2:end)+h(:,1:end-1))/2 - (dt/dy)*(.5)*(vh(:,2:end)-vh(:,1:end-1));

% x-direction
Fxx = u.*uh + (0.5)*g.*h.^2; % X-direction first %NOTE: if .5 .* g doesn't work, remove dots from here
Fyx = v.*uh; %+ (0.5)*g.*h.^2; 

uh_onehalf_x = (0.5)*(uh(2:end,:)+uh(1:end-1,:)) - (dt/dx)*(.5)*(Fxx(2:end,:)-Fxx(1:end-1,:));
uh_onehalf_y = (0.5)*(uh(:,2:end)+uh(:,1:end-1)) - (dt/dy)*(.5)*(Fyx(:,2:end)-Fyx(:,1:end-1));

% y-dir
Fyy = v.*vh + (0.5)*g*h.^2;
Fxy = v.*uh; 

vh_onehalf_x = (0.5)*(vh(2:end,:)+vh(1:end-1,:)) - (dt/dx)*(.5).*(Fxy(2:end,:)-Fxy(1:end-1,:));
vh_onehalf_y = (0.5)*(vh(:,2:end)+vh(:,1:end-1)) - (dt/dy)*(.5).*(Fyy(:,2:end)-Fyy(:,1:end-1));

hnew = h(2:end-1,2:end-1) - (dt/dx)*(uh_onehalf_x(2:end,2:end-1)-uh_onehalf_x(1:end-1,2:end-1))-(dt/dy).*(vh_onehalf_y(2:end-1,2:end)-vh_onehalf_y(2:end-1,1:end-1));  % We're going to define boundaries, therefore we index everywhere that isn't these boundaries

Fxx_onehalf_x = uh_onehalf_x.*(uh_onehalf_x./h_onehalf_x) + (0.5)*g.*h_onehalf_x.^2;
Fyx_onehalf_x = uh_onehalf_y.*vh_onehalf_y./h_onehalf_y;
Fyy_onehalf_y = vh_onehalf_y.*(vh_onehalf_y./h_onehalf_y) + (0.5)*g.*h_onehalf_y.^2;
Fxy_onehalf_y = vh_onehalf_x.*uh_onehalf_x./h_onehalf_x;

uhnew = uh(2:end-1,2:end-1) - (dt/dx).*(Fxx_onehalf_x(2:end,2:end-1)-Fxx_onehalf_x(1:end-1,2:end-1))-(dt/dy).*(Fyx_onehalf_x(2:end-1,2:end)-Fyx_onehalf_x(2:end-1,1:end-1)) - dt.*Kx; % Minus because Kx is momentum sink
vhnew = vh(2:end-1,2:end-1) - (dt/dx).*(Fxy_onehalf_y(2:end,2:end-1)-Fxy_onehalf_y(1:end-1,2:end-1))-(dt/dy).*(Fyy_onehalf_y(2:end-1,2:end)-Fyy_onehalf_y(2:end-1,1:end-1)) - dt.*Ky;

unew = uhnew./hnew;
vnew = vhnew./hnew;

end


