%Shallow Water Model
clear
close all

%define constants
g=9.81;
nx = 250;
ny = 50;
dx=0.1; %m
dy=0.1; %m
total_time_days=1;
dt_mins=0.00001;
fover8=0.01/8; %friction coeff

%Initialize
dt= dt_mins*60; %time in seconds
total_time_seconds = total_time_days*24*3600; %time length (s)
nt = fix(total_time_seconds/dt)+1; %number of timesteps

%coordinates
x=(0:nx-1).*dx;
y=(0:ny-1).*dy;
[Y,X]=meshgrid(y,x);

H=zeros(nx,ny)-5; %bottom surface elevation 5m deep everywhere

% % elevchange=900;
% % dx=elevchange/(nx-1);
% % bedelev2=0:dx:elevchange;
% % H=(bedelev2'*ones(1,ny));

std_bump=1*dx;
amp = 5;
height=amp.*exp(-((X-mean(x)).^2+(Y-mean(y)).^2)./(2*std_bump^2));
%contourf(X,Y,height)

%Fluid is initially at rest
u=zeros(nx, ny);
v=zeros(nx, ny);

%Define h as the depth if the fluid (whereas "height" is the height of the
%upper surface)

h = height -H;
h_init = height - H; %h plus H = water surface

save_interval=5; %save output every # timesteps
save_n=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop

for n = 1:500
    
    %Compute Kx, Ky -> losses from friction
    Sx = -((H(3:end,2:end-1)-H(1:end-2,2:end-1)))/(2*dx);
    Sy = -((H(2:end-1,3:end)-H(2:end-1,1:end-2)))/(2*dy);
    Kx = (fover8.*u(2:end-1,2:end-1).*abs(u(2:end-1,2:end-1)))-g.*Sx.*h(2:end-1,2:end-1);
    Ky = (fover8.*v(2:end-1,2:end-1).*abs(v(2:end-1,2:end-1)))-g.*Sy.*h(2:end-1,2:end-1);
    
    %Call the Lax-Wendroff scheme to move forward one timestep
    [unew, vnew, hnew] = Lax_Wendroff_SolutionScheme(dx, dy, dt, g, u, v, h, Kx, Ky);
    %Update and enforce boundary conditions
    u=unew([end 1:end 1],[1 1:end end]); %periodic in x
    v=vnew([1 1:end end],[1 1:end end]); %transmissive in y
    
    v(:,[1 end]) = 0; %no flow in y-direction along boundary
    h(:,2:end-1) = hnew([end 1:end 1],:);
%     if mod(n,25)==0
%     pause(0.01)
%     %contourf(X,Y,h)
%     imagesc(h')
%     title('h')
%     colorbar
%     end
end

max(max(h))
min(min(h))
