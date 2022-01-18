%% HW4 
clear 
close all

tic
%% Load in Data
U = xlsread('LakeErieCORRECT.xlsx','U-dir Vel (m s-1)');
V = xlsread('LakeErieCORRECT.xlsx','V-dir Vel (m s-1)');
Depth = xlsread('LakeErieCORRECT.xlsx','Depth (m)');

%% Visualize
figure
quiver(V,U,10) % plot vector map for Lake Erie to show movement of flow
figure
imagesc(Depth) % plot depth
set(gca,'YDir','normal')

%% List Knowns
% Concentration Matrix
jmax = length(Depth(:,1));
kmax = length(Depth(1,:));
sediment = zeros(jmax,kmax); % initialize concentration of sediment, kg/m2
c = zeros(jmax,kmax); % initialize matrix, kg/m2s
conserve = zeros(2,1);

dtriv = [28,16]; % Position of Detroit River 
mmriv = [22,7]; % Position of Maumee River
ngriv = [20,106]; % Position of Niagara River
dtflux = 3.171e-08; % flux of Mirex, kg m-2 s-1
mmflux = 7.9274e-09; % flux of Mirex, kg m-2 s-1
ngflux = 36.1920; % flux of water, m2 s-1
setvelmir = 5e-7; % settling velocity of Mirex, m/s

% Grid Cell Size (m)
dx = 5000;
dy = 5000; 

% Timestep Size
dt = 10000; % s, don't go larger than 10,000

% Diffusivities in x and y directions
Kx = 100; % m2/s
Ky = 100; % m2/s 

% Advective Matrix
Tx = zeros(jmax,kmax); 
Ty = zeros(jmax,kmax);

% Courant number in x and y
Cx = U*dt/dx; 
Cy = V*dt/dy;

% Diffusion number in x and y
Gx = Kx*dt/dx^2; 
Gy = Ky*dt/dy^2;

alpha1 = ((Cx.^2)/6 - Cx/2 + Gx + 1/3).*Cx; 
alpha2 = (-(Cx.^2)/3 + Cx/2 + 5/6 - Cx.*Cy/2 ...
         -(Cy.^2)/2 + Cy/2 - 2*Gx - 2*Gy).*Cx; 
     
alpha3 = (-1/6 + (Cx.^2)/6 + Gx).*Cx;
alpha4 = (-Cy/2 + (Cy.^2)/2 + Gy).*Cx; 
alpha5 = (Cx.*Cy/2 + Gy).*Cx;

beta1 = ((Cy.^2)/6 - Cy/2 + Gy + 1/3).*Cy; 
beta2 = (-(Cy.^2)/3 + Cy/2 + 5/6 - Cx.*Cy/2 ...
         -(Cx.^2)/2 + Cx/2 - 2*Gx - 2*Gy).*Cy; 
     
beta3 = (-1/6 + (Cy.^2)/6 + Gy).*Cy;
beta4 = (-Cx/2 + (Cx.^2)/2 + Gx).*Cy; 
beta5 = (Cx.*Cy/2 + Gx).*Cy;

%% Solve for c
iterNum = ceil((3600*24*365)/dt)*10;
figure
for n = 1:iterNum
   % [DEAL WITH SOURCES HERE]
   % Add fluxes at correct positions
   % add it to the cell at each time, constant value each time,
   % concentration value*dt
   c(dtriv(1),dtriv(2)) =  dtflux*dt; % kg/m2
   c(mmriv(1),mmriv(2)) =  mmflux*dt;
% c(dtriv(1),dtriv(2)) = c(dtriv(1),dtriv(2))+ dtflux*dt; % kg/m2
% c(mmriv(1),mmriv(2)) = c(mmriv(1),mmriv(2))+ mmflux*dt;
   
    for j=2:jmax-2 
        for k=2:kmax-2
            % Set diffusivity flux to 0 at boundary
            Dfx = -Gx*c(j+1,k) + Gx*c(j,k);
            Dfy = -Gy*c(j,k+1) + Gy*c(j,k);
            if Depth(j,k) == 0
                Dfx = 0;
                Dfy = 0;
            elseif Depth(j+1,k) == 0 % look up, forward in rows to make sure its not land
                  Dfx = 0;
                  Dfy = -Gy*c(j,k+1) + Gy*c(j,k);
            elseif Depth(j,k+1) == 0 % look right, forward in columns to make sure its not land
                  Dfx = -Gx*c(j+1,k) + Gx*c(j,k);
                  Dfy = 0;
            end
            %%%this is done so that advection terms use the proper 
            %%%upwinding... NO NEED TO TOUCH THIS
           
            if U(j,k)>=0 && V(j,k)>=0
                Tx(j,k) = alpha1(j,k)*c(j+1,k) + alpha2(j,k)*c(j,k) + alpha3(j,k)*c(j-1,k) + alpha4(j,k)*c(j,k+1) + alpha5(j,k)*c(j,k-1) + Dfx;
                Ty(j,k) = beta1(j,k)*c(j,k+1) + beta2(j,k)*c(j,k) + beta3(j,k)*c(j,k-1)+ beta4(j,k)*c(j+1,k) + beta5(j,k)*c(j- 1,k) + Dfy;
                elseif U(j,k)<=0 && V(j,k)>=0
                Tx(j,k) = alpha1(j,k)*c(j,k) + alpha2(j,k)*c(j+1,k) + alpha3(j,k)*c(j+2,k)+ alpha4(j,k)*c(j+1,k+1) + alpha5(j,k)*c(j+1,k-1)+ Dfx;
                Ty(j,k) = beta1(j,k)*c(j,k+1) + beta2(j,k)*c(j,k) + beta3(j,k)*c(j,k-1)+ beta4(j,k)*c(j-1,k) + beta5(j,k)*c(j+1,k) + Dfy;
                elseif U(j,k)>=0 && V(j,k)<=0
                Tx(j,k) = alpha1(j,k)*c(j+1,k) + alpha2(j,k)*c(j,k) + alpha3(j,k)*c(j-1,k)+ alpha4(j,k)*c(j,k-1) + alpha5(j,k)*c(j,k+1)+ Dfx;
                Ty(j,k) = beta1(j,k)*c(j,k) + beta2(j,k)*c(j,k+1) + beta3(j,k)*c(j,k+2)+ beta4(j,k)*c(j+1,k+1) + beta5(j,k)*c(j-1,k+1) + Dfy;
                elseif U(j,k)<=0 && V(j,k)<=0
                Tx(j,k) = alpha1(j,k)*c(j,k) + alpha2(j,k)*c(j+1,k) + alpha3(j,k)*c(j+2,k)+ alpha4(j,k)*c(j+1,k-1) + alpha5(j,k)*c(j+1,k+1)+ Dfx;
                Ty(j,k) = beta1(j,k)*c(j,k) + beta2(j,k)*c(j,k+1) + beta3(j,k)*c(j,k+2)+ beta4(j,k)*c(j-1,k+1) + beta5(j,k)*c(j+1,k+1) + Dfy;
            end
            % [DEAL WITH ADVECTIVE BC HERE]
            diffTx = Tx(j-1,k) - Tx(j,k); 
            diffTy = Ty(j,k-1) - Ty(j,k);
            if Depth(j,k) == 0
                diffTx = 0;
                diffTy = 0;
            elseif Depth(j-1,k) == 0 % look forward in rows to make sure its not land
                diffTx = 0; 
                diffTy = Ty(j,k-1) - Ty(j,k);
            elseif Depth(j,k+1)== 0
                diffTx = Tx(j-1,k) - Tx(j,k);
                diffTy = 0; 
            end
            % [UPDATE C HERE]
            c(j,k) = c(j,k) + diffTx + diffTy;
            
            if Depth(j,k) == 0 % fixes so you don't divide by 0 when on land
                c(j,k) = 0;
            else
                % [DEAL WITH SINKS HERE] just subtract settling flux from every cell
                sediment(j,k) = (((c(j,k) * setvelmir)/Depth(j,k))*dt); % set mirex concentration
                c(j,k) = c(j,k) - sediment(j,k);
            end
        end  
    end

    c(ngriv(1),ngriv(2)) =  - ((ngflux*c(ngriv(1),ngriv(2)))*dt)/(dx*dy); % subtract niagra flux

% Plot Movies
    s = 1:100:10000;
    subplot(1,2,1)
    if ismember(n,s)
        pause(0.00001)
        imagesc(c)
        title("Time =" + n)
        set(gca,'YDir','normal')
        colorbar
        
        subplot(1,2,2)
        pause(0.00001)
        imagesc(sediment)
        title("Time =" + n)
        set(gca,'YDir','normal')
        colorbar
    end

end
%% 
toc