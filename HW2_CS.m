%% Exercise 2: 1D Diffusion: Predicting Soil Temperature
% G582
% Caitlin Sifuentes
% 09/01/21

clear
close all
%% Load the data
data = xlsread('SoilTdata.xlsx');
depth = data(2:5,1); % ft

% Observations
dates = data(1,2:14); % month.day
soil_temp_ob = data(2:5,2:14); % fahrenheit

% Constant
D = 0.06*10.764; % ft2/day, thermal diffusivity, constant

% Initial Conditions
soil_temp = data(2:5,2); % left column IC
soil_temp_depth = data(2:5,1);

% Variables
dt = 2; % days, timestep
dx = 2; % ft, space step 
totaltime = 280; % total run time
totaldepth = 40; % total distance in ft
s = D*(dt/(dx^2)); % S needs to be <0.5 to be stable

% Boundary Conditions
air_temp = data(8,2:14); % fahrenheit , upper boundary
air_temp_day = data(7,2:14); % time for air temp measurement
query_points_at = (1:1:totaltime/dt);
query_points_st = (1:1:totaldepth/dx);
cave_lowb = zeros(1,totaltime/dt); % initialize lower boundary
cave_lowb = cave_lowb+50; %  % lower boundary condition, fahrenheit, yearly average temp of deep caves

% Interpolate Boundary Conditions
% use interp1 
airtemp_interp = interp1(air_temp_day/dt,air_temp,query_points_at,'pchip'); % sample points, sample values, query points to int thru, using cubic (3rd order)
soiltemp_interp = interp1(soil_temp_depth/dx,soil_temp,query_points_st,'pchip')'; % interp oct soil values for y axis IC
air_time = linspace(1,totaltime,totaltime/dt);
soil_depth = linspace(1,totaldepth,totaldepth/dx);

% Plot Interpolated Data
subplot(2,1,1)
plot(air_time,airtemp_interp)
title('Upper Boundary Surface Temps')
xlabel('Timesteps (2 days)')
ylabel('Air Temperature (F)')
grid on
subplot(2,1,2)
plot(soil_depth,soiltemp_interp)
title('Leftmost Boundary October Soil Temps')
xlabel('Distance (2ft)')
ylabel('Soil Temperature (F)')
grid on 

%% Intialize Matrix M
M = zeros(round(totaldepth/dx),round(totaltime/dt));

% Fill Matrix M with Boundary and Initial Conditions
M(:,1) = soiltemp_interp;
M(1,:) = airtemp_interp;
M(20,:) = cave_lowb;

% Run FTCS Finite Diff Eq
for k = 1:totaltime/dt-1
    for j = 2:totaldepth/dx-1
        M(j,k+1) = s*M(j+1,k)+(1-2*s)*M(j,k)+s*M(j-1,k);
    end
end

figure
imagesc(M)
title('1D Soil Temperature Diffusion Model')
xlabel('Time (4days)')
ylabel('Depth (2ft)') 
a = colorbar;
a.Label.String = 'Temperature (Dg Fahrenheit)';

%% Test Statistical Similarity of Model Data and Observed Data
% "The hypothesis will be falsified if model predictions are not
% statistically similar to these data."

% RESIZE MATRIX
observed_values = data(2:5,2:14); % select observed values from data
modeled_vals = zeros(length(soil_temp_depth),length(air_temp_day)); % initialize new matrix

% Resize Modeled Matrix to exact days and depths of Observed to have equal sizes
for c = 1:length(air_temp_day)
    for r = 1:length(depth)
            modeled_vals(r,c) = M(round(depth(r)/dx),round(air_temp_day(c)/dt));
    end
end

% STATISTICAL ANALYSIS
dim = size(modeled_vals);
N = (dim(1)*dim(2)); % sample size for RMSE
rmse = (sum(sum((observed_values - modeled_vals).^2)/N)).^0.5;
stdev = std(modeled_vals);
rval = corrcoef(modeled_vals,observed_values);
rsquare = rval.*rval;

% Plot a linear regression for every depth
figure
subplot(2,2,1)
scatter(modeled_vals(1,:),observed_values(1,:))
c = polyfit(modeled_vals(1,:),observed_values(1,:),1);
y_est = polyval(c,modeled_vals);
hold on
plot(modeled_vals,y_est,'r')
rval1 = corrcoef(modeled_vals(1,:),observed_values(1,:));
rsquare1 = rval1.*rval1;
txt1 = strcat('r^2 = ',num2str(rsquare1(2,1)));
text(40,57,num2str(txt1))
hold off
title('Regression for Depth = 5ft')
xlabel('Modeled Temperature Values')
ylabel('Observed Temperature Values')
grid on

subplot(2,2,2)
scatter(modeled_vals(2,:),observed_values(2,:))
c2 = polyfit(modeled_vals(2,:),observed_values(2,:),1);
y_est2 = polyval(c2,modeled_vals);
hold on
plot(modeled_vals,y_est2,'r')
rval2 = corrcoef(modeled_vals(2,:),observed_values(2,:));
rsquare2 = rval2.*rval2;
txt2 = strcat('r^2 = ',num2str(rsquare2(2,1)));
text(40,57,num2str(txt2))
hold off
title('Regression for Depth = 10ft')
xlabel('Modeled Temperature Values')
ylabel('Observed Temperature Values')
grid on

subplot(2,2,3)
scatter(modeled_vals(3,:),observed_values(3,:))
c3 = polyfit(modeled_vals(3,:),observed_values(3,:),1);
y_est3 = polyval(c3,modeled_vals);
hold on
plot(modeled_vals,y_est3,'r')
rval3 = corrcoef(modeled_vals(3,:),observed_values(3,:));
rsquare3 = rval3.*rval3;
txt3 = strcat('r^2 = ',num2str(rsquare3(2,1)));
text(40,57,num2str(txt3))
hold off
title('Regression for Depth = 20ft')
xlabel('Modeled Temperature Values')
ylabel('Observed Temperature Values')
grid on

subplot(2,2,4)
scatter(modeled_vals(4,:),observed_values(4,:))
c4 = polyfit(modeled_vals(4,:),observed_values(4,:),1);
y_est4 = polyval(c4,modeled_vals);
hold on
plot(modeled_vals,y_est4,'r')
rval4 = corrcoef(modeled_vals(4,:),observed_values(4,:));
rsquare4 = rval4.*rval4;
txt4 = strcat('r^2 = ',num2str(rsquare4(2,1)));
text(40,57,num2str(txt4))
hold off
title('Regression for Depth = 30ft')
xlabel('Modeled Temperature Values')
ylabel('Observed Temperature Values')
grid on


