%% Exercise 3: 2D Diffusion: Fluid Flow Through Porous Media as a 2-D Diffusion Problem
% G582
% Caitlin Sifuentes & McKailey Sabaj
% 09/13/21

clear
close all
%% List Knowns
% Inputs
aq_thick = 10; % m, aquifer thickness
k = 10^-11; % m^2, rock permeability 
alpha = 0.2; % porosity
mu = 10^-3; % Pa*s, water viscosity
g = 9.81; % m/s^2, acceleration due to gravity

% Model Domain size
dx = 1000; % m , spatial step or grid spacing
domain_x = 40000; % m, 40 km
domain_y = 20000; % m, 20 km
domain_cols = domain_x/dx; % km
domain_rows = domain_y/dx; % km
bucolic_montana_property = zeros(domain_rows,domain_cols); % model domain

%% Calculate Boundary Data
upper_head = [80 60 100 80 90]; % m above datum, listed fluid heads (extrapolate)
lower_head = [90 70 60 30]; % m, (extrapolate last cell)
left_head = [80 70 90]; % m, (interpolate)
right_head = [90 50]; % m, (extrapolate last two cells)

% Interpolate and Extrapolate Data
upper_indices = round([0.48 7.68 19.04 25.44 35.52]); % indices of original, need to round to indicate specific positions
up_int = interp1(upper_indices,upper_head,1:36,'spline'); 
up_query = [37 38 39 40]; % query points that extend beyond the domain 
up_ext_ind = 1:36; % indices of interpolated data
up_ext = interp1(up_ext_ind,up_int,up_query,'spline','extrap');
upper_boundary = horzcat(up_int,up_ext);

lower_indices = round([0.8 12.48 22.88 31.84]); % indices of original 
low_int = interp1(lower_indices,lower_head,2:32,'spline'); 
low_query = [33 34 35 36 37 38 39 40]; % more qp
low_ext_ind = 2:32; % indices of interpolated data
low_ext1 = interp1(low_ext_ind,low_int,low_query,'spline','extrap');
low_ext2 = interp1(low_ext_ind,low_int,1,'spline','extrap');
low_ext = horzcat(low_ext1,low_ext2);
lower_boundary = horzcat(low_int,low_ext);

left_indices = round([0.48 8.64 18.4]); % indices of original 
left_int = interp1(left_indices,left_head,1:18,'spline');
left_query = [19 20]; % query points that extend beyond the domain 
left_ext_ind = 1:18;
left_ext = interp1(left_ext_ind,left_int,left_query,'spline','extrap');
left_boundary = horzcat(left_int,left_ext)';

right_indices = round([0.64 9.92]); % indices of original 
right_int = interp1(right_indices,right_head,1:10,'spline');
right_query = [11 12 13 14 15 16 17 18 19 20]; % query points that extend beyond the domain 
right_ext_ind = 1:10;
right_ext = interp1(right_ext_ind,right_int,right_query,'spline','extrap');
right_boundary = horzcat(right_int,right_ext)';

% Input Boundary Conditions
bucolic_montana_property(1,:) = upper_boundary;
bucolic_montana_property(:,1) = left_boundary;
bucolic_montana_property(20,:) = lower_boundary;
bucolic_montana_property(:,40) = right_boundary;

%% Calculate 2D Diffusion of Pollutant
% solve Poisson's Equation with an implicit method

% Solve for A Matrix
sqr_dom = 684; % number of unknown cells within the model domain
zero_ind = 38; % used to index position of 0 in A matrix
diag_index1 = zeros(sqr_dom,1); % intiate domain size of diagonal
A1 = diag(ones(size(diag_index1))-5); % place -4 in the diagonal
diag_index2 = zeros(sqr_dom-1,1); % intiate domain size of diagonal
for i = 1:length(diag_index2) % find all positions to place 0 in 1s & 0 diagonal
    if rem(i,zero_ind) == 0
        diag_index2(i) = 0;
    else
        diag_index2(i) = 1;
    end
end
A2 = diag(diag_index2,1); % create diag for above
A3 = diag(diag_index2,-1); % create diag for below
diag_index3 = zeros(sqr_dom-zero_ind,1); % intiate domain size of diagonal
A4 = diag(ones(size(diag_index3)),zero_ind); % create diag for above
A5 = diag(ones(size(diag_index3)),-zero_ind); % create diag for below
A = A1+A2+A3+A4+A5; % add all diagonals to create A matrix

% Solve for B Matrix
B = zeros(size(bucolic_montana_property,1)-2,size(bucolic_montana_property,2)-2); % initialize B matrix 
temp = zeros(1,4); % empty vector to run calculations
col_end = 39; % used to iterate through unknowns in model domain
row_end = 19; % used to iterate through unknowns in model domain
test_r = 17; % test well location
test_c = 19; % test well location
K = 1; % m/s

% For loop used to test if there is a boundary condition above, below,
% right, or left of the current cell, then place in 'B' matrix
% Also test whether or not the pumping rate (K) effectively works to stop
% the pollutant
for c = 2:col_end
    for r = 2:row_end
            if bucolic_montana_property(r-1,c) > 0 % up
           temp(1) = -bucolic_montana_property(r-1,c); 
            end
            if bucolic_montana_property(r+1,c) > 0 % down
           temp(2) = -bucolic_montana_property(r+1,c);
            end
            if bucolic_montana_property(r,c+1) > 0 % right 
           temp(3) = -bucolic_montana_property(r,c+1);
            end
            if bucolic_montana_property(r,c-1) > 0 % left
           temp(4) = -bucolic_montana_property(r,c-1);
            end
           B(r-1,c-1) = sum(temp);
           if c == test_c && r == test_r % set up if statement for well site
               B(r-1,c-1) = sum(temp)+K; % add pumping rate to well site
           end 
           temp = zeros(1,4); % empty vector to run calculations
    end
end
B = B';
B = reshape(B,[1,sqr_dom])';
X = A^-1*B; % calculate the unknowns
X = reshape(X,[col_end-1,row_end-1])';
figure
imagesc(X)
hold on 
%% Contours of the Potential Field
x = 1:col_end-1;
y = 1:row_end-1;
[X2,Y] = meshgrid(x,y);
contourf(X2,Y,X,20)
[dx,dy] = gradient(-X);
hold on 
quiver(X2,Y,dx,dy) % lengths of the vectors are proportional to the groundwater flow rates and the orientations of the vectors define streamlines of flow.
title('Groundwater Flow') 
xlabel('x (km)')
ylabel('y (km)')
a = colorbar;
a.Label.String = 'Hydraulic Head (m)';
hold on 

%% Plot Contaminated Well and Boundaries
cwell_x = [13.44 12 10.08];
cwell_y = [6.24	7.04 9.12];
c = scatter(cwell_x,cwell_y,'*','b');
hold on 
nodrill_x = [11.68	8	4.16	2.72	1.6	0.64	0.64	1.28	2.4	4	6.08	8	10.24	15.68	19.52	21.44	23.52	25.6	27.52	28.32	29.12	29.44	29.44	28.8	27.52	25.76	23.52	19.52	15.68	11.68];
nodrill_y = [0.32	0.96	3.04	3.84	5.76	7.84	9.92	11.84	13.76	15.52	16.32	16.8	16.64	16.48	16.16	15.52	15.2	14.24	12.96	11.84	10.4	9.28	7.84	6.08	4.8	4	3.04	1.6	0.8	0.32];
d = plot(nodrill_x,nodrill_y,'k');
streamline(X2,Y,dx,dy,cwell_x(1),cwell_y(1))
streamline(X2,Y,dx,dy,cwell_x(2),cwell_y(2))
streamline(X2,Y,dx,dy,cwell_x(3),cwell_y(3))
Suggested_WellX = test_c;
Suggested_WellY = test_r;
hold on
s = scatter(Suggested_WellX,Suggested_WellY,'r','*');





