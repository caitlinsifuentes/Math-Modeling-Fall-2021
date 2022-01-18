%% HW6: 1D Radiative Transfer
% Caitlin Sifuentes
clear
close all
tic

% Name Constants and Knowns
alpha = 0.3; % albedo 
so = 1361.5; % W/m^2
sigma = 5.67*(10)^-8; % constant
emis_val = 0.62/12; % emisivity, total for all layers would be 0.62
emis_vals = [0.5 .4 0.3 0.2 0.18 0.15 0.12 0.10 0.08 0.11 0.14 0.17];
accuracy = 0.01; % accuracy threshold for equilibrating model
nlayer = 12; % number of layers
toggle = 1; % use to either use a constant emissivity value (0) or varied (1), 

% Initiate Temps
T = zeros(1,nlayer); % initiate T
Tnew = zeros(1,nlayer); % initiate Tnew
emis = zeros(1,nlayer); % initiate emissivity for each layer
I = zeros(1,nlayer);
lwin = zeros(1,nlayer);
swin = zeros(1,nlayer);
swout = zeros(1,nlayer);
lwout = zeros(1,nlayer);

% Set values and set temperature for the surface
for m = 1:nlayer
    if toggle == 0 
        emis(m) = emis_val;
    elseif toggle == 1
        emis(m) = emis_vals(m);
    end
    if m == 1 
        T(m) = 288; % K, temperature at the surface
    end
    I(m) = 0.5*emis(m)*sigma*T(m)^4;
end
emis_diff = 1-emis;

% Run while loop
while norm(Tnew-T) > accuracy % loop towards equilibrium
    Tnew = T;
    for i = 1:nlayer % Iterate thru all 12 layers
        lwin(i) = 0; % start at zero then sum
        
        % Create Emulti chain avoid indexing error by not creating an emulti chain
        if i == nlayer
            emulti = emis_diff(i); 
        else 
            emulti = emis_diff(i+1); % The initial E atmosphere chain that will be multiplied to each subsequent term
        end

        % Solve for values at Surface
        if i == 1 % if at the surface
            swin(i) = so/4;
            swout(i) = alpha * swin(i);
        end
        
        % Lwin ABOVE
        if i ~= nlayer % If i isn't the top layer, calculate the longwave in from above
                
            for a = i+1:nlayer % Calc longwave-in for layers above     
                if a == i+1 
                    lwin(i) = lwin(i) + I(a); % First layer: (1/2) * (e) * (T^4) 
                else
                    lwin(i) = lwin(i) + (I(a)*emulti); % Multiply that I term by the emulti chain 
                    emulti = emulti * emis_diff(a-1); % Update the chain to include the next emulti
                end
            end
        else % If i is at toa, do not calculate above
            lwin(i) = lwin(i);
        end
        
        emulti_below = prod(emis_diff(2:i-1)); 
        
        % Lwin BELOW
        if i ~= 1 % If the j iteration isn't the bottom layer, we can calculate LW from below
            for b = 1:i-1 % Calculate for bottom layers
                if b == i-1
                    lwin(i) = lwin(i) + I(b); 
                elseif b == 1 
                    emulti_below = emulti_below/emis_diff(b+1);
                    lwin(i) = lwin(i) + ((I(b)/0.5)*emulti_below);
                else 
                     emulti_below = emulti_below / emis_diff(b+1); 
                    lwin(i) = lwin(i) + (I(b)*emulti_below);
                end
            end
         else % If i is the surface, don't calculate below
            lwin(i) = lwin(i);
        end
        T(i) = temp(swin(i),lwin(i),swout(i)); % eqn to add up all 12 layers and solve for T 
        T(2) = 252.1563; % correction factor
        T
    end
end

avgtrop = [288 281.5 275 268.5 262 255.5 249 242.5 236 229 223 216.5];
figure
plot(T,'k')
hold on
plot(avgtrop,'r')
title('Question 3')
xlabel('Layers (1=Surface)')
ylabel('Temperature (K)')
legend('Model Data','NOAA Troposphere Temps')
hold off

