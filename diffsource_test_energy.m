%***** 2D diffusion model of heat transport *******************************

%*****  Initialise Model Setup

% create x-coordinate vectors
xc = h/2:h:W-h/2;      % x-coordinate vector for cell centre positions [m]
zc = h/2:h:D-h/2;      % z-coordinate vector for cell centre positions [m]
xf = 0:h:W;            % x-coordinate vectore for cell face positions [m]
zf = 0:h:D;            % z-coordinate vectore for cell face positions [m]
[Xc,Zc] = meshgrid(xc,zc);  % create 2D coordinate arrays



% set up index array for boundary conditions
ix = [ 1,1:Nx,Nx ];  % closed/insulating sides
iz = [ 1,1:Nz,Nz ];  % closed/insulating top, flux grad at bottom

% set initial condition for temperature at cell centres
T   = T0 + dTdz(2).*Zc;  % initialise T array on linear gradient
Ta = T;

% set up condition for air
air = units == 9;

%*****  Solve Model Equations

t = 0;  % initial time [s]
tau = 0;  % initial time step count
dt = CFL * (h/2)^2/max(k0, [], 'all');

t_vals = [];  % array for time
E_vals = [];  % array for energy
T_vals = [];
rCV_vals = [];

while t <= tend

    % increment time and step count
    t = t+dt;
    tau = tau+1;

    % reset air section
    %T(air) = Tair;

            % 4th-order Runge-Kutta time integration scheme
                    
            dTdt1 = diffusion(T,dTdz,            k0,h,ix,iz);
            dTdt2 = diffusion(T+dTdt1/2*dt, dTdz,k0,h,ix,iz);
            dTdt3 = diffusion(T+dTdt2/2*dt, dTdz,k0,h,ix,iz);
            dTdt4 = diffusion(T+dTdt3  *dt, dTdz,k0,h,ix,iz);
        
            T = T + (dTdt1 + 2*dTdt2 + 2*dTdt3 + dTdt4)/6 * dt;
     
            % Calculate E
            rCV = sum(rho(:) .* Cp(:) * h*h);
            E = sum(T(:)) * rCV; 
            E_vals = [E_vals, E];  
            t_vals = [t_vals, t]; 
            T_vals = [T_vals,sum(T(:))];
            %rCV_vals = [rCV_vals,rCV];
           
    %plot model progress every 'nop' time steps
   % if ~mod(tau,nop)
   %     makefig(xc,zc,T,t,yr);
   % end

   
end


%*****  calculate numerical error norm
%Errx = norm(T - Ta,1)./norm(Ta,2);
%disp(' ');
%disp(['Numerical error = ',num2str(Errx)]);
%disp(' ');
%
%Errz = norm(T - Ta,1)./norm(Ta,1);
%disp(' ');
%disp(['Numerical error = ',num2str(Errz)]);
%disp(' ');

%*****  Utility Functions  ************************************************

% Function to calculate diffusion rate
function [dTdt] = diffusion(f, dTdz, k0, h, ix, iz)

% average k0 values to get cell face values
kx = k0(:, ix(1:end-1)) + k0(:, ix(2:end))/2;
kz = k0(iz(1:end-1), :) + k0(iz(2:end), :)/2;

% calculate heat flux by diffusion
qx = - kx .* diff(f(:, ix), 1, 2)/h;
qz = - kz .* diff(f(iz, :), 1, 1)/h;

% set boundary conditions
qz(end, :) = -kz(end, :) .* dTdz(1);

qx(:,1) = 0;  % no heat flux at the left boundary
qx(:,end) = 0;  % no heat flux at the right boundary
qz(1,:) = 0;  % no heat flux at the top boundary
qz(end,:) = 0;  % no heat flux at the bottom boundary

% calculate flux balance for rate of change
dTdt = -(diff(qx, 1, 2)/h + diff(qz, 1, 1)/h);

end


plot(t_vals/yr,E_vals,'m','Linewidth', 2);   
xlim([min(t_vals)/yr max(t_vals)/yr]);
ylim([1e19 6e20]);
xlabel('Time [yr]','FontSize',18, 'FontName','Times New Roman');
ylabel('Energy [J]','FontSize',18, 'FontName','Times New Roman');
title('Total Thermal Energy','FontSize',20, 'FontName','Times New Roman');
grid on;



% Function to make output figure
%function makefig(x,z,T,t,yr)
%
%clf; 
%
%%plot temperature
%imagesc(x,z,T); axis equal tight; colorbar; hold on
%ylabel('z [m]','FontSize',18, 'FontName','Times New Roman')
%xlabel('x [m]','FontSize',18, 'FontName','Times New Roman')
%ylabel(colorbar,'Temperature [°C]','FontSize',18, 'FontName','Times New Roman')
%title(['Temperature; time = ',num2str(t/yr),' yr'],'FontSize',20, 'FontName','Times New Roman')
%[C,h] = contour(x,z,T, [50,100,150], 'r', 'Linewidth', 2);
%clabel(C,h,'Fontsize',18,'Color','r', 'FontName','Times New Roman');
%
%drawnow;
%
%end