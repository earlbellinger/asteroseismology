saveplot = 1;

%% Define the spherical grid

theta = 0:pi/200:pi;                   % polar angle
phi = 0:pi/100:2*pi;                   % azimuth angle

[phi,theta] = meshgrid(phi,theta);    % define the grid

%% Calculate the Spherical Harmonic

for degree = [0 1 2 3 4 5] %20 75]

%degree = 3;
order = 0;
amplitude = 30; %15/degree;
radius = 7;

Ymn = legendre(degree,cos(theta(:,1)));
Ymn = Ymn(order+1,:)';
yy = Ymn;

for kk = 2: size(theta,1)
    yy = [yy Ymn];
end

yy = yy.*cos(order*phi);  

order2 = max(max(abs(yy)));
rho = radius + amplitude*yy/order2;

r = rho.*sin(theta);    % convert to Cartesian coordinates
x = r.*cos(phi);
y = r.*sin(phi);
z = rho.*cos(theta);

%% Plot the Spherical Harmonic on the Surface of a Sphere

f = figure('Position', [0, 0, 1000, 1000]);%450, 500]);
if saveplot; set(f, 'Visible', 'off'); end
rabs = r;
rabs(rabs>=0) = rabs(rabs>=0)+10;
rabs(rabs<0) = rabs(rabs<0)-10;
s = surf(x,y,z,rabs);

h = light               %z add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
%if degree == 2; lightangle(h, 100, 30); end
%if degree == 3; lightangle(h, 50, 30); end
%if degree == 4; lightangle(h, 0, 30); end
%if degree == 5; lightangle(h, -50, 30); end
%lightangle(h, -35, 30);
%lighting gourad 
view(10,10)         % set viewpoint
colormap winter%(2)
%colormap([[.976 .443 0]; [.02 .443 .69]]);
shading interp
%camzoom(1.5)        % zoom into scene

%% Animate the Surface

scale = [linspace(0, 1, 10) linspace(1, -1, 20) linspace(-1, 0, 10)];

path = strcat('plots/plot_sph_harm/');
if saveplot; mkdir(path); end

%for ii = 1:(length(scale)-1)
ii = 11%length(scale)/2
    
    rho = radius + scale(ii)*amplitude*yy/order2;   
   
    r = rho.*sin(theta);
    x = r.*cos(phi);       
    y = r.*sin(phi);
    z = rho.*cos(theta);
    
    s.XData = x;    % replace surface x values
    s.YData = y;    % replace surface y values
    s.ZData = z;    % replace surface z values
    
    if not(saveplot); pause(0.1); end
    
    if saveplot
        %set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.5 4]);
        %saveas(f, strcat(path, num2str(degree), '.pdf'))
        
        %set(gcf,'Units','Inches');
        %pos = get(gcf,'Position');
        %set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        %saveas(gcf, strcat(path, num2str(degree), '.png'))
        
        gcf.PaperPositionMode = 'auto';
        print(strcat(path, num2str(degree)), '-dpng', '-r0')
    end
%end

end

