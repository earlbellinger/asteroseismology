saveplot = 0;

%% Define the spherical grid

theta = 0:pi/40:pi;                   % polar angle
phi = 0:pi/20:2*pi;                   % azimuth angle

[phi,theta] = meshgrid(phi,theta);    % define the grid

%% Calculate the Spherical Harmonic

for degree = 0:3

%degree = 3;
order = 0;
amplitude = 20/(degree+1);
radius = 10;

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

f = figure
if saveplot; set(f, 'Visible', 'off'); end
s = surf(x,y,z);

light               % add a light
lighting gouraud    % preferred lighting for a curved surface
axis equal off      % set axis equal and remove axis
view(40,30)         % set viewpoint
colormap parula
camzoom(1.5)        % zoom into scene

%% Animate the Surface

scale = [linspace(0, 1, 10) linspace(1, -1, 20) linspace(-1, 0, 10)];

path = strcat('plots/anim_sph_harm/', num2str(degree), '_', ...
        num2str(order), '/');
mkdir(path);

for ii = 1:(length(scale)-1)
    
    rho = radius + scale(ii)*amplitude*yy/order2;   
   
    r = rho.*sin(theta);
    x = r.*cos(phi);       
    y = r.*sin(phi);
    z = rho.*cos(theta);
    
    s.XData = x;    % replace surface x values
    s.YData = y;    % replace surface y values
    s.ZData = z;    % replace surface z values
    
    if not(saveplot); pause(0.1); end
    
    if saveplot; saveas(f, strcat(path, num2str(ii), '.png')); end
end

end

