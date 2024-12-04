% create model setup from TIFF image

function [units,D,Nz] = ModelFromImage(filename,n_units,W,Nx)

% read in RGB image from TIFF file
img = double(importdata(filename));

% get size of image
[p,q,k] = size(img); % p= vertical, q = horizontal, k = colour resolution 

% identify units by clustering analysis
[Ic,Fc] = kmeans(reshape(img,p*q,k),9,'MaxIter',1e3,'Replicates',5);

% sort clusters for ascending SiO2 content
[~,isort]   = sort(Fc(:,1)-sum(Fc,2),'descend');
Ics = zeros(size(Ic));
for ic = 1:n_units
    Ics(Ic==isort(ic)) = ic;
end
imgc = reshape(Ics,p,q);

% interpolate from original dimensions to target model size
D  = W*p/q;
Nz = floor(Nx*p/q);

ho  = W/q;
xco = ho/2:ho:W-ho/2;
zco = ho/2:ho:D-ho/2;
[Xco,Zco] = meshgrid(xco,zco);

h   = W/Nx;
xc  = h/2:h:W-h/2;
zc  = h/2:h:D-h/2;
[Xc,Zc] = meshgrid(xc,zc);

imgi = interp2(Xco, Zco, imgc, Xc, Zc);
figure(1); clf

% Define the discrete color levels
numColors = 9; % Number of distinct colors
cmin = min(imgc(:));
cmax = max(imgc(:));
cLevels = linspace(cmin, cmax, numColors + 1); % Define boundaries for each color level

% Create a colormap with 9 colors
cmap = parula(numColors); % Use any colormap with exactly 9 colors

% Subplot 1: Original Data
subplot(2, 1, 1)
imagesc(xco, zco, imgc); 
axis equal tight;
colorbar;
caxis([cmin, cmax]); % Set the color axis range
colormap(cmap); % Apply the discrete colormap
xlabel('x (m)', 'FontSize', 12)
ylabel('z (m)', 'FontSize', 12)

% Subplot 2: Interpolated Data
subplot(2, 1, 2)
imagesc(xc, zc, imgi); 
axis equal tight;
cb = colorbar;
caxis([cmin, cmax]); % Keep the same range for consistency
colormap(cmap); % Apply the same colormap
xlabel('x (m)', 'FontSize', 12)
ylabel('z (m)', 'FontSize', 12)

% Adjust colorbar ticks for discrete appearance
cb.Ticks = linspace(cmin, cmax, numColors); % Set ticks to match the number of colors
cb.TickLabels = round(cb.Ticks, 2); % Optional: Customize tick labels

%imgi = interp2(Xco,Zco,imgc,Xc,Zc);
%figure(1); clf
%subplot(2,1,1)
%imagesc(xco,zco,imgc); axis equal tight; colorbar
%xlabel('x (m)','Fontsize',12)
%ylabel('z (m)','Fontsize',12)
%subplot(2,1,2)
%imagesc(xc,zc,imgi); axis equal tight; colorbar
%xlabel('x (m)','Fontsize',12)
%ylabel('z (m)','Fontsize',12)

units = uint8(imgi);

end