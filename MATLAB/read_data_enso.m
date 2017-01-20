function [ Lat, Lon, time, mask, X ] = read_data_enso( filename, maskname )
%read_data Read in ENSO data
%   Return meshgrid of lat,long,times and sea-surface temp

X = ncread(filename,'sst');
time = ncread(filename,'time'); %days
%time = (time-time(1))/7; % time in weeks

lat = ncread(filename,'lat');
lon = ncread(filename,'lon');

mask = ncread(maskname,'mask');
mask = mask(:,:,1);
[Lat,Lon] = meshgrid(lat,lon);

% % Do not plot continent mask
% Lat(mask==0) = NaN;
% Lon(mask==0) = NaN;
end

