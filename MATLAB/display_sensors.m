function [ ] = display_sensors( x,mask, pivot )
%display_sensors Template for displaying image data
%   Displays enso sensors on white w/ black continents

    
    sensors = zeros(360,180);
    P = zeros(size(x)); P(pivot)=1:length(pivot);    
    sensors(mask==1) = P;
       
    %mask(mask==1)=x;
    b = imagesc(mask');
    shading flat, colormap(gray), drawnow
    
    hold on
    
    %scale x for colormap
    x = x+abs(min(x));
    x = x/max(x);
    
    % associate colors with each element of x
    cmap = squeeze(ind2rgb(uint8(256*x),jet(256)));
    cmap = cmap(pivot,:);
    
    S = reshape(real(sensors)',360*180,1);
    [C,IC,~] = unique(S);
    
    % align Ilin with pivot somehow
    
    [I,J] = ind2sub(size(sensors'),IC(2:end));
    
    scatter(J,I,30,cmap,'filled','MarkerEdgeColor','k');
        
    
    axis off
end

