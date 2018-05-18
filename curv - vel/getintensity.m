function intensity = getintensity(totVert,ifiledir,DateTypeString,ComputerString,frame,bandw,arcangle)
% This function is used to track the local intensity from the movie image.
% It is called by 'connect boundary'. See more input details there.
%
% Output:
%  - intensity: averaged local intensity;
% 
% Optional intput:
%  - bandw & arcangle: define the bandwidth and angle (in arc) to average
%  local intensity. Defaul value is given as follows:

if ~exist('bandw','var')
    bandw = 100;
end
if ~exist('arcangle','var')
    arcangle = 0.05;
end

% Preprocessing image (same technique used in other functions, see e.g.
% function 'MovieMaker'.
    f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frame,'%03d'),'.jpg'));  
    gg = mat2gray(f(:,:,3));    
    [counts,binlocs] = imhist(gg);
    [pks,idx] = findpeaks(-counts);binlocs(idx(end));
    [~,mind] = max(counts);
    cutoff = find(idx < mind,1,'last');
    
    g=gg;
    if ~isempty(cutoff)
        g = g(g>=binlocs(idx(cutoff)));
    else
        g = g(:);
    end
    hx = histfit(g);
    center_bg = mean(hx(2).XData);
    SD_bg  = (max(hx(2).XData) - min(hx(2).XData))/6;
    gg = mat2gray(f(:,:,3));
    gg = imadjust(gg,[min(min(gg)) (center_bg - double(min(min(gg))))/0.9 + double(min(min(gg)))],[]); % post-processed image

[rows, columns] = size(gg);   
intensity = zeros(1,size(totVert,1));    
for i = 1:size(totVert,1)
    
    % define a triangular mask to average data
    [angle,~] = cart2pol(totVert(i,4),totVert(i,5));
    xBee = totVert(i,1); yBee = totVert(i,2);
    angle1 = angle + arcangle;
    angle2 = angle - arcangle;
    
    y1 = yBee - bandw * sin(angle1);
    x1 = xBee - bandw * cos(angle1);
    y2 = yBee - bandw * sin(angle2);
    x2 = xBee - bandw * cos(angle2);
    
    yTriangle = [yBee, y1, y2, yBee];
    xTriangle = [xBee, x1, x2, xBee];
    mask = poly2mask(xTriangle, yTriangle, rows, columns);
    
    avgregionInt = gg(mask>0.5);
    intensity(i) = mean(avgregionInt);
    
%     if mod(i,100) == 1
%         disp(i)
%     end

end