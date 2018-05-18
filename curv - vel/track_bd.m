function boundary_traj = track_bd(ifiledir,DateTypeString,ComputerString,frameNUM,resolution,edge_cutoff)

% This function track the boundary and treating them as particles
% 
% Input: 
%  - ifiledir: full path of the input file folder;
%  - DataTypeString & ComputerString: Movie ID string;
%  - frameNUM: image ID for time course analysis;
%  - resolution: (optional) smooth size of contour cordinates; default
%  value = 30;
%  - edge_cutoff: (optional) exclude points close to the edge; default = 5;
%
% Output: 
% - boundary_traj: structure variable with one field totalVertices.
% Organization of 'totalVertices' is as follows
% column 1,2: x,y position of boundary points
% column 3: signed curvature
% column 4:5 surface normal of local points. 
% 
% * Use the commented code chunk to plot a beautiful flower of contour
% evolution as shown in main text Fig.4c!!!

% mycolor =  jet(numel(frameNUM));

% Set default values for optional input parameter.
    if ~exist('resolution','var')
        resolution = 30;
    end
    if ~exist('edge_cutoff','var')
        edge_cutoff = 5;
    end

for t = 1:length(frameNUM)
    frame = frameNUM(t);
    totalVertices=[];
    [sz,totalVertices] = track_crv(ifiledir,DateTypeString,ComputerString,frame,resolution);
    
    % Use this part, if you want to exlude all points that are too close to the edge.
    if exist('edge_cutoff','var')
        keep=(totalVertices(:,1)>edge_cutoff) & (totalVertices(:,2)>edge_cutoff) & (totalVertices(:,2)<sz(1)-edge_cutoff) & (totalVertices(:,1)<sz(2)-edge_cutoff);
        totalVertices=totalVertices(keep,:);
    end
    
    boundary_traj(t).totalVertices=totalVertices;
    disp(t)
    
%     if t == 1
%         f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frame,'%03d'),'.jpg'));
%         gg = mat2gray(f(:,:,3));    
%         [counts,binlocs] = imhist(gg);
%         [pks,idx] = findpeaks(-counts);binlocs(idx(end));
%         [~,mind] = max(counts);
%         cutoff = find(idx < mind,1,'last');
%     
%         g = gg;
%         if ~isempty(cutoff)
%             g = g(g>=binlocs(idx(cutoff)));
%         else
%             g = g(:);
%         end
%         hx = histfit(g);
%         center_bg = mean(hx(2).XData);
%         gg = mat2gray(f(:,:,3));
%         gg = imadjust(gg,[min(min(gg)) (center_bg - double(min(min(gg))))/0.9 + double(min(min(gg)))],[]);
%         
%         imshow(gg)
%     end
%     plot(totalVertices(:,1),totalVertices(:,2),'Linewidth',5,'Color',mycolor(t,:));hold on  
          
end
% hold off