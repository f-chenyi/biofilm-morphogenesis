function [boundary_traj]=connect_boundary(boundary_traj,ifiledir,DateTypeString,ComputerString,imgid)

% This function links the tracked interface point by the simplest nearest
% neighbor assumption, and extract local intensity from the image. By doing
% so, we are able to quantify the colocalization between morphological
% features and local curvature and velocity.
%
% Input:
%  - boundary_traj: basic struture variable from 'track_bd';
%  - ifiledir: full path of the input file folder;
%  - DataTypeString & ComputerString: Movie ID string;
%  - imgid: image ID used to extract local intensity
% 
% Output:
%  - boundary_traj: enriched structure variable; first 5 column same as the
%  input; additional columns are
% column 6: particle ID for the starting frame; for other frames, it
% records the index of linked particle from the previous frame;
% 
% column 7: the distance from the particle it is linked (in the previous
% frame). (distance)= (local velocity) x (time interval between frames)
%
% column 8: the particle index it is linked to in the next frame;
% 
% column 9: angle used to specify the expected moving direction;
%
% * Note: use commented code chuck with '***' to add constrains to the
% linking process. The scheme is to first pick out the 10 nearest
% candidates in the next frame, and then link the current point to the one
% closest to the expected moving direction. Due to geometrical continuity
% requirement, the expected moving direction is set by the average of the
% neighboring 5 points.


  tempvertice = boundary_traj(1).totalVertices(3:5:end-2,:); % [1,2,3,4,5] collectively form a new spot
  tempvertice(:,6) = 1:size(tempvertice,1); 
  tempvertice(:,7) = zeros(size(tempvertice,1),1);
  
  [tempvertice(:,9),~] = cart2pol(boundary_traj(1).totalVertices(3:5:end-2,4),boundary_traj(1).totalVertices(3:5:end-2,5));
%   % *** modified version
%   [theta,~] = cart2pol(boundary_traj(1).totalVertices(1:5*size(tempvertice,1),5),boundary_traj(1).totalVertices(1:5*size(tempvertice,1),4));
%   tempvertice(:,9) = mean(reshape(theta,[size(tempvertice,1) 5]),2); 
 
  boundary_traj(1).totalVertices = tempvertice;
  
  
for frame=1:size(boundary_traj,2)-1
    
    totalVertices=boundary_traj(frame).totalVertices;
    totalVertices2=boundary_traj(frame+1).totalVertices;
    totalVertices2(:,6)=0;
    totalVertices2(:,7)=0;
    
    for n=1:size(totalVertices,1)
        temp=totalVertices(n,:);
        [X,~]=meshgrid(temp(1,1:2),1:size(totalVertices2,1));
        temp2=totalVertices2(:,1:2)-X;
        [temp2(:,3),temp2(:,4)]=cart2pol(temp2(:,1),temp2(:,2));
        
        [C,ind] =  min(temp2(:,4));
        
%         % *** modified version
%         sortdist = sort(temp2(:,4));
%         indexset = find(temp2(:,4)<=sortdist(10));
%         [~,i] = min(abs(sin(temp2(indexset,3)-temp(1,9))));
%         ind = indexset(i); C = temp2(ind,4);
               
            if totalVertices2(ind,6)==0 % if the target point has not been linked, record it
                totalVertices2(ind,6)=n;
                totalVertices2(ind,7)=C;
                totalVertices(n,8)=ind;
            else if C < totalVertices2(ind,7) % if the target point has been linked to another particle, need to compare which one makes more sense. 
                removeidx = totalVertices2(ind,6);
                totalVertices2(ind,6)=n;
                totalVertices2(ind,7)=C; 
                totalVertices(n,8)=ind;
                totalVertices(removeidx,8)=0;
                end
            end

    end

    
% Use following code if you want to take care of the missing particles. 
%     if ~isempty(totalVertices(:,8)==0)
%         for n=1:size(totalVertices,1) 
%             if totalVertices(n,8)==0
%                 idxset = find(totalVertices2(:,6)==0); % find all the particles that are not linked yet in the second frame
%                 temp=totalVertices(n,:);
%                 [X,~]=meshgrid(temp(1,1:2),1:length(idxset));
%                 temp2=totalVertices2(idxset,1:2)-X;
%                 [temp2(:,3),temp2(:,4)]=cart2pol(temp2(:,1),temp2(:,2));
%                 idxtemp = find(abs(temp2(:,3)-temp(1,9)) < angle_cut & temp2(:,4) < cutoff);               
%                 if ~isempty(idxtemp)
%                     [dist,i] = min(temp2(idxtemp,4));
%                     idxassigned = idxset(idxtemp(i));
%                     totalVertices(n,8)=idxassigned;
%                     totalVertices(idxassigned,6)=n;
%                     totalVertices(idxassigned,6)=dist;
%                 end              
% %            temp2(:,5)=1:size(totalVertices2,1);
% %            temp2=sortrows(temp2,4);
% %            ind=floor(mean(temp2(1:average_over,5)));
% %            if find(temp2(:,5)==ind)<cutoff
% %                totalVertices(n,8)=ind;
% %            end
%             
% %             totalVertices2(n,6)=lost;
% %             lost=lost+1;
% 
%             end
%         end
%     end
    

    tempvertices2 = totalVertices2(totalVertices2(:,6)~=0,:);
    tempvertices2(:,8) = zeros(size(tempvertices2,1),1);
    [theta2,~] = cart2pol(totalVertices2(:,4),totalVertices2(:,5));
    tempvertices2(:,9) = theta2(totalVertices2(:,6)~=0);
    
    boundary_traj(frame+1).totalVertices=tempvertices2;
    boundary_traj(frame).totalVertices=totalVertices;
    
    disp(frame)
end


% % The codes below are used for figures and visuliaztion
% % ** Note: Since in the current version, the 'getintensity' function is not
% % optimized in terms of time cost, we did not include the local intensity
% % information in the output. However, it is used in an example to
% % demonstrate the colocalization of the features and indulations at the
% % edge.
% 
%  totalVertices=boundary_traj(1).totalVertices;
%  ptotalVertices2=boundary_traj(2).totalVertices;
% 
% % colocalization figure demo 
%  figure,
%  intensity = getintensity(ptotalVertices2,ifiledir,DateTypeString,ComputerString,imgid); 
%  plotyy(1:numel(ptotalVertices2),intensity,1:numel(ptotalVertices2),ptotalVertices2(:,3))
%  hold on
%  plot(locs2, 0.05*ones(1,numel(locs2)),'ro' )
%  plot(locs1, 0.04*ones(1,numel(locs1)),'bo' )
%  hold off
% 
% % velocity v.s. curvature figure and display the pearson correlation on
% % screen. Calibration and time interval need to be adjusted for particular
% % use...
% figure,
% high_calib = 6.917; tseg = 150;
% plot(ptotalVertices2(:,3)/high_calib,ptotalVertices2(:,7)*high_calib/tseg,'bo');hold on;
% p=polyfit(ptotalVertices2(:,3)/high_calib,ptotalVertices2(:,7)*high_calib/tseg,1);
% xp = linspace(min(ptotalVertices2(:,3)/high_calib),max(ptotalVertices2(:,3)/high_calib),10);
% plot(xp,xp*p(1)+p(2),'k--','LineWidth',4);hold off
% corr(ptotalVertices2(:,3),ptotalVertices2(:,7),'type', 'pearson')
% 
% 
% % particle flow image (showing how the boundary moves)
% figure,
% for n=1:2:size(totalVertices,1)
%     [ind]=find(totalVertices2(:,6)==totalVertices(n,6));
%      if isempty(ind)~= 1
%          for m=1:max(size(ind))
%             p1 = [totalVertices(n,1) totalVertices(n,2)]';p2 = [totalVertices2(ind(m),1) totalVertices2(ind(m),2)]';
%             dp = p2 - p1;
%             quiver(p1(1),p1(2),dp(1),dp(2),0,'LineWidth',4,'Color','r','MaxHeadSize',10);
%          end
%      end
% end