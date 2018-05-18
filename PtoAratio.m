function PtoAratio(ifiledir,DateTypeString,ComputerString,frameNUM,analysismode)

% This function calculates the asphericity parameter alpha of the biofilm
% contour, defined as (Perimeter)^2/(4*pi*Area). alpha equals to 1 for
% perfect circle.
% 
% Input:
%  - ifiledir: full path of the input file folder;
%    e.g.: ifiledir = 'H:\Backup Image\';
% 
%  - DataTypeString & ComputerString: Movie ID string;
%
%  - frameNUM: (optional) image ID for time course analysis; default
%  value = 1:5:288;
% 
%  - analysis mode: 0 or 1, indicating the biofilm with(1)/without(0)
%  mutations. Default value = 0.
%
% * Note: For informative analysis, the current code requires that at least
% half of the edge contains no mutation. Specific codes used to cut out the
% available half (see comments below starting with ***) should be
% ajusted according to the position (top, bottom, left,or right of the
% biofilm).

% Set input values
    if ~exist('analysismode','var')
        analysismode = 0;
    end

for frame = 1:length(frameNUM)

  f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frameNUM(frame),'%03d'),'.jpg'));
  g = f(:,:,3);g = mat2gray(g);
  imgvar = imadjust(g);
  
  [c_row,c_col,sz] = FindCenter(DateTypeString,ComputerString,30);
  
  imgm = medfilt2(imgvar); 
  [counts,x] = imhist(imgm,30);
  x = x(find(counts));
  counts = counts(find(counts));
  [~, Tlocs] = findpeaks(-counts);
  Tval = x(Tlocs);
 
  imgbw = im2bw(imgm,Tval(end));  
  imgfill = imfill(1-imgbw);     
  imgc1= imclose(imgfill,strel('disk',10));
  
  switch analysismode
      case 0
          % When no mutation appears (which is the case for low agar
          % conc.), it's still possible the expanding biofilm exceeds the
          % field of view at later stage. Under this condition, we extract
          % the part within the field of view, and compute asphericity
          % based on the cut.
          
          % deal with the upper and bottom boundary if the biofilm is out
          % of field.
          if ~isempty(imgc1(2,2000:4000) ~=0 )
              idstart1 = find(imgc1(2,1000:5000) ~=0,1,'first');
              idend1 = find(imgc1(2,1000:5000) ~=0,1,'last');
              imgc1(1,idstart1:idend1)=1;
          end
          if ~isempty(imgc1(sz(1)-1,2000:4000) ~=0)
              idstart2 = find(imgc1(sz(1)-1,1000:5000) ~=0,1,'first');
              idend2 = find(imgc1(sz(1)-1,1000:5000) ~=0,1,'last');
              imgc1(sz(1),idstart2:idend2)=1;
          end

          imgc1(1,:)=0; imgc1(sz(1),:)=0;
          imgc1 = imfill(imgc1);

%           imshow(imgc1)

          % pick out the colony biofilm from all possible connected
          % components of a binary image.
          [L, ~] = bwlabel(imgc1);
          Area = cell2mat(struct2cell(regionprops(L,'Area')));
          Perimeter = cell2mat(struct2cell(regionprops(L,'Perimeter')));
          s = regionprops(L,'centroid');
          Centroids = cat(1,s.Centroid);
          Area = Area(Centroids(:,1) > 0.4*sz(2) & Centroids(:,1) < 0.6*sz(2));
          Perimeter = Perimeter(Centroids(:,1) > 0.4*sz(2) & Centroids(:,1) < 0.6*sz(2));
          [~,Ind] = max(Area);

          % decide how much of the circle should be cut
          theta_cut = 0; area_cut =0; peri_cut = 0;

          if ~isempty(idstart1) % exceeds the top boundary
              area_cut = area_cut + (idend1 - idstart1 +1)*c_row/2;
              peri_cut = peri_cut + (idend1 - idstart1 +1);
              [theta_vec1,~] = cart2pol(idstart1 - c_col, c_row - 1);
              [theta_vec2,~] = cart2pol(idend1 - c_col, c_row - 1);
              theta_cut = theta_cut + abs(theta_vec1 - theta_vec2);
          end

        if ~isempty(idstart2)   % exceeds the bottom boundary
          area_cut = area_cut + (idend2 - idstart2 +1)*(sz(1)-c_row)/2;
          peri_cut = peri_cut + (idend2 - idstart2 +1);
          [theta_vec1,~] = cart2pol(idstart2 - c_col, c_row - sz(1));
          [theta_vec2,~] = cart2pol(idend2 - c_col, c_row - sz(1));
          theta_cut = theta_cut + abs(theta_vec1 - theta_vec2);
        end

        theta = 2*pi - theta_cut;
        area = Area(Ind) - area_cut;
        perim = Perimeter(Ind) - peri_cut;

        PAratio(frame) = (perim*2*pi/theta)^2/(4*pi*(area*2*pi/theta));
  
  
    case 1
        % When mutation happens (which is usually the case when agar conc.
        % is high), but at least half of the biofilm is good, then we
        % calculate the asphericity based on this usable half.
        %
        % The idea is mannually cut the biofilm into two parts (line a),and
        % then calculate the cut length across the biofilm, which should be
        % deducted from the perimeter later. After this processing, we pick
        % out the connected regions from the given half of the image (part
        % c) and identify the biofilm region. Calculation of asphericity is
        % then conducted on this half-flower.
        % 
        % *** Below is an example to use upper half of the biofilm. For
        % other cases code line with comment '% part x' should be adjusted
        % correspondingly.
       
        imgc1(c_row-1:c_row+3,:)=0;     % part a
        imgc1(1,:)=0; imgc1(sz(1),:)=0;
        imgc1 = imfill(imgc1);
        % imshow(imgc1)
        
        peri_cut = 0;
        temp = imgc1(c_row,:);          % part b
        temp = imopen(temp,strel('disk',50));
        peri_cut = length(find(temp == 1));     
        
        [L, ~] = bwlabel(imgc1);
        Area = cell2mat(struct2cell(regionprops(L,'Area')));
        Perimeter = cell2mat(struct2cell(regionprops(L,'Perimeter')));
        s = regionprops(L,'centroid');
        Centroids = cat(1,s.Centroid);
        Area = Area(Centroids(:,2) < c_row);            % part c
        Perimeter = Perimeter(Centroids(:,2) < c_row);  % part c
        [~,Ind] = max(Area);
        
        PAratio(frame) = ((Perimeter(Ind)-peri_cut)*2)^2/(4*pi*(2*Area(Ind)));
  end
  
  clear Ind Area Perimeter
  disp(frame)
end

% plot(frameNUM,PAratio,'b-','LineWidth',2); hold off