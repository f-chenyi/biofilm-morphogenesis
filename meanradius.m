function [imgvar,Radius,imgm,Tval] = meanradius(ifiledir,DateTypeString,ComputerString,frame,c_row,c_col,thetaMin,thetaMax)

% This function calculates the mean radius of the biofilm and returns the
% processed image needed for further analysis. See 'imgtest' for details
% about the input.

% ifiledir = 'H:\Backup Image\';
%% Step 1: Image processing
%  1.1 - Enhance raw image
    f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frame,'%03d'),'.jpg'));   
    g = f(:,:,3); g = mat2gray(g); sz = size(g);
    imgvar = imadjust(g);
    imgm = medfilt2(imgvar);
    % ** Checkmark: need to refine method here! 

% 1.2 - Pick out biofilm
    [counts,x] = imhist(imgm,30);   % determine the bg cutoff intensity
    x = x(find(counts));
    counts = counts(find(counts));
    [~, Tlocs] = findpeaks(-counts);
    Tval = x(Tlocs);
    if isempty(Tval)
        error('Unexpected error when transforming bw image');
    end


    imgbw1 = im2bw(imgm,Tval(end)); % post-binarization processing
    imgfill1 = 1-imfill(1-imgbw1);
    se = strel('disk',10);
    imgc = imclose(imgfill1, se);
     % figure, imshow (imgc)
 
%% Step 2: Biofilm periphery analysis  

    [B,L] = bwboundaries(1-imgc,'noholes');
    
    % The centeroids of the identified candidates are confined in the 
    % middle of the image to exclude possible processing flaws at the edge.
    centroid = regionprops(L,'Centroid');
    c = cell2mat(struct2cell(centroid));
    c = reshape(c,[2 length(centroid)]);
    filtind = find(c(1,:) > sz(2)/3 & c(1,:)< 2*sz(2)/3);
    
    % Find the longest border, which supposedly should be the boundary of
    % the film. 
    for k = 1:length(filtind)
        BorderLength(k) = length(B{filtind(k)});
    end
    [~, LongInd] = max(BorderLength);
    Vertices = B{filtind(LongInd)};
      % imshow(g)   % doublecheck by visualization
      % hold on
      % plot(Vertices(:,2),Vertices(:,1),'b')
      % hold off
      
    Verticesnew = Vertices(Vertices(:,1)>5 & Vertices(:,2)>5& Vertices(:,1) < (sz(1)-5) & Vertices(:,2) < (sz(2)-5),:);
    thetaver = (cart2pol(Verticesnew(:,2)-c_col,Verticesnew(:,1)-c_row));
    thetaver(thetaver<0) = thetaver(thetaver<0) + 2*pi;
    Verticesnew = Verticesnew(thetaver > thetaMin & thetaver < thetaMax,:);
    Dist = sqrt((Verticesnew(:,1)-c_row).^2+(Verticesnew(:,2)-c_col).^2);
    Dist = sort(Dist); Radius = mean(Dist);
