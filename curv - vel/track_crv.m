function [sz,totalVertices] = track_crv(ifiledir,DateTypeString,ComputerString,frame,resolution)

% This is a function to track local curvature. It is called by 'track_bd'.
% See input details there. 
% Need citation: the core function is from online source
% Note after the trakcing, the boudnary points' coordinated is smaller by a
% factor of 1/scale; the curvature is larger by a factor of 1/scale. 

% get binarized image (same usage in 'meanradius' and other functions)
    f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frame,'%03d'),'.jpg'));  
    g = f(:,:,3);g = mat2gray(g);
    imgvar = imadjust(g);

    imgm = medfilt2(imgvar);
    [counts,x] = imhist(imgm,30);
    x = x(find(counts));
    counts = counts(find(counts));
    [~, Tlocs] = findpeaks(-counts);
    Tval = x(Tlocs);
    if isempty(Tval)
        error('Unexpected error when transforming bw image');
    end

    imgbw1 = im2bw(imgm,Tval(end));  
    imgfill1 = 1-imfill(1-imgbw1);     
    se = strel('disk',10);
    imgc = imclose(imgfill1, se);


% use binary image to extract biofilm contour (same usage in 'meanradius')
     sz = size(g);
    [B,L] = bwboundaries(1-imgc,'noholes'); 
    centroid = regionprops(L,'Centroid');
    c = cell2mat(struct2cell(centroid));
    c = reshape(c,[2 length(centroid)]);
    
    filtind = find(c(1,:) > sz(2)/3 & c(1,:)< 2*sz(2)/3);
    
    for k = 1:length(filtind)
        BorderLength(k) = length(B{filtind(k)});
    end
    [~, LongInd] = max(BorderLength);
    Vertices = B{filtind(LongInd)};
    
% smooth the vertice. default resolution is 30, which actually means
% that I lose resolution. This is becuase I only care about global
% information. Consistently, prefer to use global moving function. 
    kappa=[];
    kappa(:,1) = smooth(Vertices(1:2:end,1));
    kappa(:,2) = smooth(Vertices(1:2:end,2));
    
%     kappa(:,1)=smooth(Vertices(1:2:end,1),'sgolay',2);
%     kappa(:,2)=smooth(Vertices(1:2:end,2),'sgolay',2);
    clear Vertices
    Vertices(:,[1,2]) = kappa(1:resolution:end, [2,1]);
    totalVertices(:,1) = spline(1:resolution:size(kappa,1),Vertices(:,1),1:size(kappa,1));
    totalVertices(:,2) = spline(1:resolution:size(kappa,1),Vertices(:,2),1:size(kappa,1)); %coarse-grained vertices

    totalVertices(:,3)=LineCurvature2D(totalVertices); % local curvature
    totalVertices(:,[4,5])=LineNormals2D(totalVertices); % prientation of the local curvature
%       imshow(g)
%       hold on
%       plot(totalVertices(:,1),totalVertices(:,2),'b-','Linewidth',4)
%       plot([totalVertices(:,1) totalVertices(:,1)+10000*totalVertices(:,3).*totalVertices(:,4)]',[totalVertices(:,2) totalVertices(:,2)+10000*totalVertices(:,3).*totalVertices(:,5)]','c','Linewidth',4);
%       hold off
