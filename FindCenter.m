function [c_row,c_col,sz] = FindCenter(ifiledir,DateTypeString,ComputerString,imgNUM)

% This function calculates the center of biofilms. See 'imgtest' for
% details about the input. See 'meanradius' for details about the image
% processing.

% ifiledir = 'H:\Backup Image\';

    f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(imgNUM,'%03d'),'.jpg'));    %desktop
    g = f(:,:,3); g = mat2gray(g); sz = size(g);
    imgvar = imadjust(g);
    
    imgm = medfilt2(imgvar);
    [counts,x] = imhist(imgm,30);
    x = x(find(counts));
    counts = counts(find(counts));
    [~, Tlocs] = findpeaks(-counts);
    Tval = x(Tlocs);
    
    % level = graythresh(imgm); % compare with matlab 'graythresh' function
    % just for interst...
    
    imgbw = im2bw(imgm,Tval(end));  
    imgfill = 1-imfill(1-imgbw);
    
    % Center is the point within biofilm region, which has the largest 
    % distance to the biofilm-background borderline.
    D = bwdist(imgfill);
    [~,I] = max(D(:));
    [c_row, c_col] = ind2sub(sz,I);