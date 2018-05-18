function rin = track_inner(g,c_row,c_col,sz,theta,RMin,RMax,THRESH)

% This function calculates the radius of inner region using intensity
% values of the image as criterion. The function is called in 'imgtest'.
% 
% Since we used transimitted light images for analysis, we expect regions 
% with larger pattern density and height (e.g., cell death zone) to be 
% darker, i.e. with smaller gray scale intensities.
% 
% Input: 
%  - g: processed image data represented by a numeric matrix of 0 - 1
%  intensity values;
%  - c_row, c_col: cordinates of the biofilm center;
%  - sz: matrix size of g;
%  - theta: polar angle of the analyzed region;
%  - RMin, RMax: range of radius for analysis;
%  - THRESH: cutoff intensity level.
% 
% Output:
%  - rin: outer radius of the region determined by 'THRESH'.


    if ~exist('THRESH','var')
        THRESH = 0.3;
        display('Warning: Input argument THRESH is null. Default value 0.3 is used.')
    end
    
    Nr = (RMax-RMin)+1;
    radius = linspace(RMin,RMax,Nr);
    colInd = c_col + floor(radius'*cos(theta));
    rowInd = c_row + floor(radius'*sin(theta));
       
    pxvalg = g(:);
    IndTotal = sub2ind(sz,rowInd,colInd);
    PxValTotal = pxvalg(IndTotal);
    PxValmean = mean(PxValTotal,2);
    
    rin = radius(find(smooth(PxValmean,ceil(Nr/10))<THRESH,1,'last')); % 'smooth' is used for better performance
%   figure, imshow(g);hold on;viscircles([c_col c_row],rin);hold off


    
    
    