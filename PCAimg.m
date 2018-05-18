function PCAimg(ifiledir,ofiledir,DateTypeString,ComputerString,SquareSize,vTHRESH,frameNUM,imgNUM)

% This function segment the raw transmission images into small squares, and
% calculate the global pattern orientation of biofilm.
% 
% Input:
%  - i/ofiledir: full path of the input/output file folder;
% 
%  - DataTypeString & ComputerString: Movie ID string;
% 
%  - SquareSize: window size of the segments;
% 
%  - vTHRESH: (optional) a two element vector. vTHRESH(1) specifies the
%  threshold value of analyzing the 2D Fourier spectrum, above which
%  features/pattern can be identified; vTHRESH(2) should be a negative
%  value, which will be assigned to a segment as its local pattern
%  orientation parameter 'S', if no identifiable features are found. Note
%  that according to our definition, 'S' goes from 0 to 1 for patterned
%  region. default vector = [0.1 -0.3];
%
%  - imgNUM: (optional) image ID to define the center of biofilm;
%  default value = 25;
% 
%  - frameNUM: (optional) image ID for time course analysis; default
%  value = 5:5:290;


%% Step 1: Set default values for optional input parameters

  if ~exist('imgNUM','var')
      imgNUM = 25;
  end
  if ~exist('frameNUM','var')
      frameNUM = 5:5:290;
  end
  if ~exist('vTHRESH','var')
      vTHRESH = [0.1 -0.3];
  end
  THRESH = vTHRESH(1);
  RANDMVAL = vTHRESH(2);
  
  % This chunk of code is for outputting movies showing the evolution of
  % pattern orientation over the entire biofilm
  
      % MovieString = strcat(ComputerString, DateTypeString, '_',num2str(SquareSize),'_Order.avi');
      % waveMV = VideoWriter(MovieString);
      % waveMV.FrameRate = 2;
      % open(waveMV); 

%% Step 2: Initialization

    [c_row,c_col,sz] = FindCenter(DateTypeString,ComputerString,imgNUM);  % Find Center of the biofilm
    szNew = floor(sz/SquareSize);      % Coarse grained size
    cRowNew = floor(c_row/SquareSize)+1; cColNew = floor(c_col/SquareSize)+1;
    EigenPCA = zeros(szNew);    % To store 1st PCA eigenvalues
    OrderPara = zeros(szNew);   % To store local pattern orientation parameter 'S'
    
    % coordinate matrix within a single square segment
    CordinateMatrix1 = repmat(1:SquareSize, SquareSize,1);  % x-cordinate
    CordinateMatrix2 = CordinateMatrix1';           %y-cordinate
    
    % coarse grained coordinate matrix for all segments
    CordinateMatrix3 = repmat(1:szNew(2), szNew(1),1);  % x-cordinate
    CordinateMatrix4 = repmat([1:szNew(1)]', 1,szNew(2));           %y-cordinate
    DistanceMatrix = sqrt((double(CordinateMatrix3) - double(cColNew)).^2 + (double(CordinateMatrix4) - double(cRowNew)).^2);
    
 
    for frame = 1:length(frameNUM)
        
        f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frameNUM(frame),'%03d'),'.jpg'));    %desktop
        g = double(rgb2gray(f))/255;
 
%% Step 3: Analysis
    for RowInd = 1:szNew(1)
        for ColInd = 1:szNew(2)
        
        % 3.1: use 2D fft to determine wheter the segment has features and
        % to generate K vector 
            imgElement = g((RowInd-1)*SquareSize + 1:(RowInd)*SquareSize,(ColInd-1)*SquareSize + 1:(ColInd)*SquareSize);
            imgElement = imadjust(imgElement);
            imgfft = fft2(double(imgElement));
            
            % For high quality images, the pattern edge might be identified by
            % the built-in function 'edge'; however, this is not the case for
            % our movies...
            
%             imshow(imgElement)
%             figure, imshow(edge(imgElement,'sobel'))
            
            ampfft = mat2gray(abs(fftshift(imgfft)));
            % For visualization of what the 2d fft spectrum looks like.
%             figure,imshow(ampfft)
            
            ampdescend = sort(ampfft(:),'descend');     
            if ampdescend(2)<THRESH
                OrderPara(RowInd,ColInd) = RANDMVAL;
                continue
            end
            
            % identify the candidate direction K
            [~,LargestK] = max(ampfft(:)); ampfft(LargestK) = 0;
            [~,SecLargestK] = max(ampfft(:));   
            
            % Note that k direction with the 1st & 2nd largest fft
            % amplitude should be refelction points of each another.
            [k0y,k0x] = ind2sub(size(imgElement),LargestK(1));
            [k1y,k1x] = ind2sub(size(imgElement),SecLargestK(1));
            Kvector = [k1x-k0x, k1y-k0y];
    
    % 3.2: PCA Analysis
            PCAdata(:,3) = imgElement(:);
            PCAdata(:,1) = CordinateMatrix1(:);
            PCAdata(:,2) = CordinateMatrix2(:)';
    
            [coeff,~,~] = pca(PCAdata);
            PrincipalVector1 = coeff(1:2,1);
            PrincipalVector2 = coeff(1:2,2);
            
            if abs(dot(PrincipalVector1,Kvector)) < abs(dot(PrincipalVector2,Kvector))
                PrincipalVector = PrincipalVector1;
            else
                PrincipalVector = PrincipalVector2;
            end
            
            PosVector = double([(ColInd-0.5)*SquareSize - c_col ; (RowInd-0.5)*SquareSize - c_row]);
            OrientVal = dot(PrincipalVector, PosVector)/(norm(PrincipalVector)*norm(PosVector));
            OrderPara(RowInd,ColInd) = (OrientVal)^2;
        
        end
    end
            
%% Step 4: characterization of the pattern orientation

            % Characterization 4.1: heat map
%             imagesc(OrderPara,[RANDMVAL 1])
%             colormap pink
%             savefig(gcf, strcat(ofiledir,ComputerString,'-',DateTypeString,'-OrderPara-',num2str(frameNUM(frame)),'.fig'))
%             FRM = getframe;
%             writeVideo(waveMV,FRM);
            
            
            % Characterization 4.2: average orientation vs r
            Rmax = min(szNew(1)-cRowNew, cRowNew);
            Radius = 1:Rmax-1;
            Orient = zeros(1,length(Radius));
%             StdOrient = zeros(1,length(Radius));
            
            for q = 1:length(Radius)
                if isempty(find(DistanceMatrix > Radius(q)-0.5 & DistanceMatrix < Radius(q)+0.5 & OrderPara >= 0,1)) || numel(find(DistanceMatrix > Radius(q)-0.5 & DistanceMatrix < Radius(q)+0.5 & OrderPara >= 0)) < 0.1*numel(find(DistanceMatrix > Radius(q)-0.5 & DistanceMatrix < Radius(q)+0.5))
                    continue
                end
                Orient(q) = mean(OrderPara(find(DistanceMatrix > Radius(q)-0.5 & DistanceMatrix < Radius(q)+0.5 & OrderPara >= 0)));
%                 StdOrient(q) = std(OrderPara(find(DistanceMatrix > Radius(q)-0.5 & DistanceMatrix < Radius(q)+0.5 & OrderPara >= 0 )));
            end
            clear data
            data(:,1) = Radius*SquareSize;
            data(:,2) = Orient;
            
            save(strcat(ofiledir,'0.7datahist_',num2str(frameNUM(frame),'%03d'),'.mat'),'')
            
            % visualization:  
%             plot(Radius*SquareSize, Orient, 'b-','Linewidth',6)
%             xlabel('r (pixels)','Fontsize',16);ylabel('<cos^2\theta>','Fontsize',16)
%             ax=gca; ax.LineWidth = 2; ax.FontSize = 16; box on;

%             saveas(gcf, strcat(ofiledir,ComputerString,'-',DateTypeString,'-rdep-',num2str(frameNUM(frame)),'.fig'))
            
            
            % If you want to save raw data for each frame, and generate
            % fancier plots, you can insert code HERE.
    end
% close(waveMV)
    
    