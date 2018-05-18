function imgtest(ifiledir,ofiledir,DateTypeString,ComputerString, IsOutScreen, OptStructure,vThreshold,imgNUM, frameNUM)

% This function calculate the global wavelength of the biofilm from the
% movie. It also yeilds imformation about biofilm expansion and pattern propagation.

% Input:
%  - i/ofiledir: full path of the input/output file folder;
% 
%  - DataTypeString & ComputerString: Movie ID string;
%
%  - IsOutScreen: 1 if the colony biofilm in the last frame of the movie 
%  exceeds the field of view; otherwise, this parameter is set to 0;
%
%  - OptStructure: (optional) This parameter specifies the particular
%  region for analysis, for example to exclude mutation sectors from
%  analysis. The structure variable contains (1) vThetaMin , (2) vThetaMax,
%  (together defines the range of angle), and (3) vNitv (sampling frequency
%  for 'fft').
% 
%  - vThreshold: (optional) a three element vector. vThreshold(1) 
%  defines the cutoff level in the Fourier spectrum for peak identification;
%  vThreshold(2&3) define the threshold intensities (from 0 to 1) of the 
%  center cell death zone & the herringbone pattern region. The values
%  should be adjusted to match direct counting and observation by eyes. 
%  Default value = [5.5 0.2 0.45].
%
%  - imgNUM: (optional) image ID to define the center of biofilm;
%  default value = 30;
% 
%  - frameNUM: (optional) image ID for time course analysis; default
%  value = 1:5:288;
%
%
% Output:
%  - waveNUM: wave number N_s versus R information;
%  - waveWVI: peak prominence of the identified periodic modes;
%  - DataGrowth: dynamics of colony expansion and pattern propagation.
%
%  All data is saved to 'ofiledir' as txt files. See [Step 1] for more
%  details about the output data structure.


%% Step 1: Preprocessing

% Set default values for optional input parameter.
    if ~exist('imgNUM','var')
        imgNUM = 30;
    end
    if ~exist('frameNUM','var')
        frameNUM = 1:5:288;
    end
    if ~exist('OptStructure','var')
        OptStructure = struct('vThetaMin',0,'vThetaMax',2*pi,'vNitv',4000);
    end
    if ~exist('vThreshold','var')
        vThreshold = [5.5 0.2 0.45];
    end

% Use function 'FindCenter' to calculate the center of the colony biofilm.
% It also returns the size (in pixels) of the images for further use.
    [c_row,c_col,sz] = FindCenter(ifiledir,DateTypeString,ComputerString,imgNUM);

% Set parameter values for wavelength analysis.  
    switch IsOutScreen
        case 0
            RMAX = min(c_row,sz(1)-c_row);  
            Nitv = OptStructure.vNitv;
            thetaMin = OptStructure.vThetaMin;
            thetaMax = OptStructure.vThetaMax;
            
        case 1
            RMAX = min(c_col,sz(2)-c_col);
            Nitv = 1000; 
            thetaMin = -0.6; thetaMax = 0.6;
    end
    
    PeakCriterion = vThreshold(1);
    THRESH1 = vThreshold(2); THRESH2 = vThreshold(3);
    K = 0:(Nitv/2);
    theta = linspace(thetaMin,thetaMax,Nitv);
    radius = 100:10:RMAX-100;   % Rmin = 100 pixel to avoid sigularity at
                                % the center

% Construct the output data files.
% - DataWVNUM & DataWVI
%       column 1: time t (in frame number); row 1: radius R (in pixel);
%       DataWVNUM(i,j): N_s at t(i) and R(j)
% - DataGrowth
%       column 1: time t (in frame number); column 2: biofilm radiuss R_s
%       column 3&4: winkle-delamination pattern radius R_f;
%       column 5: center cell death region radius;
%       column 6 : herringbone pattern radius;

    DataWVNUM = zeros(length(frameNUM)+1,length(radius)+1);
    DataWVNUM2 = zeros(length(frameNUM)+1,length(radius)+1);
    DataWVI = zeros(length(frameNUM)+1,length(radius)+1);
    DataGrowth = zeros(length(frameNUM),6);
    DataWVNUM(2:end,1)=frameNUM';DataWVNUM(1,2:end)=radius;
    DataWVNUM2(2:end,1)=frameNUM';DataWVNUM2(1,2:end)=radius;
    DataWVI(2:end,1)=frameNUM';DataWVI(1,2:end)=radius;
    DataGrowth(:,1)=frameNUM';                                       
    

%% Step 2: Wavelength analysis

for frame = 1:length(frameNUM)      % t loops
    
    % Use function 'meanradius' to prepocess the image and identify <R_f>.
    % The returned 'imgvar' is the processed image for wavelength analysis
    [imgvar,Rout,imgm,Tval] = meanradius(ifiledir,DateTypeString,ComputerString,frameNUM(frame),c_row,c_col,thetaMin,thetaMax);
    wavenum = zeros(1,length(radius));
    wvNUM1 = zeros(1,length(radius)); wvSTR1 = zeros(1,length(radius));
    % imshow(imgvar) 

    for j=1:length(radius)          % R loops

        circlePixels(1,:) = c_row + floor(radius(j)*sin(theta));
        circlePixels(2,:) = c_col + floor(radius(j)*cos(theta));

        imgvar1 = imgvar(:);    
        Idx = sub2ind(sz,circlePixels(1,:),circlePixels(2,:));
        valPixels = imgvar1(Idx);   % intensities along a circle

        % We then tested the periodicity of the pattern with Fourier
        % transformation, and verified it by direct counting (and in some
        % cases, autocorrelation function).
        
        % method 1: Fourier Transform
        Y = fft(valPixels-mean(valPixels));
        P2 = abs(Y/Nitv);
        P1 = P2(1:Nitv/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        if max(P1) >= PeakCriterion*max(P1(floor(0.8*length(P1)):end))
            
            [pks,locs] = findpeaks(P1,1);
            [pkMax, idMax] = max(pks);
            
            % pks(idMax) = 0;   % This line is used for identifying
            % multiple peaks. See below.
            
            if length(idMax)~= 1
               continue
            end
            
            wvNUM1(j) = K(locs(idMax))*2*pi/(thetaMax-thetaMin);
            wvSTR1(j) = pkMax;  % peak prominence meaured by peak value

            % For characterizing the second largest peak appeared in the
            % Fourier power spectrum...
            
            % if max(pks) >= PeakCriterion*max(P1(floor(0.8*length(P1)):end))
            %   [pk2Max, id2Max] = max(pks); 
            %   wvNUM2(j) = K(locs(id2Max))*2*pi/(thetaMax-thetaMin);
            %   wvSTR2(j) = pk2Max;
            % end           
        end

        % method 2: Autocorrelation
        % This method is quite unstable for most of our Movies. But it
        % could be useful for other images. Parameter values
        % adjustification are needed.
        
        % cor = xcorr(valPixels - mean(valPixels),'unbiased');
        % corx = cor((length(cor)-1)/2+1:(length(cor)-1)/2+length(theta))/var(valPixels);
        % [pks1,~] = findpeaks(corx,'MinPeakProminence',.05,'MinPeakDistance',10);
        % if length(findpeaks(corx)) > 250
        %                 continue
        % end
        % wavenum(j) = length(pks1);

    end

    DataWVNUM(frame+1,2:end) = wvNUM1;
    DataWVI(frame+1,2:end) = wvSTR1;
    % DataWVNUM2(frame+1,2:end) = wavenum;
    
    %% Step 3: Growth dynamics analysis
    
    % 3.1 - outer radius of biofilm
        RoutFilm(frame) = Rout;       

    % 3.2 - outer radius of the W-D pattern region
        [~ , wnMax] = max(wvNUM1);
        iind = find(wvNUM1==0 & radius>radius(wnMax),1,'first');
        if isempty(iind)
            RoutPattern(frame) = Rout;
        else
            RoutPattern(frame) = radius(iind);  %outer radius of wrinkle pattern
        end
        RoutPattern2(frame) = radius(wnMax);

    % 3.3 - characterizing center region
    
        RinPattern1(frame) = track_inner(imgvar,c_row,c_col,sz,theta,radius(1),radius(end),THRESH1);
        RinPattern2(frame) = track_inner(imgvar,c_row,c_col,sz,theta,radius(1),radius(end),THRESH2);

    disp(frame)

end



%%  Step 4: Export data to file

DataGrowth(:,2) = RoutFilm';
DataGrowth(:,3) = RoutPattern';DataGrowth(:,4) = RoutPattern2';
DataGrowth(:,5) = RinPattern1';DataGrowth(:,6) = RinPattern2';

% ofiledir = 'C:\Users\Andrew\Desktop\data\test\';

fwvNUM = fopen(strcat(ofiledir,ComputerString,'-',DateTypeString,'-WV.txt'),'w');
[DataRow, DataCol] = size(DataWVNUM);
for l = 1:DataRow
    for k = 1:DataCol
        fprintf(fwvNUM,'%10.8f  %',DataWVNUM(l,k));
    end
    fprintf(fwvNUM,'\n');
end
fclose(fwvNUM);

% ** If multiple peaks identification is used...
% fwvNUM2 = fopen(strcat(ofiledir,ComputerString,'-',DateTypeString,'-WV2.txt'),'w');
% [DataRow, DataCol] = size(DataWVNUM2);
% for l = 1:DataRow
%     for k = 1:DataCol
%         fprintf(fwvNUM2,'%10.8f  %',DataWVNUM2(l,k));
%     end
%     fprintf(fwvNUM2,'\n');
% end
% fclose(fwvNUM2);

fwvi = fopen(strcat(ofiledir,ComputerString,'-',DateTypeString,'-WVSTR.txt'),'w');
[DataRow, DataCol] = size(DataWVI);
for l = 1:DataRow
    for k = 1:DataCol
        fprintf(fwvi,'%10.8f  %',DataWVI(l,k));
    end
    fprintf(fwvi,'\n');
end
fclose(fwvi);

f3R = fopen(strcat(ofiledir,ComputerString,'-',DateTypeString,'-3R.txt'),'w');
[DataRow, DataCol] = size(DataGrowth);
for l = 1:DataRow
    for k = 1:DataCol
        fprintf(f3R,'%10.8f  %',DataGrowth(l,k));
    end
    fprintf(f3R,'\n');
end
fclose(f3R);
