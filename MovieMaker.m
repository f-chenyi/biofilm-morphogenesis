function MovieMaker(ifiledir,ofiledir,DateTypeString,ComputerString,frameNUM,agarinfo,straininfo)

% This function is used to assemble the SI movies.
%
% Input:
%  - i/ofiledir: full path of the input/output file folder;
%    e.g.: ifiledir = ''H:\Backup Image\'';
%
%  - DataTypeString & ComputerString: Raw movie info string;
%
%  - frameNUM: image ID for the processed movie
%
%  - (optional) agarinfo: string describing agar concentration
%    e.g.: straininfo = 'Agar conc. = 0.6%';
%
%  - (optional) straininfo: string describing the strain
%    e.g.: straininfo = 'Strain: \Delta\it{BC}';
% 
% Output: 
%  - avi movie (compressed to <10 mb) saved to 'ofiledir'

    MovieString = strcat(ofiledir,ComputerString, '_',DateTypeString, '_moviecomp.avi');
    waveMV = VideoWriter(MovieString,'Motion JPEG AVI');
    waveMV.FrameRate = 10;waveMV.Quality = 35;  % setting frame rate and movie quality
    open(waveMV);

    hour = 4;
    minute = 45; % setting initial time
    
    
for frame = 1:1:frameNUM
    
    f = imread(strcat(ifiledir,ComputerString,'\',DateTypeString,'\DSC_0',num2str(frame,'%03d'),'.jpg'));
    [counts,binlocs] = imhist(f(:,:,3));
    [pks,idx] = findpeaks(-counts);binlocs(idx(end));
    [~,mind] = max(counts);
    cutoff = find(idx < mind,1,'last');

    % pick out the largest intensity peak (which is the background)
    g = f(:,:,3);
    if ~isempty(cutoff)
        g = g(g>=binlocs(idx(cutoff)));
    else
        g = g(:);
    end
    
    % fit background intensity distribution with Gaussian
    hx = histfit(g);
    center_bg = mean(hx(2).XData);
    SD_bg  = (max(hx(2).XData) - min(hx(2).XData))/6;
    gg = f(:,:,3);
    imshow(gg,[min(min(f(:,:,3))) (center_bg - double(min(min(f(:,:,3)))))/0.9 +   double(min(min(f(:,:,3))))]) % set bg intensity = 0.9 for all frames in the movie
    saveas(gcf,strcat(num2str(frame,'%03d'),'.svg'),'svg')
    continue
    % add time 
    minute = minute + 15; 
    if minute >= 60
        hour = hour + 1;
        minute = 0;
    end
    text(5000,3800,strcat(num2str(hour,'%02d'),'h',num2str(minute,'%02d'),'min'),'Fontsize',24,'Fontweight','bold');
    
    % add annotations to movie if applicable
    if exist('agarinfo','var')
        text(100,200,agarinfo,'Fontsize',24,'Fontweight','bold');    % agar conc.
    end
    
    if exist('straininfo','var')
       text(100,400,straininfo,'Fontsize',24,'Fontweight','bold','Interpreter','tex');
    end
  

    FRM = getframe;
    writeVideo(waveMV,FRM);
end

close(waveMV)