function plot_kymo(ifiledir)
% This function plot the kymograph of the merging of multiple blisters
% Input:
%  - ifiledir: full path of the input folder
% 
% * Note: Names of the txt files are 'kymoX.txt'. X = 1,2,3 represents the
% number of individual bundles. For example, 'kymo3.txt' contains positioal
% information when three blisters are separated. The txt files are exported
% directly from ImageJ point measurements, so column 1 = time ID, column 3
% = x, column 7 = y.
% 
% The shown code here is an example of one merging event of three blisters

% load data
addpath(ifiledir);
temp3 = load('kymo3.txt');temp2 = load('kymo2.txt');temp1 = load('kymo1.txt');

% designate time
Endtime = size(temp3,1)/3 + size(temp2,1)/2 + size(temp1,1);    % total time
Mergetime1 = size(temp3,1)/3; Mergetime2 = size(temp3,1)/3 + size(temp2,1)/2; % time when the 1st and 2nd merging evenet occurs
time = 0:Endtime-1;

% use initial coordinates to define the projection axis
ox = [temp3(1,6) temp3(2,6) temp3(3,6)];oy = [temp3(1,7) temp3(2,7) temp3(3,7)];
p = polyfit(ox,oy,1); dir = [1 p(1)]; dir =dir/norm(dir);
origin = [mean(ox) mean(oy)];
dist1 = zeros(Endtime,1);dist2 = zeros(Endtime,1);dist3 = zeros(Endtime,1);

% Chunk 1: 3 bundles 
    for k = 1:Mergetime1
            vec1 = [temp3(3*k-2,6) temp3(3*k-2,7)]-origin;
            vec2 = [temp3(3*k-1,6) temp3(3*k-1,7)]-origin;
            vec3 = [temp3(3*k,6) temp3(3*k,7)]-origin;
            dist1(k) = dot(vec1,dir);
            dist2(k) = dot(vec2,dir);
            dist3(k) = dot(vec3,dir);
    end

% Chunk 2: 2 bundles
    for k = Mergetime1+1:Mergetime2
        vec1 = [temp2(2*(k-Mergetime1)-1,6) temp2(2*(k-Mergetime1)-1,7)]-origin;
        vec2 = [temp2(2*(k-Mergetime1),6) temp2(2*(k-Mergetime1),7)]-origin;
        dist1(k) = dot(vec1,dir);
        dist2(k) = dot(vec1,dir);
        dist3(k) = dot(vec2,dir);
    end

% Chunck 3: 1 bundle (3 blisters merged into 1)
    for k =Mergetime2+1:Endtime
        vec = [temp1(k-Mergetime2,6) temp1(k-Mergetime2,7)]-origin;
        dist1(k) = dot(vec,dir);
        dist2(k) = dot(vec,dir);
        dist3(k) = dot(vec,dir);
    end
dist1 = dist1*6.917/1000;dist2 = dist2*6.917/1000;dist3 = dist3*6.917/1000;     % calibration: 6.917 micron/pixel
plot(time,smooth(dist1),'k',time,smooth(dist2),'k',time,smooth(dist3),'k');hold off

% ---------------------------------------------------
% The code below is used for similar purpose, but plot three merging events
% of two blisters. Modify the 'load data' code chunk for particular file
% names.
%
% traj = 3:1;calib = 6.917;color = ['r','k','b'];
% for i = 1:numel(traj)
%     %load data
%     bpdata = load(strcat('kymograph',num2str(traj(i)),'.txt'));
%     apdata = load(strcat('kymograph',num2str(traj(i)),' - after.txt'));
%     
%     % set time
%     Endtime = size(bpdata,1)/2 + size(apdata,1);
%     Mergetime = size(bpdata,1)/2;
%     time = 0:Endtime-1;
%     time = time/2;
%     
%     % set origin and projection direction
%     origin = [(bpdata(1,6)+ bpdata(2,6))/2 (bpdata(1,7)+ bpdata(2,7))/2];
%     dir = [bpdata(1,6) bpdata(1,7)]-origin;
%     dir =dir/norm(dir);
%     
%     % calculate distance on the projected line
%     dist1 = zeros(Endtime,1);
%     dist2 = zeros(Endtime,1);
%     for k = 1:Mergetime
%         vec1 = [bpdata(2*k-1,6) bpdata(2*k-1,7)]-origin;
%         vec2 = [bpdata(2*k,6) bpdata(2*k,7)]-origin;
%         dist1(k) = dot(vec1,dir);
%         dist2(k) = dot(vec2,dir);
%     end
%     for k =Mergetime+1:Endtime
%         vec = [apdata(k-Mergetime,6) apdata(k-Mergetime,7)]-origin;
%         dist1(k) = dot(vec,dir);
%         dist2(k) = dot(vec,dir);
%     end
%     dist1=dist1*calib/1000;dist2=dist2*calib/1000;
%     plot(time,smooth(dist1),color(i),time,smooth(dist2),color(i));hold on
%     %scatter(time,dist1-dist2)
%     %hold on
%     clear bpdata apdata time dir dist1 dist2
% end
% hold off