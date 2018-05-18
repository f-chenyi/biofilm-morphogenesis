function ave_theta = mean_angle(theta,ave_method)

% This function calculates the average directions of each row of the input
% matrix (specified by angles).
% 
% Input:
%  - theta: numeric matrix;
%  - ave_method: 0 or 1; This function provides two averaging options. 0
%  calculate the arithmetic mean of angle values. 1 calculates the argument
%  of mean(e^(i*theta)).


switch ave_method
    case 0
        sz = size(theta);
        for row = 1:sz(1)
            thetaold = theta(row,1);
            for col = 2:sz(2)
                thetanew = theta(row,col);
                thetaold = angle(exp(1i*thetaold)+exp(1i*thetanew))*2;
            end
            ave_theta(row,1)=thetaold/sz(2);
        end
    case 1
        cnum = exp(1i*theta);
        ave_cnum = mean(cnum,2);
        ave_theta = angle(ave_cnum);
end