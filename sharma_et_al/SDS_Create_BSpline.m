function [B,hx] = SDS_Create_BSpline(support_size,type)

% Author:  Samir Sharma
% Created: June 2011
% Last Modified: December 18th, 2011

% Inputs:
%   support_size:   1D support size (in voxels)
%                   (can be a two-element vector, one for each dimension)
%   type:           'linear' or 'cubic'

% Outputs:
%   B:              1D B-Spline for each dimension
%   hx:             knot size (in voxels)


for aa = 1:numel(support_size)
    if support_size(aa) == 1
        B{aa} = 1;
        hx(aa) = 0;
    else
        %* BEGIN: Preliminaries *
        if strcmp(type,'linear')
            hx(aa) = round(1/2*(support_size(aa) - 1));
        elseif strcmp(type,'cubic')
            hx(aa) = round(1/4*(support_size(aa) - 1));
        end
        t1 = linspace(0,1,hx(aa)+1);
        t2 = linspace(1,2,hx(aa)+1);
        t2(1) = [];
        %* END: Preliminaries *

        %* BEGIN: Create B-Spline *
        if strcmp(type,'linear')
            x1 = t1;  
            tmp_B = [x1 fliplr(x1(1:end-1))];
        elseif strcmp(type,'cubic')
            x1 = 2/3 - (1 - abs(t1)./2) .* (t1.^2);
            x2 = 1/6 * (2-abs(t2)).^3;
            x = [x1 x2];                    % concatenate the two pieces
            tmp_B = [fliplr(x(2:end)) x];   % concatenate flipped and non-flipped 
        end
        B{aa} = tmp_B; clear tmp_B;
        %* END: Create B-Spline *
    end
end
