function [basis] = SDS_BSpline_Basis(B,hx,full_res)

% Author:  Samir Sharma
% Created: June 2011
% Last Modified: December 18th, 2011

% Inputs:
%   B:  1D B-Spline (can be 1D for two dimensions)
%   hx: knot size (in voxels) (can have two elements, one for each
%   dimension)
%   full_res: full resolution of the image (in voxels)

% Outputs:
%   basis: the B-Spline set

for bb = 1:numel(B)
    if B{bb} == 1
        basis{bb} = ones(1,full_res(bb));
    else
        %* BEGIN: Preliminaries *
        extent = 3*hx(bb) - 1 + full_res(bb);    % full extent which will be cropped
        num_vecs = floor((((extent - 1) - 1) / hx(bb)) + 1);
        last_starting = (num_vecs - 1)*hx(bb) + 1;
        support_size = numel(B{bb});
        tmp_basis = zeros(num_vecs, last_starting + support_size - 1);
        %* END: Preliminaries *


        %* BEGIN: Create basis *
        for aa = 1:num_vecs
            tmp_basis(aa,(aa-1)*hx(bb)+1:(aa-1)*hx(bb)+support_size) = B{bb};
        end
        tmp_basis = tmp_basis(:,3*hx(bb)+1:3*hx(bb)+full_res(bb));

        % get rid of zero vectors (this should only happen in the linear case)
        indx = [];
        for aa = 1:size(tmp_basis,1)
            if sum(abs(tmp_basis(aa,:))) == 0
                indx = [indx aa];
            end
        end
        tmp_basis(indx,:) = [];
        %* END: Create basis *
        basis{bb} = tmp_basis; clear tmp_basis;
    end
end    