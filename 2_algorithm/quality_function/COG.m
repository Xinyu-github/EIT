function [varargout] = COG(A, varargin)
%COG find the center of gravity of a vector or a matrix
% [Loc1] = COG(A) finds the center of gravity of the vector A
% [Loc1, Loc2] = COG(A) finds the center of gravity of the 2-D Matrix A
% [Loc1, Loc2, Loc3,...] = COG(A) finds the center of gravity of a
% multi-dimensional arrays.
% [Loc1, ...] = COG(A, coordVec1, ...) finds the center of gravity of the
% matrix A, using input coordinates. The length of each vector coordVec
% must be equal to the size of the equivalent dimension. For example, if A
% is a 2x3 matrix, then coordVec1 is a 2-element vector, and coordVec2 is a
% 3-element vector.
%
% Note that if you use the input of the coordinates, you should include
% vectors even for singleton dimensions.
% For example, if A is a 1x9 vector, you should write:
% [Loc1] = COG(A,1,coordVec1)
% where coordVec1 is a 9-element vector.
%
% Written by Noam Greenboim
% MatlabHowTo.com
%

S = size(A);
nDim = length(S);  % number of dimensions

if isempty(varargin)
    Vflag = false;
else
    Vflag = true;  % use input coordinates
    if length(varargin) ~= nDim
        error('Input matrix dimensions and number of input coordinates must match')
    end
end

varargout = cell(1,nDim);  % preallocation

for i=1:nDim  % for every dimension
    if S(i) == 1  % check if singleton dimension
        if Vflag
            varargout{i} = varargin{i};
        else
            varargout{i} = 1;
        end
        continue
    end
    M = shiftdim(A,i-1);
    Sx = zeros(S(i),1);   % preallocation
    for j=1:S(i)
        Sx(j,1) = sum(M(j,:));
    end
    Wx = sum(Sx);   % total weight
    if Vflag  % use input coordinates
        x_loc = varargin{i};
    else
        x_loc = 1:S(i);
    end
    varargout{i} = (x_loc(:)'*Sx)/Wx;
end

if nargout<=1 &&  sum(S>1)==1  % if only one non-singleton dimension
    varargout = varargout(S>1);
end
