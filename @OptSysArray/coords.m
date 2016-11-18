function [X,Y] = coords(OSA,local)

% COORDS: Returns a set of vectors representing the coordinates of the
% pixels in the grid.
%
% NOTE: The result is the coordinates in the current domain.  That is, if
% the domain is 'x', the result is x and y coordinate vectors.  If the
% domain is 'k', the result is kx and ky coordinate vectors.
%
% USAGE:
% [X,Y] = coords(OptSysArray,[local])
%
% INPUTS:
% OptSysArray: The object in question.
% local: optional bool if you want to have origin suppressed.

% OUTPUTS:
% X,Y: grids giving the x|kx and y|ky grid coordinates.
%
% Written by: Johanan L. Codona, Steward Observatory: CAAO
% July 1, 2002: First test version.
% July 12, 2002: Adapted to the OptSysArray class.
% 20090407: JLCodona New version for AOSim2.
% 20090417 JLCodona: Added Offset support for segmented pupils.

if(nargin>1 && local)
    ORIGIN = OSA.origin();  % ALWAYS include the origin.
else
    ORIGIN = OSA.origin() + OSA.Offset;
end

nx = OSA.nx; ny = OSA.ny;
dx = OSA.dx; dy = OSA.dy;

%% X:  NOTE: x is dim 2.
n0 = OptSysArray.middlePixel(nx);
X = [1-n0:nx-n0];

% x-space...
%     X = X * dx + ORIGIN(2) + AOG.Offset(2);
X = X * dx + ORIGIN(2);

%% Y (dim 1)
if(nargout>1)
    n0 = midp(ny);
    Y = [1-n0:ny-n0];
    
    %        Y = Y * dy + ORIGIN(1) + AOG.Offset(1);
    Y = Y * dy + ORIGIN(1);
end

end

function p = isodd(n)
p = (mod(n,2)==1);
return;
end


function org = midp(n)
if(isodd(n))
    org = (n-1)/2+1;
else
    org = n/2+1;
end

end
