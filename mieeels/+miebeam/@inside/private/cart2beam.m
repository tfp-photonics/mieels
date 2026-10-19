function data = cart2beam( ~, pos, bimp, varargin )
%  CART2BEAM - Convert field positions to beam coordinates.
%
%  Usage for obj = miebeam.inside :
%    data = cart2beam( obj, pos, bimp, PropertyPairs )
%  Input
%    pos        :  field positions of size n×3
%    bimp       :  impact parameters  
%  Output
%    data.r     :  radial distance
%    data.z     :  vertical distance
%    data.rad   :  unit vector in radial direction
%    data.azi   :  unit vector in azimuthal directions


%  impact parameters
imp = bimp( : ) * [ 1, 0 ];
%  Cartesian coordinates
x = pos( :, 1 ) - imp( :, 1 ) .'; 
y = pos( :, 2 ) - imp( :, 2 ) .'; 

%  radial and vertical distance
data.r = hypot( x, y );
data.z = pos( :, 3 );
%  unit vectors
r = max( data.r, 1e-10 );
data.rad = cat( 3,   x ./ r, y ./ r );
data.azi = cat( 3, - y ./ r, x ./ r );
