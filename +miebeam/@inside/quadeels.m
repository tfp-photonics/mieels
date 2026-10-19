function quad = quadeels( obj, b )
%  QUADEELS - Quadrature point for EELS.
%
%  Usage for obj = miebeam.inside :
%    quad = quadeels( obj, b )
%  INPUT
%    b      :  impact parameter
%  Output
%    quad   :  quadrature angles and weights

%  diameter and impact parameter
diameter = obj.mie.diameter;
assert( isscalar( b ) );

if b == 0
  %  quadrature points on sphere
  [ pos, w ] = cap( diameter, pi, obj.nquad .* [ 2, 1 ] );
else
  %  impact angle
  a = 0.5 * diameter;
  t0 = 0.5 * pi - acos( b / a );
  %  integration points
  [ t, w1 ] = misc.lgwt( obj.nquad( 1 ), 0.5 * pi - t0, 0.5 * pi + t0 );
  [ u, w2 ] = misc.lgwt( obj.nquad( end ), 0, 1 );
  %  adapt azimuthal range
  u1 =     pi - asin( cot( t0 ) * cot( t ) );
  u2 = 2 * pi + asin( cot( t0 ) * cot( t ) );
  u = u1 + ( u2 - u1 ) * u .';
  w = ( u2 - u1 ) .* sin( t ) .* w1 * w2 .';
  t = repmat( t, 1, size( u, 2 ) );
  %  integration positions
  [ t, u ] = deal( t( : ), u( : ) );
  pos = 0.5 * diameter * [ cos( u ) .* sin( t ), sin( u ) .* sin( t ), cos( t ) ];
  %  add quadrature points on cap
  [ pos1, w1 ] = cap( diameter, 0.5 * pi - t0, obj.nquad );
  [ pos, w ] = deal( [ pos1; pos ], [ w1; w( : ) ] );
  %  rotate positions
  [ sint, cost ] = deal( sin( t0 ), cos( t0 ) );
  rot = [ 1, 0, 0; 0, cost, sint; 0, -sint, cost ];
  pos = pos * rot .';
  %  upper and lower halfsphere
  pos = [ pos; pos .* [ 1, 1, -1 ] ];
  w = [ w; w ];
  %  rotate positions
  pos = pos( :, [ 2, 1, 3 ] );
end

%  set output
[ quad.u, quad.t ] = cart2sph( pos( :, 1 ), pos( :, 2 ), pos( :, 3 ) );
quad.t = 0.5 * pi - quad.t;
quad.w = w;


function [ pos, w ] = cap( diameter, t0, nquad )
%  CAP - Discretize sphere cap.

%  integration points
[ t, w1 ] = misc.lgwt( nquad( 1 ), 0, t0 );
[ u, w2 ] = misc.lgwt( nquad( end ), 0, 2 * pi );
%  grid of integration points
w = reshape( sin( t ) .* w1 * w2 .', [], 1 );
%  integration positions
[ t, u ] = ndgrid( t, u );
[ t, u ] = deal( t( : ), u( : ) );
pos = 0.5 * diameter * [ cos( u ) .* sin( t ), sin( u ) .* sin( t ), cos( t ) ];
