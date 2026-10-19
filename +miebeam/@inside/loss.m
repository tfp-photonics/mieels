function P = loss( obj, bimp, k0, varargin )
%  LOSS - Energy loss probability for inside beam.
%
%  Usage for obj = miebeam.inside :
%    P = loss( obj, bimp, k0, PropertyPairs )
%  Input
%    bimp     :  impact parameters
%    k0       :  wavenumber of light in vacuum
%  Output
%    P.surf   :  surface loss
%    P.rad    :  radiative loss
%    P.bulk   :  bulk loss

%  assert inside beams
mie = obj.mie;
a = 0.5 * mie.diameter;
assert( all( bimp < a ) );

%  prefactor 
miebeam.atomicunits;
fac = 4 * fine ^ 2 / ( bohr * hartree * k0 );
%  allocate output
P = struct( 'bimp', bimp, 'k0', k0 );
[ P.surf, P.rad, P.bulk ] = deal( zeros( size( bimp ) ) );
%  bulk loss
P.bulk = 2 * sqrt( a ^ 2 - bimp .^ 2 ) * bulkloss( obj, k0 );

%  loop over impact parameters
for i1 = 1 : numel( bimp )
  %  solve Mie equations
  sol = solve( obj, bimp( i1 ), k0 ); 
  %  radiative loss
  P.rad( i1 ) = fac / k0 ^ 2 * sum( abs( sol.b2 ) .^ 2 + abs( sol.c2 ) .^ 2 );
  %  surface loss
  data = work( obj, sol );
  P.surf( i1 ) = fac * sum( data.w .* data.work );
end
