function q = qinc( obj, bimp, k0, varargin )
%  QINC - Multipole coefficients for aloof electron beam
%
%  Usage for obj = miebeam.aloof :
%    q = qinc( obj, bimp, k0, PropertyPairs )
%  Input
%    bimp     :  impact parameter
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    sig      :  sign of exponential in current distribution
%  Output
%    q        :  multipole coefficients, Eq. (7)

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'sig', 1 );
%  parse input
parse( p, varargin{ : } );

%  assert aloof beams
mie = obj.mie;
assert( all( bimp > 0.5 * mie.diameter ) );

%  material properties
mat = mie.mat2;
[ eps, mu, Z ] = deal( mat.eps( k0 ), mat.mu( k0 ), mat.Z( k0 ) );
gamma = 1 ./ sqrt( 1 - eps * mu * obj.vel ^ 2 );

%  table of spherical degrees and orders
[ l, m ] = sphtable( obj.lmax );
%  beam shape coefficients
switch isempty( obj.M )
  case 0
    [ M, N ] = deal( obj.M, obj.N );
  otherwise
    [ M, N ] = beamcoeffs( obj, k0 );
end
%  deal with sign
if p.Results.sig == -1
  [ M, N ] = deal( conj( M ), conj( N ) );
end

%  modified Bessel function
[ m1, b1 ] = ndgrid( m, bimp );
K = besselk( m1, k0 * b1 / ( obj.vel * gamma ) );
%  multipole coefficients, Eqs. (87,88)
fac = 1i * mu * k0 ./ sqrt( l .* ( l + 1 ) ) .* K;
a =   fac .* N / ( obj.vel * gamma );
b = - fac .* M .* m / Z;

%  set output
q = struct( 'a', a, 'b', b, 'k0', k0, 'mode', 'aloof' );
