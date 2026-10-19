function P = loss( obj, bimp, k0, varargin )
%  LOSS - Energy loss probability for aloof geometry.
%
%  Usage for obj = miebeam.aloof :
%    P = loss( obj, bimp, k0 )
%  Input
%    bimp     :  impact parameters
%    k0       :  wavenumber of light in vacuum
%  Output
%    P.surf   :  EELS probability, Eq. (72)
%    P.rad    :  photon emission probability, Eq. (77)

%  assert aloof beams
mie = obj.mie;
assert( all( bimp > 0.5 * mie.diameter ) );

%  spherical harmonics
[ l, m ] = multipole.sphtable( obj.lmax );
%  material properties
mat = mie.mat2;
[ eps, mu, n ] = deal( mat.eps( k0 ), mat.mu( k0 ), mat.n( k0 ) );
gamma = 1 ./ sqrt( 1 - eps * mu * obj.vel ^ 2 );
%  modified Bessel function
q = k0 / obj.vel;
[ m1, b1 ] = ndgrid( m, bimp );
K = besselk( m1, q * b1 / gamma );
%  beam coefficients M,N of Stamatopoulos
switch isempty( obj.M )
  case 0
    [ M, N ] = deal( obj.M, obj.N );
  otherwise
    [ M, N ] = beamcoeffs( obj, k0 );
end
N = N / ( n * obj.vel * gamma );
%  Mie coefficients
[ a, b ] = miebeam.miecoefficients( mie, k0, 'ltab', l );

%  atomic units
miebeam.atomicunits;
%  EELS probability, Eq. (71), convert to units of (1/eV)
fac = 4 * fine ^ 2 / ( bohr * hartree * k0 );
surf = fac * sum( abs( K ) .^ 2 ./ ( l .* ( l + 1 ) ) .*  ...
       ( abs( m .* M ) .^ 2 .* real( b ) + abs( N ) .^ 2 .* real( a ) ), 1 );
%  photon loss probability
rad = fac * sum( abs( K ) .^ 2 ./ ( l .* ( l + 1 ) ) .*  ...
                       ( abs( m .* M .* b ) .^ 2 + abs( N .* a ) .^ 2 ), 1 );

%  set output
P = struct( 'surf', surf, 'rad', rad, 'bulk', 0 * rad, 'bimp', bimp, 'k0', k0 );


