function t = transition( obj, k0 )
%  TRANSITION - Transition matrix for Mie solver.
%
%  Usage for obj = miebeam.miesolver :
%    t = transition( obj, k0 )
%  Input
%    k0       :  wavenumber of light in vacuum
%  Output
%    t        :  transition matrix for [ b1, c1, b2, c2 ]

a = 0.5 * obj.diameter;
%  material properties at particle inside and outside
[ mat1, mat2 ] = deal( obj.mat1, obj.mat2 );
[ k1, Z1 ] = deal( mat1.k( k0 ), mat1.Z( k0 ) );  x1 = a * k1;
[ k2, Z2 ] = deal( mat2.k( k0 ), mat2.Z( k0 ) );  x2 = a * k2;
%  angular degrees
ltab = multipole.sphtable( obj.lmax );

%  Riccati-Bessel functions
[ j, jp ] = multipole.riccatibessel( ltab, x1, 'j' ); 
[ h, hp ] = multipole.riccatibessel( ltab, x2, 'h' ); 
%  allocate transition matrix
t = zeros( 4, 4, numel( ltab ) );

%  continuity of transverse electric field
t( 1, 1, : ) = + j;
t( 1, 3, : ) = - h;
t( 2, 2, : ) = - 1i * jp ./ x1;
t( 2, 4, : ) = + 1i * hp ./ x2;
%  continuity of transverse magnetic field
t( 3, 2, : ) = - j / Z1;
t( 3, 4, : ) = + h / Z2;
t( 4, 1, : ) = - 1i * jp ./ ( x1 * Z1 );
t( 4, 3, : ) = + 1i * hp ./ ( x2 * Z2 );
%  invert coefficient matrices
t = pageinv( t );
