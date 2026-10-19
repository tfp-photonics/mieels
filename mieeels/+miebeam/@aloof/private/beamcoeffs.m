function [ M, N ] = beamcoeffs( obj, k0 )
%  BEAMCOEFFS - Beam coefficients M,N of Stamatopoulos.
%
%  Usage for obj = miebeam.aloof :
%    [ M, N ] = beamcoeffs( obj, k0, PropertyPairs )
%  Input
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    imat   :  index to embedding medium
%    sig    :  sign of exponential in current distribution
%  Output
%    M,N    :  beam coefficients of Stamatopoulos et al., JOSA B 42, 1620 (2025)

%  material at outside or inside
mie = obj.mie;
mat = mie.mat2;
%  Lorentz contacion factor
vel = obj.vel;
[ eps1, mu1, n1 ] = deal( mat.eps( k0 ), mat.mu( k0 ), mat.n( k0 ) );
gamma1 = 1 ./ sqrt( 1 - eps1 * mu1 * vel ^ 2 );

%  tables for n! and n!!
n = 2 * obj.lmax + 1;
tab1 = cumprod( [ 1, 1 : n ] ) .';
tab2 = ( 0 : n ) .';
i1 = rem( tab2, 2 ) == 0;
tab2(  i1 ) = arrayfun( @( n ) prod( 2 : 2 : n ), tab2(  i1 ), 'uniform', 1 );
tab2( ~i1 ) = arrayfun( @( n ) prod( 1 : 2 : n ), tab2( ~i1 ), 'uniform', 1 );
%  factorial and double factorial functions
fac1 = @( n ) tab1( n + 1 );
fac2 = @( n ) tab2( max( 1, n + 1 ) );

%  table of spherical degrees and orders
[ l, m ] = multipole.sphtable( obj.lmax );
m2 = abs( m );
%  beam coefficient M, Eq. (29)
fac = sqrt( ( 2 * l + 1 ) / pi ) .* sqrt( fac1( l - m2 ) ./   ...
  fac1( l + m2 ) ) .* fac2( 2 * m2 - 1 ) ./ ( n1 * gamma1 * vel ) .^ m2;
M = fac .* gegenbauerC( l - m2, m2 + 0.5, 1 / ( n1 * vel ) );
%  deal with negative degrees
ind = m < 0;
M( ind ) = ( -1 ) .^ m( ind ) .* M( ind );

%  use expression after Eq. (43) for evaluation of N
[ Li1, i1 ] = ismember( [ l, m + 1 ], [ l, m ], 'rows' );  i1 = i1( Li1 );
[ Li2, i2 ] = ismember( [ l, m - 1 ], [ l, m ], 'rows' );  i2 = i2( Li2 );
%  coefficients of Eq. (38)
c1 = 0.5 * sqrt( ( l( Li1 ) - m( Li1 ) ) .* ( l( Li1 ) + m( Li1 ) + 1 ) );
c2 = 0.5 * sqrt( ( l( Li2 ) + m( Li2 ) ) .* ( l( Li1 ) - m( Li2 ) + 1 ) );
%  beam coefficient N
N = zeros( size( M ) );
N( Li1 ) =   c1 .* M( i1 );
N( Li2 ) = - c2 .* M( i2 ) + N( Li2 );
