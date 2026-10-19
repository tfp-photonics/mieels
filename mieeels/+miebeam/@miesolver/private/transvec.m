function [ M, N ] = transvec( obj, pos, k1, varargin )
%  TRANSVEC - Transverse vector functions for Helmholtz equation.
%
%  Usage for obj = miebeam.miesolver :
%    [ M, N ] = transvec( obj, pos, k0, PropertyPairs )
%  Input
%    pos    :  positions
%    k1     :  wavenumber
%  PropertyName
%    fun    :  'j' for Bessel or 'h' for Hankel function
%  Output
%    M,N    :  transverse vector function

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'fun', 'j' );
%  parse input
parse( p, varargin{ : } );

%  convert positions to spherical coordinates
[ u, t, r ] = cart2sph( pos( :, 1 ), pos( :, 2 ), pos( :, 3 ) );
t = 0.5 * pi - t;
%  spherical degrees and orders
base = multipole.base( obj.lmax );
[ ltab, mtab ] = deal( base.tab.l, base.tab.m );
%  Bessel function and vector spherical harmonics
x = r * k1;
x( abs( x ) < 1e-10 ) = 1e-10;
[ z, zp ] = multipole.riccatibessel( ltab, x, p.Results.fun ); 
[ xm, xe, x0 ] = multipole.vsh( ltab, mtab, t, u );
xe = 1i * xe;

%  dummy indices for tensor class
[ i, k, m ] = deal( 1, 2, 3 );
%  convert arrays to tensor class
[ z, zp ] = deal( misc.tensor( z, [ m, i ] ), misc.tensor( zp, [ m, i ] ) );
[ xm, xe, x0 ] = deal( misc.tensor( xm, [ m, i, k ] ),  ...
                       misc.tensor( xe, [ m, i, k ] ), misc.tensor( x0, [ m, i, k ] ) );
%  prefactor
fac = misc.tensor( sqrt( ltab .* ( ltab + 1 ) ), m );
x = misc.tensor( x, i );
%  expansion functions
M = z * xm;
N = ( fac * z * x0 - 1i * zp * xe ) ./ x;

%  set output
M = double( M, [ m, i, k ] );
N = double( N, [ m, i, k ] );
