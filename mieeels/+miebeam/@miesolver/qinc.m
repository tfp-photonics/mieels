function q = qinc( obj, fun, k0, varargin )
%  QINC - Inhomogeneity for incoming fields.
%
%  Usage for obj = miebeam.miesolver :
%    q = qinc( obj, fun, k0, PropertyPairs )
%  Input
%    fun      :  incoming fields [e,h]=fun(pos,k0)
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    quad     :  quadrature points and weights
%  Output
%    q        :  inhomogeneity for incoming fields

%  table of spherical degrees and orders
base = multipole.base( obj.lmax );
tab = base.tab;

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'quad', quadsph( base, 'vector' ) );
%  parse input
parse( p, varargin{ : } );

%  vector functions
quad = p.Results.quad;
[ xm, xe ] = multipole.vsh( tab.l, tab.m, quad.t, quad.u );
xe = 1i * xe;
%  positions where fields are requested
pos = 0.5 * obj.diameter *  ...
  [ cos( quad.u) .* sin( quad.t ), sin( quad.u ) .* sin( quad.t ), cos( quad.t ) ];
[ e, h ] = fun( pos, k0 );
e = reshape( e, size( pos, 1 ), 3, [] );
h = reshape( h, size( pos, 1 ), 3, [] );

%  dummy indices for tensor class
[ l, i, m, k ] = deal( 1, 2, 3, 4 );
%  convert vector spherical harmonics to tensor class
xm = misc.tensor( xm, [ l, i, k ] );
xe = misc.tensor( xe, [ l, i, k ] );
%  convert electromagnetic fields to tensor class
e = misc.tensor( quad.w .* e, [ i, k, m ] );
h = misc.tensor( quad.w .* h, [ i, k, m ] );
%  perform integration
em = double( sum( conj( xm ) * e, i, k ), [ l, m ] );
ee = double( sum( conj( xe ) * e, i, k ), [ l, m ] );
hm = double( sum( conj( xm ) * h, i, k ), [ l, m ] );
he = double( sum( conj( xe ) * h, i, k ), [ l, m ] );

%  set output
q.k0 = k0;
q.e = permute( cat( 3, em, ee ), [ 3, 1, 2 ] );
q.h = permute( cat( 3, hm, he ), [ 3, 1, 2 ] );
