function [ e, h ] = fields( obj, sol, pos )
%  FIELDS - Scattered electromagnetic fields.
%
%  Usage for obj = miebeam.miesolver :
%    [ e, h ] = fields( obj, sol, pos )
%  Input
%    sol      :  Mie solution
%    pos      :  positions where scattered fields are requested
%  Output
%    e,h      :  scattered electromagnetic fields

%  allocate output
m = size( sol.b1, 2 );
[ e, h ] = deal( zeros( [ size( pos ), m ] ) );
%  indices to insides and outside positions
[ ~, ~, r ] = cart2sph( pos( :, 1 ), pos( :, 2 ), pos( :, 3 ) );
ind = r < 0.5 * obj.diameter;

%  wavenumber of light in vacuum
k0 = sol.k0;
%  dummy indices for tensor class
[ i, k, l, m ] = deal( 1, 2, 3, 4 );

%  inner fields
if nnz( ind )
  %  transverse vector functions
  k1 = obj.mat1.k( k0 );
  [ M, N ] = transvec( obj, pos( ind, : ), k1, 'fun', 'j' );
  [ M, N ] = deal( misc.tensor( M, [ l, i, k ] ), misc.tensor( N, [ l, i, k ] ) );
  %  field coefficients
  b1 = misc.tensor( sol.b1, [ l, m ] );
  c1 = misc.tensor( sol.c1, [ l, m ] );
  %  electromagnetic fields
  e1 =           sum( b1 * M + c1 * N, l );
  h1 = 1i * k1 * sum( b1 * N - c1 * M, l );
  %  update fields
  e( ind, :, : ) = double( e1, [ i, k, m ] );
  h( ind, :, : ) = double( h1, [ i, k, m ] );
end

%  outer fields
if nnz( ~ind )
  %  transverse vector functions
  k2 = obj.mat2.k( k0 );
  [ M, N ] = transvec( obj, pos( ~ind, : ), k2, 'fun', 'h' );
  [ M, N ] = deal( misc.tensor( M, [ l, i, k ] ), misc.tensor( N, [ l, i, k ] ) );
  %  field coefficients
  b2 = misc.tensor( sol.b2, [ l, m ] );
  c2 = misc.tensor( sol.c2, [ l, m ] );
  %  electromagnetic fields
  e2 =           sum( b2 * M + c2 * N, l );
  h2 = 1i * k2 * sum( b2 * N - c2 * M, l );
  %  update fields
  e( ~ind, :, : ) = double( e2, [ i, k, m ] );
  h( ~ind, :, : ) = double( h2, [ i, k, m ] );
end
