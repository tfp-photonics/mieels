function q = qinc( obj, bimp, k0, varargin )
%  QINC1 - Multipole coefficients for inside electron beam
%
%  Usage for obj = miebeam.inside :
%    q = qinc( obj, bimp, k0, PropertyPairs )
%  Input
%    bimp     :  impact parameter
%    k0       :  wavenumber of light in vacuum
%  PropertyName
%    sig      :  sign of exponential in current distribution
%  Output
%    q        :  multipole coefficients for Mie solver

%  assert inside beams
mie = obj.mie;
assert( all( bimp < 0.5 * mie.diameter ) );

%  table of spherical degrees and orders
base = multipole.base( mie.lmax );
tab = base.tab;
%  allocate output
q = struct( 'k0', k0, 'mode', 'inside' );
[ q.e, q.h ] = deal( zeros( 2, numel( tab.l ), numel( bimp ) ) );

for i1 = 1 : numel( bimp )
  %  quadrature points
  b = bimp( i1 );
  quad = quadeels( obj, b );
  q1 = mie.qinc( @( pos, k0 ) fields( obj, pos, b, k0, 'imat', 1, varargin{ : } ), k0, 'quad', quad );
  q2 = mie.qinc( @( pos, k0 ) fields( obj, pos, b, k0, 'imat', 2, varargin{ : } ), k0, 'quad', quad );
  q.e( :, :, i1 ) = q1.e - q2.e;
  q.h( :, :, i1 ) = q1.h - q2.h;
end
