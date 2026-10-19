function P = loss( obj, bimp, k0, varargin )
%  LOSS - Energy loss probability for aloof and inside beams.
%
%  Usage for obj = miebeam.solver :
%    P = loss( obj, bimp, k0, PropertyPairs )
%  Input
%    bimp     :  impact parameters
%    k0       :  wavenumber of light in vacuum
%  Output
%    P.surf   :  surface loss
%    P.rad    :  radiative loss
%    P.bulk   :  bulk loss

%  allocate output
P = struct( 'bimp', bimp, 'k0', k0 );
[ P.surf, P.rad, P.bulk ] = deal( zeros( size( bimp ) ) );
%  index to aloof beams
ind = bimp > 0.5 * obj.diameter;

%  aloof beams
if nnz( ind )
  P1 = loss( obj.aloof, bimp( ind ), k0 );
  [ P.surf( ind ), P.rad( ind ), P.bulk( ind ) ] = deal( P1.surf, P1.rad, P1.bulk );
end
%  inside beams
if nnz( ~ind )
  P2 = loss( obj.inside, bimp( ~ind ), k0 );
  [ P.surf( ~ind ), P.rad( ~ind ), P.bulk( ~ind ) ] = deal( P2.surf, P2.rad, P2.bulk );
end
