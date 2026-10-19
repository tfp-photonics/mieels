function [ a, b, c, d ] = miecoefficients( obj, k0, varargin )
%  MIECOEFFICIENTS - Mie coefficients for sphere.
%
%  Usage :
%    [ a, b, c, d ] = miebeam.miecoefficients( obj, k0, PropertyPairs )
%  Input
%    k0         :  wavevector of light in vacuum
%  PropertyName
%    ltab       :  table of angular degrees
%    excitation :  'out' or 'in' for excitation from outside or inside
%  Output
%    a,b        :  Mie coefficients for outside field
%    c,d        :  Mie coefficients for inside field

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addOptional( p, 'ltab', [] );
addParameter( p, 'excitation', 'out' );
%  parse input
parse( p, varargin{ : } );

%  wavenumbers and impedances
[ k1, Z1 ] = deal( obj.mat1.k( k0 ), obj.mat1.Z( k0 ) );
[ k2, Z2 ] = deal( obj.mat2.k( k0 ), obj.mat2.Z( k0 ) );

[ z1, z2 ] = deal( 0.5 * k1 * obj.diameter, 0.5 * k2 * obj.diameter );
%  Riccati-Bessel functions
%  compute Riccati-Bessel functions
[ j1, zjp1 ] = multipole.riccatibessel( 1 : obj.lmax, z1, 'j' );
[ j2, zjp2 ] = multipole.riccatibessel( 1 : obj.lmax, z2, 'j' );
[ h1, zhp1 ] = multipole.riccatibessel( 1 : obj.lmax, z1, 'h' );
[ h2, zhp2 ] = multipole.riccatibessel( 1 : obj.lmax, z2, 'h' );

switch p.Results.excitation
  case 'out'
    %  Mie coefficients for outside field, Hohenester Eq. (12)
    a = ( Z2 * z1 * j1 .* zjp2 - Z1 * z2 * j2 .* zjp1 ) ./  ...
        ( Z2 * z1 * j1 .* zhp2 - Z1 * z2 * h2 .* zjp1 );
    b = ( Z2 * z2 * j2 .* zjp1 - Z1 * z1 * j1 .* zjp2 ) ./  ...
        ( Z2 * z2 * h2 .* zjp1 - Z1 * z1 * j1 .* zhp2 );
    %  Mie coefficients for inside field, Eq. (E.13)
    c = ( Z2 * z2 * j2 .* zhp2 - Z2 * z2 * h2 .* zjp2 ) ./  ...
        ( Z2 * z1 * j1 .* zhp2 - Z1 * z2 * h2 .* zjp1 ) * k1 / k2;
    d = ( Z2 * z2 * h2 .* zjp2 - Z2 * z2 * j2 .* zhp2 ) ./  ...
        ( Z2 * z2 * h2 .* zjp1 - Z1 * z1 * j1 .* zhp2 ) * k1 / k2;
  case 'in'
    %  Mie coefficients for outside field
    a = ( Z1 * z1 .* h1 .* zjp1 - Z1 * z1 .* j1 .* zhp1 ) ./  ...
        ( Z2 * z1 .* j1 .* zhp2 - Z1 * z2 .* h2 .* zjp1 ) * k2 / k1;
    b = ( Z1 * z1 .* j1 .* zhp1 - Z1 * z1 .* h1 .* zjp1 ) ./  ...
        ( Z2 * z2 .* h2 .* zjp1 - Z1 * z1 .* j1 .* zhp2 ) * k2 / k1;
    %  Mie coefficients for inside field
    c = ( Z1 * z2 .* h2 .* zhp1 - Z2 * z1 .* h1 .* zhp2 ) ./  ...
        ( Z2 * z1 .* j1 .* zhp2 - Z1 * z2 .* h2 .* zjp1 );
    d = ( Z1 * z1 .* h1 .* zhp2 - Z2 * z2 .* h2 .* zhp1 ) ./  ...
        ( Z2 * z2 .* h2 .* zjp1 - Z1 * z1 .* j1 .* zhp2 ); 
end

if ~isempty( p.Results.ltab )
  [ a, b ] = deal( a( p.Results.ltab ), b( p.Results.ltab ) );
  [ c, d ] = deal( c( p.Results.ltab ), d( p.Results.ltab ) );
end

