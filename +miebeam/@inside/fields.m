function [ e, h ] = fields( obj, pos, bimp, k0, varargin )
%  FIELDS - Fields for infinite electron beam.
%
%  Usage for obj = miebeam.inside :
%    [ e, h ] = fields( obj, pos, bimp, k0, PropertyPairs )
%  Input
%    pos    :  field positions of size n×3
%    bimp   :  impact parameters
%    k0     :  wavenumber of light in vacuum
%  PropertyName
%    sig    :  sign in exponential of current distribution
%    imat   :  material index
%  Output
%    e,h    :  electromagnetic fields of infinite electron beam

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'sig', 1 );
addParameter( p, 'imat', 1 );
%  parse input
parse( p, varargin{ : } );

%  wavenumber of electron beam 
q = k0 / obj.vel;
%  material properties
mat = [ obj.mie.mat2, obj.mie.mat1 ];
mat = mat( p.Results.imat );
%  Lorentz contraction factor
[ eps, mu ] = deal( mat.eps( k0 ), mat.mu( k0 ) );
gamma = 1 / sqrt( 1 - eps * mu * obj.vel ^ 2 );

%  modified Bessel functions
data = cart2beam( obj, pos, bimp );
K0 = besselk( 0, q * data.r / gamma );
K1 = besselk( 1, q * data.r / gamma );
%  electric field in radial and z-direction
fac = exp( 1i * p.Results.sig * q * data.z );
er = - fac .* K1 * p.Results.sig;
ez =   fac .* K0 * 1i / gamma;
%  magnetic field in azimuthal direction
%    see Garcia de Abajo, Rev. Mod. Phys. 82, 209 (2010), Eq. (6)
ht = - fac .* K1;

%  assemble fields
e = cat( ndims( data.rad ), er .* data.rad,     ez );
h = cat( ndims( data.azi ), ht .* data.azi, 0 * ez );
%  multiply with prefactors
e = q / ( 2 * pi * obj.vel * gamma * eps ) * e;
h = q / ( 2 * pi * gamma ) * h;
%  change order of output arrays
e = permute( e, [ 1, 3, 2 ] );
h = permute( h, [ 1, 3, 2 ] );
