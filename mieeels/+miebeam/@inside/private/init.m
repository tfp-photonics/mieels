function obj = init( obj, mat1, mat2, diameter, vel, varargin )
%  INIT - Initialize Mie solver for EELS and inside beams.
%
%  Usage for obj = miebeam.inside :
%    obj = init( obj, mat1, mat2, diameter, vel, PropertyPairs )
%  Input
%    mat1       :  material at sphere inside
%    mat2       :  material at sphere outside
%    diameter   :  diameter of sphere
%    vel        :  electron velocity in units of speed of light
%  PropertyName
%    lmax       :  cutoff for angular degrees
%    nquad      :  number of quadrature points for penetrating beams
%    nz         :  integration points for work integral
%    zlim       :  integration limits in units of sphere radius
%    phiout     :  half aperture collection angle of spectrometer

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'lmax', 10 );
addParameter( p, 'nquad', 40 );
addParameter( p, 'nz', 201 );
addParameter( p, 'zlim', 20 );
addParameter( p, 'phiout', 1e-2 );
%  parse input
parse( p, varargin{ : } );

%  initialize Mie solver and save electron velocity
obj.mie = miebeam.miesolver( mat1, mat2, diameter, p.Results.lmax );
obj.vel = vel;
%  save optional input parameters
for name = [ "lmax", "nquad", "nz", "zlim", "phiout" ]
  obj.( name ) = p.Results.( name );
end

