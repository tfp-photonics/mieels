function obj = init( obj, mat1, mat2, diameter, vel, varargin )
%  INIT - Initialize Mie solver for EELS and aloof electron beams.
%
%  Usage for obj = miebeam.aloof :
%    obj = init( obj, mat1, mat2, diameter, vel, PropertyPairs )
%  Input
%    mat1       :  material at sphere inside
%    mat2       :  material at sphere outside
%    diameter   :  diameter of sphere
%    vel        :  electron velocity in units of speed of light
%  PropertyName
%    lmax       :  cutoff for angular degrees

%  set up parser
p = inputParser;
p.KeepUnmatched = true;
addParameter( p, 'lmax', 10 );
%  parse input
parse( p, varargin{ : } );

%  initialize Mie solvers
obj.mie = multipole.miesolver( mat1, mat2, diameter, 'lmax', p.Results.lmax );
%  electron velocity
obj.vel = vel;
obj.lmax = p.Results.lmax;

%  initialize Gegenbauer polynomials for nondispersive embedding medium
%  for constant material parameters the evaluation for k0=0 should work
try
  [ obj.M, obj.N ] = beamcoeffs( obj, 0 );
end
