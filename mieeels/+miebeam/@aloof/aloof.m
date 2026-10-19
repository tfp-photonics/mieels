classdef aloof
  %  Mie solution for EELS and aloof electron beams.
  
  properties
    mie       %  Mie solver
    vel       %  electron velocity in units of speed of light
    lmax      %  cutoff for angular degrees
  end

  properties (Hidden)
    M,N       %  Gegenbauer polynomials for aloof beams
  end
  
  methods
    function obj = aloof( varargin )
      %  Initialize Mie solver for EELS and aloof electron beams.
      %
      %  Usage :
      %    obj = miebeam.aloof( mat1, mat2, diameter, vel, PropertyPairs )
      %  Input
      %    mat1       :  material at sphere inside
      %    mat2       :  material at sphere outside
      %    diameter   :  diameter of sphere
      %    vel        :  electron velocity in units of speed of light
      %  PropertyName
      %    lmax       :  cutoff for angular degrees
      obj = init( obj, varargin{ : } );
    end  
  end
end
