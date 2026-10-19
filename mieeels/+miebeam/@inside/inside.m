classdef inside
  %  Mie solution for EELS and inside beams.
  
  properties
    mie       %  Mie solver 
    vel       %  electron velocity in units of speed of light
    lmax      %  cutoff for angular degrees
  end

  properties (Hidden)
    nquad     %  number of quadrature points
    nz        %  integration points for work integral
    zlim      %  integration limits in units of sphere radius
    phiout    %  half aperture collection angle of spectrometer
  end
  
  methods
    function obj = inside( varargin )
      %  Initialize Mie solver for EELS and inside beams.
      %
      %  Usage :
      %    obj = miebeam.inside( mat1, mat2, diameter, vel, PropertyPairs )
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
      obj = init( obj, varargin{ : } );
    end  
  end
end
