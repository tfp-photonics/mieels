classdef solver
  %  Mie solution for EELS and aloof and inside electron beams.
  
  properties
    aloof     %  Mie solution for aloof electron beams
    inside    %  Mie solution for inside electron beams
    diameter  %  diameter of sphere
  end

  
  methods
    function obj = solver( varargin )
      %  Initialize Mie solver for EELS.
      %
      %  Usage :
      %    obj = miebeam.solver( mat1, mat2, diameter, vel, PropertyPairs )
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
      obj.aloof  = miebeam.aloof(  varargin{ : } );
      obj.inside = miebeam.inside( varargin{ : } );
      obj.diameter = obj.aloof.mie.diameter;
    end  
  end
end
