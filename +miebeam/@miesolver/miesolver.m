classdef miesolver
  %  Mie solver.

  properties
    mat1          %  material at sphere inside
    mat2          %  material at sphere outside
    diameter      %  diameter of sphere
    lmax          %  angular degree cutoff
  end

  methods
    function obj = miesolver( varargin )
      %  Initialize Mie solver.
      %
      %  Usage :
      %    obj = miebeam.miesolver( mat1, mat2, diameter, lmax )
      %  Input
      %    mat1         :  material at sphere inside
      %    mat2         :  material at sphere outside
      %    diameter     :  diameter of sphere
      %    lmax         :  angular degree cutoff
      [ obj.mat1, obj.mat2, obj.diameter, obj.lmax ] = deal( varargin{ : } );
    end
  end
end