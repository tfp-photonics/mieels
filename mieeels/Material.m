classdef Material < matlab.mixin.Heterogeneous
  %  Material properties of linear materials.
  
  properties
    epsfun    %  function for permittivity
    mufun     %  function for permeability
              %    functions return eps, mu in units of free-space values
  end
  
  properties (Hidden) 
    zsqrt = @sqrt  %  complex square root function for resonance modes
  end
  
  methods
    function obj = Material( varargin )
      %  Set material properties.
      %
      %  Usage :
      %    obj = Material( epsfun, mufun )
      %  Input
      %    epsfun   :  constant value or function handle for permittivity
      %    mufun    :  constant value or function handle for permeability
      %                  val = fun( ene ) internally uses energies in eV
      obj = init( obj, varargin{ : } );
    end

    function x = eps( obj, k0 )
      %  EPS - Permittivity.
      x = obj.epsfun( k0 );
    end

    function x = mu( obj, k0 )
      %  MU - Permeability.
      x = obj.mufun( k0 );
    end

    function x = n( obj, k0 )
      %  N - Refractive index.
      x = sqrt( obj.mufun( k0 ) ) .* sqrt( obj.epsfun( k0 ) );
    end

    function x = Z( obj, k0 )
      %  Z - Impedance.
      x = obj.zsqrt( obj.mufun( k0 ) ) ./ obj.zsqrt( obj.epsfun( k0 ) );
    end

    function x = k( obj, k0 )
      %  K - Wavenumber of light in medium.
      x = k0 .* obj.zsqrt( obj.mufun( k0 ) ) .* obj.zsqrt( obj.epsfun( k0 ) );
    end
  end  %  end methods
  
  
  methods (Access = private)
    function obj = init( obj, eps, mu )
      %  INIT - Initialize material properties.
      
      %  material functions use eV, internal functions use k0 (1/nm)
      eV2nm = 1 / 8.0655477e-4;
      
      %  permittivity
      if isa( eps, 'numeric' )
        obj.epsfun = @( k0 ) eps;
      else
        obj.epsfun = @( k0 ) eps( eV2nm / ( 2 * pi ) * k0 );
      end      
     %  permeability
      if isa( mu, 'numeric' )
        obj.mufun = @( k0 ) mu;
      else
        obj.mufun = @( k0 ) mu( eV2nm / ( 2 * pi ) * k0 ); 
      end
    end       
  end  %  end private methods
  
end  %  end classdef
