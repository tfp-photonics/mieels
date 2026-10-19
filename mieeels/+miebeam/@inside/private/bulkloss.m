function P = bulkloss( obj, k0 )
%  BULKLOSS - Bulk losses for EELS.
%
%  Usage for obj = miebeam.inside :
%    P = bulkloss( obj, k0 )
%  Input
%    k0   :  wavenumber of light in vacuum
%  Output
%    P    :  bulk losses per unit length

%  electron rest mass in eV and reduced Compton wavelength
mass = 0.51e6;
compton = 386e-6;
%  ratio of cutoff wavenumber of spectrometer and Q
eV2nm = 1 / 8.0655477e-4;
ene = eV2nm * k0 / ( 2 * pi );
cutoff = sqrt( 1 + ( mass * obj.vel * obj.phiout / ene ) ^ 2 );
  
%  dielectric function
mie = obj.mie;
epsz = mie.mat1.eps( k0 );
%  bulk loss, Garcia de Abajo, RMP 2010, Eq. (18)
P = 1 / ( 137 * pi * compton * mass * obj.vel ^ 2 ) *  ...
  imag( ( obj.vel ^ 2 - 1 / epsz ) .*  ...
   log( ( cutoff ^ 2 - obj.vel ^ 2 * epsz ) ./ ( 1 - obj.vel ^ 2 * epsz ) ) );
