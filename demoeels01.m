%  DEMOEELS01 - EELS and radiative loss probabilities for fixed wavelength.

%  nanosphere
diameter = 160;
a = 0.5 * diameter;

%  dielectric functions
mat1 = Material( 1, 1 );
% mat2 = Material( epstable( 'gold.dat' ), 1 );
mat2 = Material( 4, 1 );

k0 = 2 * pi / 600;

%  imapct parameters and electron velocity
bimp1 = linspace( 0, a - 2, 14 );
bimp2 = linspace( a + 5, 3 * a, 20 );
bimp = [ bimp1, bimp2 ];
vel = 0.4;

%  Mie solver
lmax = 10;
mie = miebeam.solver( mat2, mat1, diameter, vel, 'lmax', lmax );
%  EELS excitation
P = loss( mie, bimp, k0 );

%%  final figure
figure

ind = bimp < a;

opts = { 'LineWidth', 1, 'MarkerSize', 6 };

h1 = plot( bimp( ~ind ), P.surf( ~ind ), '+-', opts{ : } );  hold on
h2 = plot( bimp( ~ind ), P.rad( ~ind ),  'o-', opts{ : } );  

set( gca, 'ColorOrderIndex', 1 );
plot( bimp( ind ), P.surf( ind ), '+-', opts{ : } );  
plot( bimp( ind ), P.rad( ind ),  'o-', opts{ : } );  
h3 = plot( bimp( ind ), P.bulk( ind ), 'k-', 'LineWidth', 1 );

xl = xline( a, '-', 'Radius' );
xl.LabelVerticalAlignment = 'bottom';
xl.FontSize = 12;

legend( [ h1, h2, h3 ], 'surf', 'rad', 'bulk' )

xlabel( 'Impact parameter (nm)' );
ylabel( 'Loss probability (1/eV)' );
