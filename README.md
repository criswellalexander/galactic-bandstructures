# galactic-bandstructures

ExtraGalacticBandStructuresv2.m: Primary code. Creates galactic disk, places dwarf galaxies, calls integrators, outputs and saves figures.
  Note: In order to run the program on your own computer, you will need to change the path called when defining the 'savelocale' variable to a location of your choice.
  Note: Surface density profiles are not produced in the current version. This is intentional.
  
integrosdg.m: First integrator. Directs dwarf galaxy motion based on galactic gravitational potential.
integro3dstars.m: Second integrator. Directs star particle motion based on galactic gravitational potential and influence of dwarf galaxies.
