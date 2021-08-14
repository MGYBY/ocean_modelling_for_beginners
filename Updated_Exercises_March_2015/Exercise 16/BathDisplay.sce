// This script produces a graph of the bathymetry created with "BathCreator.f95".

clf; scf(0); a=gcf(); a.figure_size= [1000,350];
h1=read("topo.dat",-1,153); // read input data
x = (0:150)'; y = (0:50)'; // location vectors  
hzero = max(h1,0.0);
plot3d(x,y,-0.2*hzero(2:52,2:152)',-60,85,' ',flag=[7,2,3]);

