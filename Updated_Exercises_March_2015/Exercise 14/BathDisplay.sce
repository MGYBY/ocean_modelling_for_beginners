// This script produces a graph of the bathymetry created with "BathCreator.f95".

clf; scf(0); a=gcf(); a.figure_size= [500,350];
h1=read("topo.dat",-1,103); // read input data
x = (0:0.1:10)'; y = (0:0.1:5)'; // location vectors  
hzero = max(h1,0.0);
plot3d(x,y,-hzero(2:52,2:102)',-70,5,' ',flag=[7,2,3]);

