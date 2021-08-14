//******************************************
// This is the Scilab script for Exercise 3.
//
// Author: J. Kaempf, 2015 (updated)
//*******************************************

clf; scf(0); a=gcf(); a.figure_size=[600,300]; 

y=read("output.txt",-1,4); // read input data

for i=2:180 // animation loop

time = int(y(i,1)/60); //time in minutes

z(1:2) = -50.; // equilibrium depth of object
  
if y(i-1,2)>z(1); col = 2; else; col = 5; end; // set different colours

drawlater;
plot2d(y(i-1:i,1)/60,y(i-1:i,2),col);
c = gce(); c.children.thickness = 2;
b = gca(); b.data_bounds=[0,-100;30,0]; b.font_size = 2;
 b.auto_ticks = ["off","off","on"]; b.sub_ticks = [3,3]; 
  b.x_ticks = tlist(["ticks", "locations","labels"],..
 [0 5 10 15 20 25 30], ["0" "5" "10" "15" "20" "25" "30"]);
  b.y_ticks = tlist(["ticks", "locations","labels"],..
 [-100 -80 -60 -40 -20 0], ["-100" "-80" "-60" "-40" "-20" "0"]);  

plot2d(y(i-1:i,1)/60,z(1:2),1);//,'000','',[0,-100,30,0],[2,4,2,5]);
c = gce(); c.children.thickness = 2;
  
xstring(10, -97,"time (minutes)");  // add x label
c = gce(); c.font_size = 3;
xstring(-3, -35,"z (m)");
c = gce(); c.font_size = 3; c.clip_state = "off";

if col == 5; // add title
 title("Too light",'Fontsize',3, 'color','red');
else;
 title("Too heavy",'Fontsize',3,'color','blue');
end;
   
drawnow;

// conversion to GIF files (optional)
//if i < 10 then
 // xs2gif(0,'ex100'+string(i)+'.gif')
//else
 // if i < 100 then
   // xs2gif(0,'ex10'+string(i)+'.gif')
  //else
  //  xs2gif(0,'ex1'+string(i)+'.gif')
 // end
//end

xpause(1d4); // use this to slow down animation

end

