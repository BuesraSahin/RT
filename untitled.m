clear 
close all

radius = 20;

front_lens = sqrt(radius^2 - (x_front-radius).^2); 
figure(1)
plot(40-x_optaxis,raymatrix','r') %Rays
hold on

%Plano-convex lens
line([40-thickness 40-thickness], [max(front_lens) -max(front_lens)],'color','black')
plot(40-x_front,front_lens,'black',40-x_front,-front_lens,'black') %Lens front surface 
plot(x_optaxis,zeros(1,length(x_optaxis)),'k--') %Optical axis
hold off