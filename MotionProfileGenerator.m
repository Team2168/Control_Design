close all;
clear all;

%inputs:
q0 = 0;
q1 = 10;
v0 = 0;
v1 = 0;
a0 = 0;
a1 = 0;
t0 = 0;
t1 = 3.5;


v_max=4;
a_max=5;
j_max=30;
v_min = -v_max;
a_min = -a_max;
j_min = -j_max;

T=t1-t0;
h=q1-q0;

A0 = q0;
A1 = v0;
A2 = 1/2 * a0;
A3 = 1/(2*T^3) * (20*h - (8*v1 + 12*v0)*T - (3*a0 - a1)*T^2);
A4 = 1/(2*T^4) * (-30*h +(14*v1 + 16*v0)*T + (3*a0 - 2*a1)*T^2);
A5 = 1/(2*T^5) * (12*h - 6*(v1+v0)*T + (a1-a0)*T^2);

P = [A5 A4 A3 A2 A1 A0];
Pd = polyder(P);
Pdd = polyder(Pd);

time = t0:1/100:t1;

figure(1)
plot(time,polyval(P,time))
title('Quintic Spline Position')

figure(2)
plot(time,polyval(Pd,time))
title('Quintic Spline Deriv')

figure(3)
plot(time,polyval(Pdd,time))
title('Quinitic Spline Accel')

tic
%% S-curves
%Case 1 Where the max velocity traversed is equal to v_max
%calculate the time ranges
if (((v_max - v0) * j_max ) < a_max^2) %max accel is not reached

    Tj1 = sqrt((v_max - v0)/j_max);
    Ta = 2*Tj1;
else
    Tj1 = a_max/j_max;
    Ta = Tj1 + ((v_max - v0)/a_max);
end

if (((v_max - v1) * j_max) < a_max^2) %min accel not reached (a_min = -a_max)
    Tj2 = sqrt((v_max -v1)/j_max);
    Td = Tj2 + (v_max - v1) / a_max;
else
    Tj2 = a_max/j_max;
    Td = Tj2 + ((v_max -v1)/a_max);
end

%duration of constant velocity
Tv = (q1-q0)/v_max - (Ta/2)*(1+v0/v_max) - (Td/2)*(1+(v1/v_max));

t1=Ta+Tv+Td;
T=t1-t0;

time = t0:1/100:t1;

Ta
Tv
Td
Tj1
Tj2

%compute actual max/min accel and vel
a_lima = j_max*Tj1
a_limd = -j_max*Tj2
v_lim = v0 + (Ta - Tj1)*a_lima

%Computation of Trajectory for q1 > q0

for i=1:length(time)
    %Accel
    if (time(i) <= Tj1)
        pos(i) = q0 + v0*time(i) + j_max*time(i)^3/6;
        vel(i) = v0 + j_max*time(i)^2/2;
        acc(i) = j_max*time(i);
        jerk(i) = j_max;
    elseif (time(i) > Tj1 && time(i) <= Ta - Tj1) 
        pos(i) = q0 + v0*time(i) + (a_lima/6)*(3*time(i)^2 - 3*Tj1*time(i) + Tj1^2);
        vel(i) = v0 + a_lima*(time(i) - (Tj1/2));
        acc(i) = j_max*Tj1;
        jerk(i) = 0;
    elseif (time(i)>  Ta - Tj1 && time(i) <= Ta)
        pos(i) = q0 + (v_lim + v0)*Ta/2 - v_lim*(Ta-time(i)) - j_min*((Ta-time(i))^3/6);
        vel(i) = v_lim + j_min*((Ta-time(i))^2/2);
        acc(i) = -j_min*(Ta-time(i));
        jerk(i) = j_min;
    elseif (time(i) > Ta && time(i) <= Ta + Tv)
        pos(i) = q0 + (v_lim + v0)*(Ta/2)+v_lim*(time(i)-Ta);
        vel(i) = v_lim;
        acc(i) = 0;
        jerk(i) = 0;
    elseif (time(i) > T-Td && time(i) <= T-Td+Tj2)     
        pos(i) = q1 - (v_lim +v1)*(Td/2) + v_lim*(time(i)-T+Td)-j_max*(((time(i)-T+Td)^3)/6);
        vel(i) = v_lim - j_max*((((time(i)-T+Td)^2))/2);
        acc(i) = -j_max*(time(i)-T+Td);
        jerk(i) = j_min;
    elseif (time(i) > T-Td+Tj2 && time(i) <= T-Tj2)    
        pos(i) = q1 - (v_lim+v1)*(Td/2)+v_lim*(time(i)-T+Td) + (a_limd/6)*(3*(time(i)-T+Td)^2 - 3*Tj2*(time(i)-T+Td) + Tj2^2);
        vel(i) = v_lim+a_limd*(time(i)-T+Td-Tj2/2);
        acc(i) = a_limd;
        jerk(i) = 0;
    elseif (time(i) > T-Tj2 && time(i) <=T)    
        pos(i) = q1-v1*(T-time(i))-j_max*((T-time(i))^3/6);
        vel(i) = v1+j_max*((T-time(i))^2/2);
        acc(i) = -j_max*(T-time(i));
        jerk(i) = j_max;        
    end
end

toc


figure(4)
title('Position')
plot(time,pos)

figure(5)
title('vel')
plot(time,vel)

figure(6)
title('accel')
plot(time,acc)

%%
% clear all
% close all
% Double S Curve From the Beginning


%inputs
q0_hat = 0;
q1_hat = 3;
q2_hat = 5;
v0_hat = 0;
v1_hat = 0;
v2_hat = 0;

%Constraints
vmax_hat = 8;
amax_hat =  10; 
jmax_hat =  30;


[time,pos,vel,acc,jerk] = doubleSCurve(q0_hat,q1_hat,v0_hat,v1_hat,vmax_hat,amax_hat,jmax_hat);
[time2,pos2,vel2,acc2,jerk2] = doubleSCurve(q1_hat,q2_hat,v1_hat,v2_hat,vmax_hat,amax_hat,jmax_hat);
time2 = time2 + time(end);
time = [time time2];
pos = [pos pos2];


figure(100)
title('Position')
plot(time,pos)

figure(200)
title('Position2')
plot(time2,pos2)
figure(201)
title('vel')
plot(time2,vel2)
% 
% figure(102)
% title('accel')
% plot(time,acc)
% 
% 

%% Quintic Splines From RBE500


q0 = 0;
q1 = 3;
q2 = 7;
q3 = 10;
v0 = 0;
v3 = 0;


[t, p, pdot] = quinticSplines(q0,q1,q2,q3,v0,v3);

q0 = 0;
q1 = 4;
q2 = 4.5;
q3 = 5;
v0 = 0;
v3 = 0;

[t2, p2, pdot2] = quinticSplines(q0,q1,q2,q3,v0,v3);
%plot position
% figure(1)
% plot(t,p)
% title('Position Trajectory')
% xlabel('seconds')
% ylabel('meters')
% 
% 
% %plot velocity
% figure(2)
% plot(t,pdot)
% title('Velocity Trajectory')
% xlabel('seconds')
% ylabel('meters/second')
% 
% figure(3)
% plot(p,p2)
% title('global')

%% B-Splines

savedata()

