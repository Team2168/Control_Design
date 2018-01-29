function [time,pos,vel,acc,jerk] = doubleSCurve(q0_hat, q1_hat, v0_hat, v1_hat, vmax_hat, amax_hat, jmax_hat)
%% Double S Curve From the Beginning

vmin_hat = -vmax_hat;
amin_hat = -amax_hat;
jmin_hat = -jmax_hat;


%Correct input values for sign
sigma = sign(q1_hat - q0_hat);

q0 = sigma*q0_hat;
q1 = sigma*q1_hat;
v0 = sigma*v0_hat;
v1 = sigma*v1_hat;
t0=0; % Assumes t0 = 0, others a translation in time is needed, see section 5.1

v_max = ((sigma+1)/2)*vmax_hat + ((sigma-1)/2)*vmin_hat;
v_min = ((sigma+1)/2)*vmin_hat + ((sigma-1)/2)*vmax_hat;

a_max = ((sigma+1)/2)*amax_hat + ((sigma-1)/2)*amin_hat;
a_min = ((sigma+1)/2)*amin_hat + ((sigma-1)/2)*amax_hat;

j_max = ((sigma+1)/2)*jmax_hat + ((sigma-1)/2)*jmin_hat;
j_min = ((sigma+1)/2)*jmin_hat + ((sigma-1)/2)*jmax_hat;


%Assuming that vmax and amax are reached, compute the time intervals
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

disp('Take 1')
Ta
Tv
Td
Tj1
Tj2


% Case 2: Where the max velocity is not the velocity limit
if(Tv < 0)

    Tj1 = a_max/j_max;
    Tj2 = a_max/j_max;
    
    delta = (a_max^4/j_max^2) + 2*(v0^2+v1^2) +a_max*(4*(q1-q0) -2*(a_max/j_max)*(v0+v1));
    
    Ta = ((a_max^2/j_max) - 2*v0 + sqrt(delta))/(2*a_max);
    Td = ((a_max^2/j_max) - 2*v1 + sqrt(delta))/(2*a_max);
    Tv = 0;
    disp('Take 2')
Ta
Tv
Td
Tj1
Tj2
    
    %Is there an acceleration time?
    if(Ta < 0 || Td < 0)
        if(Ta < 0)
            Td = 2*((q1-q0)/(v1+v0));
            Tj2 = (j_max*(q1-q0) - sqrt(j_max*(j_max*(q1-q0)^2 + (v1+v0)^2*(v1-v0))))/(j_max*(v1+v0));
            Ta = 0;
        end
        
         if(Td < 0)
            Ta = 2*((q1-q0)/(v1+v0));
            Tj1 = (j_max*(q1-q0) - sqrt(j_max*(j_max*(q1-q0)^2 - (v1+v0)^2*(v1-v0))))/(j_max*(v1+v0));
            Td = 0;
         end
    else
        if(Ta < 2*Tj1 || Td < 2*Tj1)
            % Go back to Step Tv < 0 calc
            % Progressively decrease a_max
            % a_max = lamda*a_max, where 0 < lamda < 1;
        end
    end
end

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

pos = sigma*pos;
vel = sigma*vel;
acc = sigma*acc;
jerk = sigma*jerk;

end