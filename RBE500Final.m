%% RBE500 Harrilal Final Report

function RBE500Final
question1
question2
question3b
question3d
question4
question5
end

function question1

A = [ 1 0 0 0 0 0;
    1 1 1 1 1 1;
    1 2 4 8 16 32;
    1 3 9 27 81 243;
    0 1 0 0 0 0
    0 1 6 27 108 405];

B = [ 0; 0.2; 0.8; 1; 0; 0.1];

%solve for coeffs of polynomial
y = A\B


%reverse y to put in form for polyval
r=wrev(y)
t = 0:.1:3;

p = polyval(r,t);

%plot position
figure(1)
plot(t,p)
title('Position Trajectory')
xlabel('seconds')
ylabel('meters')

k = polyder(r);
pdot = polyval(k,t);
%plot velocity
figure(2)
plot(t,pdot)
title('Velocity Trajectory')
xlabel('seconds')
ylabel('meters/second')
%<sprintf('\L')>
end

function question2

H0_1 = Tx(5)*Rx(-60)*Tz(4)*Ry(-90)*Rz(90)*Tx(-4*sqrt(3))

H1_0 = inv(H0_1)

end

%

function question3b

figure(3)
[T,Y] = ode45(@equation3,[0 12],[0 0]);
plot(T,Y(:,1),'-',T,Y(:,2),'-.'); grid on;
legend('Y', 'dY');
title('Question 3B ODE Solution')

end

function question3d
num = [15];
den = [1 6 17 54 72];
sys = tf(num,den)

figure(4)
step(sys)
end

function question4

syms theta1 d2 theta3 theta4

H0_1 = Ty(3)*Rzrad(theta1)
H1_2 = Ty(5)*Tx(d2)
H2_3 = Rx(-90)*Tz(8)*Rzrad(theta3)
H3_4 = Rx(90)*Tx(4)*Rzrad(theta4)
H4_5 = Tx(2)

H0_5 = H0_1*H1_2*H2_3*H3_4*H4_5

end



function question5
num = [1];
den = [1 4 6];
sys = tf(num,den)


sys_cl = feedback(sys,1)


%PD controller
Kp = 4;
Kd = 1;
C = pid(Kp,0,Kd)
T = feedback(C*sys,1)

figure(5)
t = 0:0.01:4;
u = exp(-2*t).*sin(4*t);
lsim(sys,u,t)   % u,t define the input signal



end




%supporting functions
function dy = equation3(t, y)
dy(1,:) = y(2);
dy(2,:) = 5*sin(3*t) -6*y(2) - 8*y(1);
end

%Function which performs the 2 dimensional rotation
%Input argument is theta in Degrees
function A=R2D(theta)

A=[cosd(theta),-sind(theta); 
    sind(theta), cosd(theta)];

end

%Function which performs the 3 dimensional rotation
%about X-Axis
%Input argument is theta in Degrees
function A=Rx(theta)

A=[1,0,0,0; 
    0, cosd(theta),-sind(theta), 0; 
    0 sind(theta), cosd(theta), 0;
    0,0,0,1];

end

%Function which performs the 3 dimensional rotation
%about y-Axis
%Input argument is theta in Degrees
function A=Ry(theta)

A=[cosd(theta), 0 , sind(theta),0; 
    0, 1, 0, 0; 
    -sind(theta), 0, cosd(theta), 0;
    0,0,0,1];

end

%Function which performs the 3 dimensional rotation
%about Z-Axis
%Input argument is theta in Degrees
function A=Rz(theta)

A=[cosd(theta), -sind(theta), 0, 0; 
    sind(theta) cosd(theta), 0, 0; 
    0,0,1, 0;
    0,0,0,1];

end

%Function which performs the 3 dimensional rotation
%about Z-Axis
%Input argument is theta in radians or symbolic
function A=Rzrad(theta)

A=[cos(theta), -sin(theta), 0, 0; 
    sin(theta) cos(theta), 0, 0; 
    0,0,1, 0;
    0,0,0,1];

end

%Function which performs the 3 dimensional Translation
%Translation is along X-axis
%Input argument is distance in meters
function A=Tx(d)

A=[1, 0, 0, d;
   0, 1, 0, 0;
   0, 0, 1, 0;
   0, 0, 0, 1];

end

%Function which performs the 3 dimensional Translation
%Translation is along y-axis
%Input argument is distance in meters
function A=Ty(d)

A=[1, 0, 0, 0;
   0, 1, 0, d;
   0, 0, 1, 0;
   0, 0, 0, 1];

end

%Function which performs the 3 dimensional Translation
%Translation is along Z-axis
%Input argument is distance in meters
function A=Tz(d)

A=[1, 0, 0, 0;
   0, 1, 0, 0;
   0, 0, 1, d;
   0, 0, 0, 1];

end

