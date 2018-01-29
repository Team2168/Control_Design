
%Mech Power (lbf*distance/time * Watt/lbs-ft/s conversion)

close all
n=0:200;
mechPower05in=2000*0.5*0.0833333./n*(1.355817948329);
mechPower1in=2000*1*0.0833333./n*(1.355817948329);
mechPower15in=2000*1.5*0.0833333./n*(1.355817948329);
mechPower2in=2000*2*0.0833333./n*(1.355817948329);


plot(n,mechPower05in,n,mechPower1in,n,mechPower15in,n,mechPower2in);
legend('0.5in Stroke', '1in Stroke', '1.5in Stroke', '2in stroke')
title('Mechanical Power Needed to move 2000lb vs Time for various distances @ 100% Eff')
xlabel('Time for Complete Actuation (Seconds)')
ylabel('Power Required for Actuation (Watts)')
hold on
line([0 200],[21 21])
hold off