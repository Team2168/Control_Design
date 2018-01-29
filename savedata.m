function savedata( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes heree
x = 0:.1:1;
A = [x; exp(x)];

fileID = fopen('C:/Program Files/MATLAB/exp.txt','w');
fprintf(fileID,'%6s %12s\n','x','exp(x)');
fprintf(fileID,'%6.2f %12.8f\n',A);
fclose(fileID);

end

