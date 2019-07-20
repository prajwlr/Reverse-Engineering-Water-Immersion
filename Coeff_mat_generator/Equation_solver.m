clear variables
close all
clc

n = input('input value of n\n')
nsq = n*n
coeff_matrix = zeros(nsq)

fid = fopen('coeff.txt','r')
if fid<0
    error('Some error\n\n')
else
    Tarr = textscan(fid , " %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ")
end

for i=1:nsq
    for j = 1:nsq
        coeff_matrix(i,j)=Tarr{j}(i)
    end
end

rank(coeff_matrix)

V = zeros(nsq,1)

vf = {4000,5000,6000,7000,2650.42,4678.20,5693.63,6709.05,1550.07,4253.05,5464.59,6528.77,7000,1133.97,3304.22,5309.40}

for i=1:16
    V(i,1)=(vf{i}(1))/1000;
end

ansarr = inv(coeff_matrix)*V

for i=1:16
    if(abs(ansarr(i,1))<500)fprintf("0\n")
    else fprintf("1\n")
    end
end
