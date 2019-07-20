clear variables
close all
clc

n = input('input value of n\n')
nsq = n*n
coeff_matrix = zeros(nsq)

fid = fopen('coeff.txt','rt')
if fid<0
    error('Some error\n\n')
end

oneline = fgets(fid)

test_row = zeros(1,nsq)
tstr=''
ptr=1

test_row = str2num(oneline)

coeff_matrix(ptr,:)=test_row

rank(coeff_matrix)