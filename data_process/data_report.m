function [  ] = data_report(  )
close all
clc
cd ../data
load 'ellipse_uniform.mat'

N = size(bounce_array,2);

flag_count = zeros(N,1);
flag_sb_c = flag_count;
flag_en_c = flag_count;
flag_cl_c = flag_count;
flag_mc_c = flag_count;
flag_de_c = flag_count;
flag_ld_c = flag_count;

data_warn = 0;
data_dang = 0;

for i=1:N
    flag_count(i) = sum(bounce_array(i).flags);
    flag_sb_c(i) = bounce_array(i).flag_sb;
    flag_en_c(i) = bounce_array(i).flag_en;
    flag_cl_c(i) = bounce_array(i).flag_cl;
    flag_mc_c(i) = bounce_array(i).flag_mc;
    flag_de_c(i) = bounce_array(i).flag_de;
    flag_ld_c(i) = bounce_array(i).flag_ld;
    
    if flag_count(i) > 1
        data_dang = data_dang+1;
    elseif flag_count(i) == 1
        data_warn = data_warn+1;
    end
end

figure
subplot(2,1,1)
stem(flag_count)
subplot(2,1,2)
hold on
stem(flag_sb_c)
stem(flag_en_c)
stem(flag_cl_c)
stem(flag_mc_c)
stem(flag_de_c)
stem(flag_ld_c)
legend 'Sb' 'En' 'Cl' 'Mc' 'De' 'ld'

prt_str = ['Total number of data %d\n clean %2.2f\n warn: %2.2f \n dang: %2.2f\n'];
fprintf(prt_str,N,(N-(data_warn+data_dang))/N*100,data_warn/N*100,data_dang/N*100)
end

