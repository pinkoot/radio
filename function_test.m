function PSP=function_test()    %declaration of function;
%function_test is the function which generates maximum-length sequence;
v=[-1 -1 -1 -1 -1 -1 -1 -1 -1]; % declaration of register's vector 
                                %(n=9) for GNSS, initial conditions for 
                                %MATLAB is -1;
for i= 1:511                  %the beginning cycle for all numbers ;       
    PSP(i)=v(7);              %read off the output from the seven register;
    f=v(5)*v(9);                %multiplication 
    v(2:9)=v(1:8);              %
    v(1)=f;
end
PSP=-PSP;

 