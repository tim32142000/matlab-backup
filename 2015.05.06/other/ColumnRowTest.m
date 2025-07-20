clear
tic
A=zeros(10000);
for i=1:10000
    for j=1:10000
        A(j,i)=100000*j+i;
    end
end
toc