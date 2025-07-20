lyap=zeros(1,1000);
j=0;
for(r=0:0.001:1)
    xn1=rand(1);
    %xn1 = 0;
    lyp=0;
    j=j+1;
    for(i=1:10000)
        xn=xn1;
        %logistic map
        %xn1=r*xn*(1-xn);
        xn1=r*sin(pi * xn);
        %wait for transient
        if(i>300)
            % calculate the sum of logaritm
            %lyp=lyp+log(abs(r-2*r*xn1));
            lyp=lyp+log(abs(pi * r * cos(pi*xn)));
        end
    end
    %calculate lyapun
    lyp=lyp/10000;
    lyap(j)=lyp;
end
r=0:0.001:1;
plot(r,lyap);