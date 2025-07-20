tau=5.4*10000;
plot3(x(1:1000:10000000-2*tau),x(1+tau:1000:10000000-tau),x(1+2*tau:1000:10000000))
xlabel('x')
ylabel('y')
zlabel('z')