total=1;

ts=100;

x=ones(ts,1);
y=ones(ts,1);

rx=total/ts;
ry=-total/ts;


for i=2:timestep
    x(i)=x(i-1)*(1+rx);
    y(i)=y(i-1)*(1+ry);
end

value=(x+y)/2;