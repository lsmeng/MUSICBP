function plotAll1(ret)
x=ret.x;

sr=ret.sr;
[n m]=size(x);
t=(1:m)/sr;

hold on;
if isfield(ret,'ploty')
    ploty=ret.ploty;
else
    ploty=1:n;
end
if isfield(ret,'plotcolor');
    plotcolor=ret.plotcolor;
else
    plotcolor='b';
end
for i=1:n
    
plot(t,x(i,1:m)*ret.scale+ploty(i),plotcolor);
xlabel('Time (s)');
ylabel('Station Index');
box on;
end