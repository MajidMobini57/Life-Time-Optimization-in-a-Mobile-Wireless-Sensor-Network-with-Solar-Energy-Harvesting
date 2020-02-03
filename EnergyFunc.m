function B=EnergyFunc(xn,yn,xs,ys,B_max)     %Figure=0 figure off & Figure=1 figure on

%% Initialization 
efs=.36;
emp=.01;
l=10;
Eelec=.1;

%% Main Loop
N   = numel(xn);                       %Number of Nodes

%% Distance
    d=zeros(1,N);
    for i=1:N
        d(1,i)=norm([xn(i)-xs yn(i)-ys],1);
    end
%% Energy
d0=sqrt(efs/emp);
for i=1:N
    if d(i)<d0
        Ec=(l*Eelec)+(l*efs*d(i)^2);        %Equation for Energy
        B(i,1)=Ec;
    elseif d(i)>=d0
        Ec=(l*Eelec)+(l*emp*d(i)^4);
        B(i,1)=Ec;
    end
end
B=(B./max(B)).*B_max;
end

