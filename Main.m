clear ; close all; clc; format long;
%% SYSTEM PARAMETERS
rng(1)
K = 80;
F = 5;
N = 10;
S_max = 10;
Si_r = 0;
Sd_r = 0.05;
tr_r = 0.46;
r_r =0.25;
B_max = 1;% B_min = 0.001;
rho = 0.98;
w1 = 1; w2 = 0;
%% SET THE UPPER LIMIT OF repeat > 1 FOR REPEATED LIFETIME COMPUTATIONS
for repeat = 1:1
%% INITIAL COMPUTATIONS
% Min_lifetime_no = S_max/(B_max*(1/N))+1;
% Max_lifetime_no = S_max/(B_min*(1/N))+1;
% % % 
% % % M_corr=zeros(K*F,K*F);
% % % for i=1:K*F
% % %     for j=1:K*F
% % %         M_corr(i,j) = 1-(abs(i-j)*(1-rho));
% % %     end
% % % end
% % % M_corr = M_corr.*(M_corr>=0);
% % % for i = 1:K*F
% % %     for j = max(K*F-1,i):K*F
% % % if i ~= j 
% % %     M_corr(i, j) = 2 * sin(pi * M_corr(i, j) / 6);
% % % M_corr(j, i) = 2 * sin(pi * M_corr(j, i) / 6);
% % % end
% % %     end 
% % % end
% % % rng(2)
% % % Ch = chol(M_corr); G = randn(N,K*F); Br = G * Ch;
% % % Br = normcdf(Br)*(B_max-B_min) + B_min;
Min = -1; Max =  1;                    %Mahdode Harkat

% Nodes Data
xn  = unifrnd(Min,Max,[1,N]);   %Nodes Position
yn  = unifrnd(Min,Max,[1,N]);   %Nodes Position

% Sink Data
xs=0;
ys=0;

Br=EnergyFunc(xn,yn,xs,ys,B_max);

for f = 1:F
    B(:,f) = Br(:,1);
end
rng(3)
Si = S_max*ones(N,1).*(1-Si_r*rand(N,1)); 
Si_no = Si;
tr1 = floor(tr_r*K*F);
Activity = []; Energy = [];
% z= w1*(Si)+w2*(Si-B);
%% OPTIMIZATION EVENTS (or runs)
BBr=[Br];
BB=[B];
RR=[];
rr=SolarRecharge(K,F,N);
for k = 1:K
    S=[Si];
    X=[];
    SSi=[];
    R=[];
    for f=1:F
        r=rr(:,k*f);
        cvx_begin quiet
        variables Xo(N,1) So(N,1) z
        minimize(z)
        subject to
        So-z<=0;
        Xo >= 0;
        sum(Xo) == 1;
        So == Si - B(:,f).*Xo  + r.*(1-Xo);          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % z== w1*(Si)+w2*(Si-B);
        cvx_end
        X=[X Xo]; Si=So;SSi=[SSi So];R=[R r.*(1-Xo)];
    end
    % R = zeros(N,K*F);
    % aa=rand(N,1);
    % aaa=aa/sum(aa);
    % R(:,tr1) = r_r*S_max*aaa;
    % R(:,tr1)
    % sum(R(:,tr1))
    % % a=X(:,4);
    % % b=sort(a);
    % % A=zeros(1,length(a));
    % % for i=1:length(a)/2 
    % %     
    % %     a1=find(a==b(i));
    % %     a2=find(a==b(length(a)+1-i));
    % %     A(a1)=b(length(a)+1-i);
    % %     A(a2)=b(i);
    % % end
    % %  R = zeros(N,K*F);
    % %  R(:,tr1) = r_r*S_max.*A;
    % R = zeros(N,K*F);
    % A=ones(10,1)-X(:,2);
    % R(:,tr1)=(r_r*S_max).*(A/sum(A));
    for f = 1:F 
         S(:,f+1) = S(:,f) - B(:,f).*X(:,f) + r.*X(:,f);
    end
    S = S.*(S>=0);
    S=min(S_max,S);
    C=0.05;
    xn=xn+C*randn(1,N);
    yn=yn+C*randn(1,N);

    % Clipping
    xn=max(xn,Min);
    xn=min(xn,Max);
    yn=max(yn,Min);
    yn=min(yn,Max);

    Br=EnergyFunc(xn,yn,xs,ys,B_max);

    for f = 1:F
        B(:,f) = Br(:,1);%(S(:,f)-S(:,f+1)+R(:,f+(k-1)*F))./X(:,f); 
    end
    Activity = [Activity X];
    Energy = [Energy S(:,(1:F))];
    Si = S(:,F+1);
    BB=[BB B];
    BBr=[BBr Br];
    RR=[RR R];
end
save DATA
%% FINAL COMPUTATIONS
Sd = Sd_r*S_max;
TD = ceil(find(Energy <= Sd, 1)/N);
if isempty(TD) == 1; TD = K*F;
    errordlg('Death instant beyond K*F. Re-define system parameters.', 'WARNING');
end; td(repeat) = TD;
x_no = ones(N,1)*(1/N);
Energy_no = [Si_no]; Energy_greedy = [Si_no];
for t = 1:K*F
    Energy_no(:,t+1) = Energy_no(:,t)    +RR(:,t)- BB(:,t) .* x_no ;
    g_greedy         = Energy_greedy(:,t)+RR(:,t) -BB(:,t); 
    g_greedy=min(S_max,g_greedy);

    x_greedy(:,t) = g_greedy >= max(g_greedy); 
    Energy_greedy(:,t+1) = Energy_greedy(:,t) - BB(:,t).* x_greedy(:,t) + RR(:,t);
end
Energy_greedy = Energy_greedy(:,(1:K*F));
Energy_no = Energy_no(:,(1:K*F));
td_no(repeat) = ceil(find(Energy_no<=Sd,1)/N);
TD_greedy = ceil(find(Energy_greedy<=Sd,1)/N);
if isempty(TD_greedy) == 1; TD_greedy = K*F;
errordlg('Death instant beyond K*F. Re-define system parameters.', 'WARNING');
end; td_greedy(repeat) = TD_greedy;
Energy_no = Energy_no.*(Energy_no>=0); 
Energy_greedy = Energy_greedy.*(Energy_greedy>=0); 
%% OUTPUT DATA AND GRAPHS
figure(1); plot(Activity'); hold on; plot(sum(Activity/N),'y','lineWidth',5);
legend('optimized: solid lines','Average: thick line','Location','Best');
axis([1 K*F 0 1.02*max(max(Activity))]); xlabel('Frame index');
ylabel('Activity levels'); title('Proposed optimization method'); hold off;
figure(2); plot(Energy','lineWidth',2); hold on; plot(Energy_no','--');
plot(td,Sd,'ok','lineWidth',5); plot(td_no,Sd,'or','lineWidth',5); grid on;
axis([1 K*F 0 S_max ]); xlabel('Frame index'); ylabel('Residual energies');
legend('optimized: solid lines', 'non-optimized: dashed lines','Location','Best');
title('Proposed optimization method');
hold off;
figure(3); plot(x_greedy'); hold on; plot(x_no','--','lineWidth',2);
axis([1 K*F 0 1.02*max(max(x_greedy))]);
title('Greedy algorithm');
xlabel('Frame index'); ylabel('Activity levels'); hold off;
figure(4); plot(Energy_greedy','lineWidth',2); hold on; plot(Energy_no','--');
plot(td_greedy,Sd,'ok','lineWidth',5); plot(td_no,Sd,'or','lineWidth',5);grid on;
axis([1 K*F 0 S_max ]); xlabel('Frame index'); ylabel('Residual energies');
legend('Greedy: solid lines', 'non-optimized: dashed lines','Location','Best');
title('Greedy algorithm');
hold off;
figure(5)
Xg = x_greedy';
Xo = Activity';
Avg_X_greedy = sum(Xg(1:td_greedy(repeat),:))/(td_greedy(repeat));
Avg_X = sum(Xo(1:td(repeat),:))/(td(repeat));
bar(Avg_X_greedy,0.8); hold on; bar(Avg_X,0.5,'y');
legend('greedy', 'optimized','Location','Best'); xlabel('Sensor node index');
ylabel('Average activity over the lifespan'); hold off;
% figure(6); plot(Br(1:5,:)','Linewidth',3);
figure(6);
for n=1:4
    plot(BBr(n,:),'o--','Linewidth',1);hold on
end
xlim([0,K])
xlabel('Iteration index'); ylabel('Maximum consumption of 4 sensor nodes'); hold off;
disp(['Optimized lifespan: ', num2str(td)]);
disp(['Greedy optimized lifespan: ', num2str(td_greedy)]);
disp(['Run count: ', num2str(repeat)]);
% disp(['System parameters: ','K=',num2str(K),' F=',num2str(F),' N=',num2str(N),' S_max=',num2str(S_max),' Si_r=',num2str(Si_r),' Sd_r=',num2str(Sd_r),' tr_r=',num2str(tr_r),' r_r=',num2str(r_r),' B_min=',num2str(B_min),' B_max=',num2str(B_max),' rho=',num2str(rho),' w1=',num2str(w1),' w2=',num2str(w2)]);
disp(['Optimized mean lifespan: ', num2str(sum(td)/sum(td>0))]);
disp(['Greedy mean lifespan: ', num2str(sum(td_greedy)/sum(td_greedy>0))]);
disp(['Optimized lifespan standard deviation: ', num2str(sqrt(var(td)))]);
disp(['Greedy lifespan standard deviation: ', num2str(sqrt(var(td_greedy)))]);
disp(['Optimized mean lifespan improvement:', num2str(100*(sum(td./td_no)/sum(td>0)-1)),' %']);
disp(['Greedy mean lifespan improvement:', num2str(100*(sum(td_greedy./td_no)/sum(td>0)-1)),' %']);
end