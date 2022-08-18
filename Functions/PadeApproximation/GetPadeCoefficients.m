close all

A=[]; B=[]; C=[]; D=1;
sys = ss(A,B,C,D,'InputDelay',0.03);
sysx=pade(sys,1);
sysx2=pade(sys,2);
sysx3=pade(sys,3);
sysx4=pade(sys,4);
sysx_0=pade(sys,0);
%% Step Response
t=0:0.001:0.5;
u=ones(size(t));
out_perfect=lsim(sys,u,t);
out_dyn=lsim(sysx,u,t);
out_dyn2=lsim(sysx2,u,t);
out_dyn3=lsim(sysx3,u,t);
out_dyn4=lsim(sysx4,u,t);
figure();
subplot(1,3,1)
plot(t,u,'k');hold on;
plot(t,out_perfect,'r');
plot(t,out_dyn,'b');
plot(t,out_dyn2,'g');
plot(t,out_dyn3,'m');
plot(t,out_dyn4,'y');
legend('Input','Perfect','Order 1','Order 2','Order 3','Order 4');
set(gca,'XLim',[t(1)-0.1 t(end)+0.1]);
set(gca,'YLim',[-1.1 u(end)+0.1]);

%% Response to Sine wave
t=0:0.001:5;
u=10*sin(2*t);
out_perfect=lsim(sys,u,t);
out_dyn=lsim(sysx,u,t);
out_dyn2=lsim(sysx2,u,t);
out_dyn3=lsim(sysx3,u,t);
out_dyn4=lsim(sysx4,u,t);
subplot(1,3,2)
plot(t,u,'k');hold on;
plot(t,out_perfect,'r');
plot(t,out_dyn,'b');
plot(t,out_dyn2,'g');
plot(t,out_dyn3,'m');
plot(t,out_dyn4,'y');

%% Response to Sine wave
t=0:0.001:1;
u=10*sin(20*t);
out_perfect=lsim(sys,u,t);
out_dyn=lsim(sysx,u,t);
out_dyn2=lsim(sysx2,u,t);
out_dyn3=lsim(sysx3,u,t);
out_dyn4=lsim(sysx4,u,t);
subplot(1,3,3)
plot(t,u,'k');hold on;
plot(t,out_perfect,'r');
plot(t,out_dyn,'b');
plot(t,out_dyn2,'g');
plot(t,out_dyn3,'m');
plot(t,out_dyn4,'y');


%% Response to experimental data

Datapath=fullfile(pwd,'pp9_young','GPOPS_Opt_Input');
matFiles=dir(fullfile(Datapath,'*.mat'));
StateNames={'qa_filt','qh_filt','qa_dot','qh_dot'};


for i=1%:length(matFiles)
    
    dat=load(fullfile(Datapath,matFiles(i).name));
    figure();
    for j=1:length(StateNames)
        subplot(3,2,j)
        u=dat.(StateNames{j});
        t=(0:length(dat.(StateNames{j}))-1)./100;
        out_perfect=lsim(sys,u,t);
        out_dyn=lsim(sysx,u,t);
        out_dyn2=lsim(sysx2,u,t);
        out_dyn3=lsim(sysx3,u,t);
        out_dyn4=lsim(sysx4,u,t);
        plot(t,u,'k');hold on;
        plot(t,out_perfect,'r');
        plot(t,out_dyn,'b');
        plot(t,out_dyn2,'g');
        plot(t,out_dyn3,'m');
        plot(t,out_dyn4,'y');
        title(StateNames{j});
        xlabel('Time [s]');
        ylabel(StateNames{j})
    end
    
    subplot(3,2,5)
    x=[dat.qa_filt(1:end-1) dat.qh_filt(1:end-1) dat.qa_dot(1:end) dat.qh_dot(1:end)];
    K=[-35 -10 6 -2];
    u=x*K';
    t=(0:length(dat.(StateNames{j}))-1)./100;
    out_perfect=lsim(sys,u,t);
    out_dyn=lsim(sysx,u,t);
    out_dyn2=lsim(sysx2,u,t);
    out_dyn3=lsim(sysx3,u,t);
    out_dyn4=lsim(sysx4,u,t);
    plot(t,u,'k');hold on;
    plot(t,out_perfect,'r');
    plot(t,out_dyn,'b');
    plot(t,out_dyn2,'g');
    plot(t,out_dyn3,'m');
    plot(t,out_dyn4,'y');
    title('Full State');
    xlabel('Time [s]');ylabel('x');
    
    
    subplot(3,2,6)
    x=[dat.qa_filt(1:end-1) dat.qh_filt(1:end-1) dat.qa_dot(1:end) dat.qh_dot(1:end)];
    K=[-5 -10 10 -0.1];
    u=x*K';
    t=(0:length(dat.(StateNames{j}))-1)./100;
    out_perfect=lsim(sys,u,t);
    out_dyn=lsim(sysx,u,t);
    out_dyn2=lsim(sysx2,u,t);
    out_dyn3=lsim(sysx3,u,t);
    out_dyn4=lsim(sysx4,u,t);
    plot(t,u,'k');hold on;
    plot(t,out_perfect,'r');
    plot(t,out_dyn,'b');
    plot(t,out_dyn2,'g');
    plot(t,out_dyn3,'m');
    plot(t,out_dyn4,'y');
    title('Full State');
    xlabel('Time [s]');ylabel('x');
end


%% Response to experimental data from event(1)

Datapath=fullfile(pwd,'pp9_young','GPOPS_Opt_Input');
matFiles=dir(fullfile(Datapath,'*.mat'));
StateNames={'qa_filt','qh_filt','qa_dot','qh_dot'};

dt=2;
for i=1%:length(matFiles)
    
    dat=load(fullfile(Datapath,matFiles(i).name));
    figure();
    for j=1:length(StateNames)
        subplot(3,2,j)
        u=dat.(StateNames{j});
        t=(0:length(dat.(StateNames{j}))-1)./100;
        ind0=find(t>=dat.event(1),1,'first');
        indf=find(t>=dat.event(2)+dt,1,'first');
        u=u(ind0:indf);
        t=t(ind0:indf);
        t=t-t(1);
        out_perfect=lsim(sys,u,t);
        out_dyn=lsim(sysx,u,t);
        out_dyn2=lsim(sysx2,u,t);
        out_dyn3=lsim(sysx3,u,t);
        out_dyn4=lsim(sysx4,u,t);
        plot(t,u,'k');hold on;
        plot(t,out_perfect,'r');
        plot(t,out_dyn,'b');
        plot(t,out_dyn2,'g');
        plot(t,out_dyn3,'m');
        plot(t,out_dyn4,'y');
        title(StateNames{j});
        xlabel('Time [s]');
        ylabel(StateNames{j})
    end
    
    subplot(3,2,5)
    x=[dat.qa_filt(1:end-1) dat.qh_filt(1:end-1) dat.qa_dot(1:end) dat.qh_dot(1:end)];
    K=[-35 -10 6 -2];
    u=x*K';
    t=(1:length(dat.qa_filt))./100;
    ind0=find(t>=dat.event(1),1,'first');
    indf=find(t>=dat.event(2)+dt,1,'first');
    u=u(ind0:indf);
    t=t(ind0:indf);
    t=t-t(1);
    out_perfect=lsim(sys,u,t);
    out_dyn=lsim(sysx,u,t);
    out_dyn2=lsim(sysx2,u,t);
    out_dyn3=lsim(sysx3,u,t);
    out_dyn4=lsim(sysx4,u,t);
    plot(t,u,'k');hold on;
    plot(t,out_perfect,'r');
    plot(t,out_dyn,'b');
    plot(t,out_dyn2,'g');
    plot(t,out_dyn3,'m');
    plot(t,out_dyn4,'y');
    title('Full State');
    xlabel('Time [s]');ylabel('x');
    
    subplot(3,2,6)
    x=[dat.qa_filt(1:end-1) dat.qh_filt(1:end-1) dat.qa_dot(1:end) dat.qh_dot(1:end)];
    K=[-5 -10 10 -0.1];
    u=x*K';
    t=(1:length(dat.qa_filt))./100;
    ind0=find(t>=dat.event(1),1,'first');
    indf=find(t>=dat.event(2)+dt,1,'first');
    u=u(ind0:indf);
    t=t(ind0:indf);
    t=t-t(1);
    out_perfect=lsim(sys,u,t);
    out_dyn=lsim(sysx,u,t);
    out_dyn2=lsim(sysx2,u,t);
    out_dyn3=lsim(sysx3,u,t);
    out_dyn4=lsim(sysx4,u,t);
    plot(t,u,'k');hold on;
    plot(t,out_perfect,'r');
    plot(t,out_dyn,'b');
    plot(t,out_dyn2,'g');
    plot(t,out_dyn3,'m');
    plot(t,out_dyn4,'y');
    title('Full State');
    xlabel('Time [s]');ylabel('x');
end





%% Predict forward with time delayed system by integrating ODE
t=0:0.001:5;
e=cos(4*t);
[ta,y]=ode15s(@approx_delay,[0 5],2.*e(1)./5); e_sim=cos(4*ta);
figure();plot(t,u,'b');hold on;plot(ta,y,'r');plot(ta,5*y-e_sim,'g');
legend('e','y','u')

% legend('Input','Perfect','Order 1','Order 2','Order 3','Order 4');
% [ta2,qa2]=ode15s(@approx_delay2,[0 5],u(1));

% s = tf('s');
% sys = exp(-0.1*s);
% sysx = pade(sys,3)
%% Test Pade with longer Time delay


A=[]; B=[]; C=[]; D=1;
sys = ss(A,B,C,D,'InputDelay',0.3);
sysx=pade(sys,1);
