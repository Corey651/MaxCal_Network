

%Jij will be 40 neurons and 4 time steps, so a 160x160 matrix
%it will be stacked, neurons 1 2... of first time then second time, etc.


%% initialize ground truth.  This is stochastic, so the produced data will
%not be identical to the attached data.

N=40;
T=4;
J=zeros(N*T);
MJ=.015;
sigJ=.015;
a=4;
s0=randi(2,N*T,1);
s0=s0-1;
s0=2*s0-1;  %either +1 or -1
Equiltime=5000;
Corrtime=500;
Numtime=500000;
for t=1:T
    for s=1:t-1
        for i=1:N
            for j=1:N
                J((t-1)*N+i,(s-1)*N+j)=(randn*sigJ*a^(-abs(t-s))+MJ*a^(-abs(t-s)));
            end
            J((t-1)*N+i,(s-1)*N+i)=20*MJ*a^(-abs(t-s));
        end
    end
end
J=J+J';
for t=1:T
    for i=1:N
        for j=1:i-1
            J((t-1)*N+i,(t-1)*N+j)=randn*sigJ+MJ;
            J((t-1)*N+j,(t-1)*N+i)=J((t-1)*N+i,(t-1)*N+j);
        end
    end
end
h=-.1*ones(N*T,1);



%% simulate ground truth to extract means and correlations

[Havg,FM,SM,~]=MCMC(J,h,N*T,s0,Equiltime,Corrtime,Numtime);  %Using MCMC.m
CC=SM-FM*FM';
FFM=FM;


%% calculate LC estimates, Eq. 7.

hUC=zeros(N*T,1);
JLC=-inv(CC);
for k=1:N*T
    JLC(k,k)=0;
end


%% Solves the Uncoupled problem exactly, Eq. 5

pH=@(L,s1,s2,s3,s4) exp(L(1).*s1+L(2).*s2+L(3).*s3+L(4).*s4+L(5).*s1.*s2+L(6).*s1.*s3+L(7).*s1.*s4+L(8).*s2.*s3+L(9).*s2.*s4+L(10).*s3.*s4);
logZ=@(L) log(pH(L,1,1,1,1)+pH(L,1,1,1,-1)+pH(L,1,1,-1,1)+pH(L,1,1,-1,-1)+pH(L,1,-1,1,1)+pH(L,1,-1,1,-1)+pH(L,1,-1,-1,1)+pH(L,1,-1,-1,-1)+pH(L,-1,1,1,1)+pH(L,-1,1,1,-1)+pH(L,-1,1,-1,1)+pH(L,-1,1,-1,-1)+pH(L,-1,-1,1,1)+pH(L,-1,-1,1,-1)+pH(L,-1,-1,-1,1)+pH(L,-1,-1,-1,-1));
JUC=zeros(N*T,N*T);
for i=1:N
    S=@(L) logZ(L)-L(1)*FM(i)-L(2)*FM(N+i)-L(3)*FM(2*N+i)-L(4)*FM(3*N+i)-L(5)*SM(i,N+i)-L(6)*SM(i,2*N+i)-L(7)*SM(i,3*N+i)-L(8)*SM(N+i,2*N+i)-L(9)*SM(N+i,3*N+i)-L(10)*SM(2*N+i,3*N+i);
    L0=zeros(10,1);
    options=optimset('MaxIter',10000,'MaxFunEvals',10000);
    LL=fminsearch(S,L0,options);
    LL=fminsearch(S,LL,options);

    hUC(i)=LL(1);
    hUC(N+i)=LL(2);
    hUC(2*N+i)=LL(3);
    hUC(3*N+i)=LL(4);
    JLC(N+i,i)=LL(5);
    JUC(N+i,i)=LL(5);
    JLC(2*N+i,i)=LL(6);
    JUC(2*N+i,i)=LL(6);
    JLC(3*N+i,i)=LL(7);
    JUC(3*N+i,i)=LL(7);
    JLC(2*N+i,N+i)=LL(8);
    JUC(2*N+i,N+i)=LL(8);
    JLC(3*N+i,N+i)=LL(9);
    JUC(3*N+i,N+i)=LL(9);
    JLC(3*N+i,2*N+i)=LL(10);
    JUC(3*N+i,2*N+i)=LL(10);
    JLC(i,N+i)=JLC(N+i,i);
    JUC(i,N+i)=JUC(N+i,i);
    JLC(i,2*N+i)=JLC(2*N+i,i);
    JUC(i,2*N+i)=JUC(2*N+i,i);
    JLC(i,3*N+i)=JLC(3*N+i,i);
    JUC(i,3*N+i)=JUC(3*N+i,i);
    JLC(N+i,2*N+i)=JLC(2*N+i,N+i);
    JUC(N+i,2*N+i)=JUC(2*N+i,N+i);
    JLC(N+i,3*N+i)=JLC(3*N+i,N+i);
    JUC(N+i,3*N+i)=JUC(3*N+i,N+i);
    JLC(2*N+i,3*N+i)=JLC(3*N+i,2*N+i);
    JUC(2*N+i,3*N+i)=JUC(3*N+i,2*N+i);
end
%% Calculate LC corrections to h using the mean-field equation. Eq. 6.
hLC=zeros(N*T,1);
for i=1:N
    for t=1:T
        hLC((t-1)*N+i)=hUC((t-1)*N+i);
        for j=1:N
            if j~=i
                for s=1:T
                    hLC((t-1)*N+i)=hLC((t-1)*N+i)-JLC((t-1)*N+i,(s-1)*N+j)*FM((s-1)*N+j);
                end
            end
        end
    end
end



%% Plot estimates vs ground truth, commented to suppress output
Joffdiag=JLC-JUC;
JJ=reshape(J,[N*N*T*T,1]);
JJP=reshape(JLC,[N*N*T*T,1]);
JJO=reshape(Joffdiag,[N*N*T*T,1]);
Iod=find(JJO~=0);
%scatter(JJ(Iod),JJP(Iod))
%scatter(JJ,JJP);
% JJod=JJ(Iod);
% JJPod=JJP(Iod);
% csvwrite('Brain_K_pair.txt',[JJod,JJPod])
% Id=find(JJO==0);
% JJd=JJ(Id);
% JJPd=JJP(Id);
% csvwrite('Brain_K_self.txt',[JJd,JJPd])
% csvwrite('Brain_h_neg_point_1.txt',[hPlef])

%% Calculate dynamics and synchrony for UC/LC vs ground truth vs just UC
Time=20000;
EQtime=500;
s=zeros(N,Time+EQtime+T);
R=randi(2,N,T);
R=R-1;
R=2*R-1;  %either +1 or -1
s(:,1:T)=R;
htemp=zeros(N,1);
Jtemp=J(3*N+1:4*N,3*N+1:4*N);
for tt=T+1:Time+EQtime+T
    for i=1:N
        htemp(i)=h(3*N+i);
        for j=1:N
            for ss=1:T-1
                htemp(i)=htemp(i)+J((ss-1)*N+j,3*N+i)*s(j,tt-T+ss);
            end
        end
    end
    [~,FM]=MCMC(Jtemp,htemp,N,s(:,tt-1),1000,1,1);
    s(:,tt)=FM;
end
M=mean(s);

htemp=zeros(N,1);
Jtemp=JLC(3*N+1:4*N,3*N+1:4*N);
for tt=T+1:Time+EQtime+T
     for i=1:N
        htemp(i)=hLC(3*N+i);
        for j=1:N
            for ss=1:T-1
                htemp(i)=htemp(i)+JLC((ss-1)*N+j,3*N+i)*s(j,tt-T+ss);
            end
        end
    end
    [~,FM]=MCMC(Jtemp,htemp,N,s(:,tt-1),1000,1,1);
    s(:,tt)=FM;
end
MLC=mean(s);

htemp=zeros(N,1);
Jtemp=JUC(3*N+1:4*N,3*N+1:4*N);
for tt=T+1:Time+EQtime+T
     for i=1:N
        htemp(i)=hUC(3*N+i);
        for j=1:N
            for ss=1:T-1
                htemp(i)=htemp(i)+JUC((ss-1)*N+j,3*N+i)*s(j,tt-T+ss);
            end
        end
    end
    [~,FM]=MCMC(Jtemp,htemp,N,s(:,tt-1),1000,1,1);
    s(:,tt)=FM;
end
MUC=mean(s);

% csvwrite('Brain_time_series.txt',[M',MMF',MPlef'])

%% Calculate synchrony variances to see how our estimates change with
%temperature (interaction strength)

%This takes a long time to run.  Be wary.  We attach the output just in
%case

beta=.1:.005:1.25;
varb=zeros(length(beta),3);  %first column is truth, second is MF+LR, third is just MF
magb=zeros(length(beta),3);
for l=300:length(beta)
    Jb=beta(l)*J;
    hb=beta(l)*h;
    
    [Havg,FMb,SMb,~]=MCMC(Jb,hb,N*T,s0,Equiltime,Corrtime,Numtime);
    CCb=SMb-FMb*FMb';
    FFMb=FMb;

    hUCb=zeros(N*T,1);
    JLCb=-inv(CCb);
    for k=1:N*T
        JLCb(k,k)=0;
    end
                
    pH=@(L,s1,s2,s3,s4) exp(L(1).*s1+L(2).*s2+L(3).*s3+L(4).*s4+L(5).*s1.*s2+L(6).*s1.*s3+L(7).*s1.*s4+L(8).*s2.*s3+L(9).*s2.*s4+L(10).*s3.*s4);
    logZ=@(L) log(pH(L,1,1,1,1)+pH(L,1,1,1,-1)+pH(L,1,1,-1,1)+pH(L,1,1,-1,-1)+pH(L,1,-1,1,1)+pH(L,1,-1,1,-1)+pH(L,1,-1,-1,1)+pH(L,1,-1,-1,-1)+pH(L,-1,1,1,1)+pH(L,-1,1,1,-1)+pH(L,-1,1,-1,1)+pH(L,-1,1,-1,-1)+pH(L,-1,-1,1,1)+pH(L,-1,-1,1,-1)+pH(L,-1,-1,-1,1)+pH(L,-1,-1,-1,-1));
    JUCb=zeros(N*T,N*T);
    for i=1:N
        S=@(L) logZ(L)-L(1)*FMb(i)-L(2)*FMb(N+i)-L(3)*FMb(2*N+i)-L(4)*FMb(3*N+i)-L(5)*SMb(i,N+i)-L(6)*SMb(i,2*N+i)-L(7)*SMb(i,3*N+i)-L(8)*SMb(N+i,2*N+i)-L(9)*SMb(N+i,3*N+i)-L(10)*SMb(2*N+i,3*N+i);
        L0=zeros(10,1);
        options=optimset('MaxIter',10000,'MaxFunEvals',10000);
        LL=fminsearch(S,L0,options);
        LL=fminsearch(S,LL,options);
        hUCb(i)=LL(1);
        hUCb(N+i)=LL(2);
        hUCb(2*N+i)=LL(3);
        hUCb(3*N+i)=LL(4);
        JLCb(N+i,i)=LL(5);
        JUCb(N+i,i)=LL(5);
        JLCb(2*N+i,i)=LL(6);
        JUCb(2*N+i,i)=LL(6);
        JLCb(3*N+i,i)=LL(7);
        JUCb(3*N+i,i)=LL(7);
        JLCb(2*N+i,N+i)=LL(8);
        JUCb(2*N+i,N+i)=LL(8);
        JLCb(3*N+i,N+i)=LL(9);
        JUCb(3*N+i,N+i)=LL(9);
        JLCb(3*N+i,2*N+i)=LL(10);
        JUCb(3*N+i,2*N+i)=LL(10);
        JLCb(i,N+i)=JLCb(N+i,i);
        JUCb(i,N+1)=JUCb(N+i,i);
        JLCb(i,2*N+i)=JLCb(2*N+i,i);
        JUCb(i,2*N+i)=JUCb(2*N+i,i);
        JLCb(i,3*N+i)=JLCb(3*N+i,i);
        JUCb(i,3*N+i)=JUCb(3*N+i,i);
        JLCb(N+i,2*N+i)=JLCb(2*N+i,N+i);
        JUCb(N+i,2*N+i)=JUCb(2*N+i,N+i);
        JLCb(N+i,3*N+i)=JLCb(3*N+i,N+i);
        JUCb(N+i,3*N+i)=JUCb(3*N+i,N+i);
        JLCb(2*N+i,3*N+i)=JLCb(3*N+i,2*N+i);
        JUCb(2*N+i,3*N+i)=JUCb(3*N+i,2*N+i);
    end
    hLCb=zeros(N*T,1);
    for i=1:N
        for t=1:T
            hLCb((t-1)*N+i)=hUCb((t-1)*N+i);
            for j=1:N
                if j~=i
                    for s=1:T
                        hLCb((t-1)*N+i)=hLCb((t-1)*N+i)-JLCb((t-1)*N+i,(s-1)*N+j)*FMb((s-1)*N+j);
                    end
                end
            end
        end
    end

    Time=500000;
    EQtime=500;
    s=zeros(N,Time+EQtime+T);
    R=randi(2,N,T);
    R=R-1;
    R=2*R-1;  
    s(:,1:T)=R;
    htemp=zeros(N,1);
    Jtemp=Jb(3*N+1:4*N,3*N+1:4*N);
    for tt=T+1:Time+EQtime+T
        for i=1:N
            htemp(i)=hb(3*N+i);
            for j=1:N
                for ss=1:T-1
                    htemp(i)=htemp(i)+Jb((ss-1)*N+j,3*N+i)*s(j,tt-T+ss);
                end
            end
        end
        [~,FMb]=MCMC(Jtemp,htemp,N,s(:,tt-1),100,1,1);
        s(:,tt)=FMb;
    end
    Mb=mean(s);
    varb(l,1)=var(Mb);
    magb(l,1)=mean(Mb);

    htemp=zeros(N,1);
    Jtemp=JLCb(3*N+1:4*N,3*N+1:4*N);
    for tt=T+1:Time+EQtime+T
        for i=1:N
            htemp(i)=hLCb(3*N+i);
            for j=1:N
                for ss=1:T-1
                    htemp(i)=htemp(i)+JLCb((ss-1)*N+j,3*N+i)*s(j,tt-T+ss);
                end
            end
        end
        [~,FMb]=MCMC(Jtemp,htemp,N,s(:,tt-1),100,1,1);
        s(:,tt)=FMb;
    end
    MLCb=mean(s);
    varb(l,2)=var(MLCb);
    magb(l,2)=mean(MLCb);

    htemp=zeros(N,1);
    Jtemp=JUCb(3*N+1:4*N,3*N+1:4*N);
    for tt=T+1:Time+EQtime+T
        for i=1:N
            htemp(i)=hUCb(3*N+i);
            for j=1:N
                for ss=1:T-1
                    htemp(i)=htemp(i)+JUCb((ss-1)*N+j,3*N+i)*s(j,tt-T+ss);
                end
            end
        end
        [~,FMb]=MCMC(Jtemp,htemp,N,s(:,tt-1),100,1,1);
        s(:,tt)=FMb;
    end
    MUCb=mean(s);
    varb(l,3)=var(MUCb);
    magb(l,3)=mean(MUCb);
end
%%}
% csvwrite('Brain_Temp_mean.txt',[beta',magb])
% csvwrite('Brain_Temp_var.txt',[beta',varb])

