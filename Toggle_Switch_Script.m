%Toggle switch MaxCal MF

    %script to simulate the Toggle Switch

hP=-4.605;
hS=7.6;
KC=-exp(1-hP-hS);
KK=[-.35,KC,-.05];

for i=1:3
    
    K=KK(i);           % Running the simulation for Each K
    T=10000000;        % The length of time
    NA=zeros(1,T+1);
    NB=zeros(1,T+1);
    NA(1)=20;
    NB(1)=0;
    U=rand(1,T);
    KBA=K;
    KAB=K;
    hbeta=hP;
    halp=hP;
    hB=hS;
    hA=hS;
    
    % Calculating the probabilities of each event
    
    p1=exp(KBA*NB(1)+halp);                     %Production of A
    p2=exp(KAB*NA(1)+hbeta);                    %Production of B
    p3=NA(1)*exp(-hA);                          %Degradation of A
    p4=NB(1)*exp(-hB);                          %Degradation of B
    p5=1;                                       %Nothing
    p6=NB(1)*exp(halp+KBA*(NB(1)-1)-hB);        %Production of A and degradation of B
    p7=NA(1)*exp(hbeta+KAB*(NA(1)-1)-hA);       %Production of B and degradation of A

    P=[p1,p2,p3,p4,p5,p6,p7];
    Z=sum(P);
    II=zeros(1,T);
    la=zeros(1,T);
    lb=zeros(1,T);
    lA=zeros(1,T);
    lB=zeros(1,T);

    for t=1:T
        Ui=U(t)*Z;
        S=double(cumsum(P));
        I=find(S>Ui,1);
        II(t)=I;
        
        % Updating the number of proteins and probabilities
        if I==1
            NA(t+1)=NA(t)+1;
            NB(t+1)=NB(t);
            p2s=p2*exp(KAB);
            p3s=NA(t+1)*exp(-hA);
            p7s=NA(t+1)*exp(hbeta+KAB*(NA(t+1)-1)-hA);
            Z=Z+(p3s-p3)+(p2s-p2)+(p7s-p7);
            p3=p3s;
            p2=p2s;
            p7=p7s;
            P(2)=p2;
            P(3)=p3;
            P(7)=p7;
            la(t)=1;
            lA(t)=NA(t);
            lB(t)=NB(t);
        elseif I==2
            NA(t+1)=NA(t);
            NB(t+1)=NB(t)+1;
            p1s=p1*exp(KBA);
            p4s=NB(t+1)*exp(-hB);
            p6s=NB(t+1)*exp(halp+KBA*(NB(t+1)-1)-hB);
            Z=Z+(p4s-p4)+(p1s-p1)+(p6s-p6);
            p4=p4s;
            p1=p1s;
            p6=p6s;
            P(1)=p1;
            P(4)=p4;
            P(6)=p6;
            lb(t)=1;
            lA(t)=NA(t);
            lB(t)=NB(t);
        elseif I==3
            NA(t+1)=NA(t)-1;
            NB(t+1)=NB(t);
            p2s=p2*exp(-KAB);
            p3s=NA(t+1)*exp(-hA);
            p7s=NA(t+1)*exp(hbeta+KAB*(NA(t+1)-1)-hA);
            Z=Z+(p3s-p3)+(p2s-p2)+(p7s-p7);
            p3=p3s;
            p2=p2s;
            p7=p7s;
            P(2)=p2;
            P(3)=p3;
            P(7)=p7;
            lA(t)=NA(t)-1;
            lB(t)=NB(t);
        elseif I==4
            NA(t+1)=NA(t);
            NB(t+1)=NB(t)-1;
            p1s=p1*exp(-KBA);
            p4s=NB(t+1)*exp(-hB);
            p6s=NB(t+1)*exp(halp+KBA*(NB(t+1)-1)-hB);
            Z=Z+(p4s-p4)+(p1s-p1)+(p6s-p6);
            p4=p4s;
            p1=p1s;
            p6=p6s;
            P(1)=p1;
            P(4)=p4;
            P(6)=p6;
            lA(t)=NA(t);
            lB(t)=NB(t)-1;
        elseif I==5
            NA(t+1)=NA(t);
            NB(t+1)=NB(t);
            lA(t)=NA(t);
            lB(t)=NB(t);
        elseif I==6
            NA(t+1)=NA(t)+1;
            NB(t+1)=NB(t)-1;
            p2s=p2*exp(KAB);
            p3s=NA(t+1)*exp(-hA);
            p7s=NA(t+1)*exp(hbeta+KAB*(NA(t+1)-1)-hA);
            p1s=p1*exp(-KBA);
            p4s=NB(t+1)*exp(-hB);
            p6s=NB(t+1)*exp(halp+KBA*(NB(t+1)-1)-hB);
            Z=Z+(p3s-p3)+(p2s-p2)+(p7s-p7)+(p1s-p1)+(p4s-p4)+(p6s-p6);
            p3=p3s;
            p2=p2s;
            p7=p7s;
            p4=p4s;
            p1=p1s;
            p6=p6s;
            P(2)=p2;
            P(3)=p3;
            P(7)=p7;
            P(4)=p4;
            P(1)=p1;
            P(6)=p6;
            la(t)=1;
            lA(t)=NA(t);
            lB(t)=NB(t)-1;
        else
            NA(t+1)=NA(t)-1;
            NB(t+1)=NB(t)+1;
            p2s=p2*exp(-KAB);
            p3s=NA(t+1)*exp(-hA);
            p7s=NA(t+1)*exp(hbeta+KAB*(NA(t+1)-1)-hA);
            p1s=p1*exp(KBA);
            p4s=NB(t+1)*exp(-hB);
            p6s=NB(t+1)*exp(halp+KBA*(NB(t+1)-1)-hB);
            Z=Z+(p3s-p3)+(p2s-p2)+(p7s-p7)+(p1s-p1)+(p4s-p4)+(p6s-p6);
            p3=p3s;
            p2=p2s;
            p7=p7s;
            p4=p4s;
            p1=p1s;
            p6=p6s;
            P(2)=p2;
            P(3)=p3;
            P(7)=p7;
            P(4)=p4;
            P(1)=p1;
            P(6)=p6;
            lb(t)=1;
            lA(t)=NA(t)-1;
            lB(t)=NB(t);
        end
    end
    csvwrite(strcat('Toggle_K_',num2str(K),'.txt'),[NA',NB'])
end


%find phase diagram from Uncoupled Max Cal, Eq. C10

Lambda=exp(hP+hS);
KK=-.4:.001:0;
NH=zeros(length(KK),1);
NL=zeros(length(KK),1);
N0=zeros(length(KK),1);
for i=1:length(KK)
    f0=@(N0) N0-Lambda*exp(KK(i)*N0);
    Ns=fsolve(f0,Lambda);
    N0(i)=Ns;
    if KK(i)<KC
        f0=@(N) [N(1)-Lambda*exp(KK(i)*N(2)),N(2)-Lambda*exp(KK(i)*N(1))];
        Ns=fsolve(f0,[Lambda,0]);
        NH(i)=Ns(1);
        NL(i)=Ns(2);
    else
        NH(i)=N0(i);
        NL(i)=N0(i);
    end
end
        
csvwrite('Toggle_Phase.txt',[KK',N0,NH,NL])        



