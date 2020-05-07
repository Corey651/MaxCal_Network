function [Havg,FM,SM,H] =MCMC(J,h,n,s0,Equiltime,Corrtime,Numtim)

    % J,h are the Lagrange multipliers
    %s0 is our initial starting point, a column of -1's and 1's
    %Equiltime is the equilibration time of the Markov chain
    %Corrtime is the correlation time
    %Numtim is the number of correlation times simulated
   
    M=Numtim*Corrtime+Equiltime;
    H=zeros(1,M+1);
   
    Fm=zeros(n,1);
    Sm=zeros(n,n);  
    
    % make sure that J is symmetric with 0 diagonals
    for i=1:n
        for j=1:i-1
            J(i,j)=J(j,i);
        end
        J(i,i)=0;
    end
    %calculate initial Hamiltonian functions
    G=zeros(n,1);
    for i=1:n
        G(i)=G(i)+h(i);
        for j=1:n
            G(i)=G(i)+J(i,j)*s0(j);
        end
    end
    for i=1:n
        H(1)=H(1)+h(i)*s0(i);
        for j=1:n
            H(1)=H(1)+.5*J(i,j)*s0(j)*s0(i);
        end
    end


    %now do MCMC for M iterations
    s=s0;
    %calculate spin flip choices
    Flip=rand(M,1)*n;
    %single flip basic choice
    Flip=floor(Flip)+1;
    MCMC1=rand(M,1);
    for t=1:Equiltime
        spin=Flip(t);
        negenergy=s(spin)*G(spin);
        change=0;
        if negenergy < 0 
            s(spin)=-s(spin);
            change=1;
            H(t+1)=H(t)-2*negenergy;
        else
            if MCMC1(t)<exp(-2*negenergy)
                s(spin)=-s(spin);
                change=1;
                H(t+1)=H(t)-2*negenergy;
            end
        end
        %now update G's
        if change==1
            for i=1:n
                if i ~= spin
                    G(i)=G(i)+2*J(i,spin)*s(spin);
                end
            end
        else
            H(t+1)=H(t);
        end
    end

    for k=1:Numtim
        zed=Equiltime+(k-1)*Corrtime;
        for t=zed+1:zed+Corrtime
            spin=Flip(t);
            negenergy=s(spin)*G(spin);
            change=0;
            if negenergy < 0 
                s(spin)=-s(spin);
                change=1;
                H(t+1)=H(t)-2*negenergy;
            else
                if MCMC1(t)<exp(-2*negenergy)
                    s(spin)=-s(spin);
                    change=1;
                    H(t+1)=H(t)-2*negenergy;
                end
            end
            %now update G's
            if change==1
                G=G+2*s(spin)*J(:,spin);
            else
                H(t+1)=H(t);
            end
        end
        Fm=Fm+s;
        Sm=Sm+s*s';
    end
    % Calculate means and variances
    Havg=H(Equiltime+Corrtime:Corrtime:end);
    FM=Fm/Numtim;
    SM=Sm/Numtim;
end