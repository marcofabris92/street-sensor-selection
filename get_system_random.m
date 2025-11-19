function [A, B, C, Ch, n, m, p, r, N, tmax, ngain, pnew, a] = get_system_random()
    n=25; % n
    m=10; % m
    p=n; % p
    pnew=7; % p_star
    r=5; % p_tilde
    A=randn(n, n); % A
    B=randn(n, m); % B
    C=randn(p, n); % C
    
    rho = max(abs(eig(A)));
    A = A/max(1.1*rho,1);
    
   
    Ch=C(1:r, :); % C_tilde
    
    N=ceil(n/r);
    tmax=round(2*(N*m+N*r)*(N*m+N*r+1)/2); % tK
    ngain=0.001; % noisegain

    % a is eta
    if abs(eigs(A, 1, 'lm'))>0.5 
        a=0.5/abs(eigs(A, 1, 'lm'));
    else
        a=1;
    end
end

