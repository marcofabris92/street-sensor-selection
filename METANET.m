function [x,y] = METANET(p,x0,u)

T = length(u)-1;
x = zeros(2*p.n,T+1);
y = zeros(p.n,T+1);
x(:,1) = x0;
for t = 1:T
    rho = x(1:p.n,t);
    v = x(p.n+1:end,t);
    y(:,t) = rho.*v;
    x(1:p.n,t+1) = rho + p.Ts*p.DL*(p.PU*y(:,t)+p.B*u(t));
    % set p.eta_tilde = 0 to restore classic METANET
    rho_hat = rho + ...
        p.eta_tilde*p.kappa_tilde*(p.PD*rho)./(rho+p.kappa_tilde);
    x(p.n+1:end,t+1) = v + ...
        p.Ts/p.tau*(p.vfree.*(exp(-1/p.a*((rho_hat/p.rho_crit).^p.a)))-v);

    % % Enable to restore classic METANET
    % v_conv = p.conv*p.Ts*p.DL*(v.*(p.PU*v-v));
    % v_anti = p.anti*p.Ts*p.DL/p.tau*(p.PD*rho-rho)./(rho+p.kappa);
    % w = 1.1;
    % for i = p.n+1:2*p.n
    %     D_v_i = x(i,t+1) + v_conv(i-p.n) - v_anti(i-p.n);
    %     if D_v_i >= 0 && D_v_i <= w*p.vfree
    %         x(i,t+1) = D_v_i;
    %     elseif D_v_i > w*p.vfree
    %         x(i,t+1) = w*p.vfree;
    %     else
    %         x(i,t+1) = 0;
    %     end
    % end
end
y(:,T+1) = x(1:p.n,T+1).*x(p.n+1:end,T+1);

x = x';
y = y';

end

