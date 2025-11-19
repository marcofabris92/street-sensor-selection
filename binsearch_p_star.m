function p_star = binsearch_p_star(sys,metric,want_obs)

%% Linear search
% p = size(sys.C,1);
% n = size(sys.C,2);
% p_star = 1;
% done = 0;
% while ~done && p_star <= n
%     fprintf(['Trying with ' num2str(p_star) ' sensor(s).\n'])
%     [Q_star,detectable] = exaustive_selection(sys,@f,metric,p_star);
% 
%     % Get the selection matrix
%     S = zeros(p_star,p);
%     for k = 1:p_star
%         S(k,Q_star(k)) = 1;
%     end
% 
%     % Check whether p_star is ok
%     if my_rank(S*sys.C) == p_star && detectable
%         done = 1;
%     else
%         p_star = p_star + 1;
%     end
% end

%% Binary search
p_star = 1;
lower = 1;
upper = min(size(sys.C,1),size(sys.C,2));

while upper > p_star
    % Computes sensor selection
    [~,detectable,observable] = exaustive_selection(sys,@f,metric,p_star);

    % Updating lower and upper bounds for p_star
    condition = want_obs*observable + (1-want_obs)*detectable;
    if condition
        upper = p_star;
        p_star = floor((lower+upper)/2);
    else
        lower = p_star;
        p_star = ceil((lower+upper)/2);
    end
end


end