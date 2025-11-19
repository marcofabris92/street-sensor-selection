function r = my_rank(M)

% r = 0;
% if cond(M) < 1e3 % 1000
%     r = rank(M);
% else
%     s = svd(M);
%     if s(1) == 0
%         return
%     end
%     Etol = 0.99*norm(s)^2;
%     done = 0;
%     Er = 0;
%     while ~done && r < min(size(M))
%         r = r + 1;
%         Er = Er + s(r)^2;
%         if Er >= Etol
%             done = 1;
%         end
%     end
% end

r = rank(M);
    
end