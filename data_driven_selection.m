function [Q_star,detectable,observable] =...
    data_driven_selection(sys,p_star)


C = sys.C;
[p,n] = size(C);

    
ff = zeros(p,1);
fftrue = zeros(p,1);
YhatQ = zeros(sys.Nsamples,1);
for j = 1:p 
    for k = 1:sys.Nsamples
        YhatQ(k) = norm(sys.y(j,sys.tsamples(k)))^2;
    end
    ff(j) = -trace(sys.Eu'*vech_1(sys.PhiTinv*YhatQ)*sys.Eu);
    sys.C = C(j,:);
    fftrue(j) = f(8,sys);
    sys.C = C;
end
%[sort(ff,'descend') sort(fftrue,'descend')]
%error_ff = norm(ff-fftrue) 
error_H2 = abs(sum(ff)-sum(fftrue));
error_H2

[~,Q_star] = sort(ff,'descend');
Q_star = sort(Q_star(1:p_star),'ascend');
%[~,Q_star_true] = sort(fftrue,'descend');
%Q_star_true = sort(Q_star_true(1:p_star),'ascend');
%[Q_star Q_star_true]

% obs(A,C,order,getObs,getDet,getLag,check)
[~,observable,detectable] = obs(sys.A,C(Q_star,:),n,0,1,0,1);

end

