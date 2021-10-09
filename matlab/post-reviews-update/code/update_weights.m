function [w]=update_weights(logw,loglike)
LogW1=logw+loglike;

maxLogW1=max(LogW1);   
diffW=LogW1-maxLogW1;
%update weights
w1=exp(diffW);
% normalize
sw=sum(w1);
if sw==0
    disp(['filter degeneracy'])
    w1=w1+realmin;
    sw=sum(w1);
end
w=w1./sw;
return
end

%%% see https://doi.org/10.1155/2018/5763461
%% to try: Log-PF: Particle Filtering in Logarithm Domain
