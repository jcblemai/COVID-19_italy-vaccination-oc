function [indn]=systematic_resampling(w,NSample)
u=rand*1/NSample;
csum=0;
j=0;
for cont_sample=1:NSample
    while csum<u
        j=j+1;
        csum=csum+w(j);
        
    end
    indn(cont_sample)=j;
    u=u+1/NSample;
end
return
end
