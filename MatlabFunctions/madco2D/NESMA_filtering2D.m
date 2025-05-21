function S_NESMA=NESMA_filtering2D(raw,txy,thresh)

[m,n,p]=size(raw);
S_NESMA=single(zeros(size(raw)));

for i=1:m
    for j=1:n
        if raw(i,j,1)>0.01%10
            rmin=max(i-txy,1);rmax=min(i+txy,m);
            smin=max(j-txy,1);smax=min(j+txy,n);
            L=length(rmin:rmax)*length(smin:smax);
            rawi=reshape(raw(rmin:rmax,smin:smax,:),[L p]);
            x(1,:)=raw(i,j,:);
            D=100*sum(abs(bsxfun(@minus,rawi,x)),2)./sum(x);
            pos=D<thresh;
            S_NESMA(i,j,:)=mean(rawi(pos==1,:),1);
        end
    end
end
end


