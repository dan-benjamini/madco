function S_NESMA=NESMA_filtering(raw,txy,tz,thresh)

[m,n,o,p]=size(raw);
S_NESMA=single(zeros(size(raw)));

for k=1:o
    disp(['NESMA filtering ... Slice # ' num2str(k) ' of '  num2str(o)]);
    for i=1:m
        for j=1:n
            if raw(i,j,k,1)>10                
                rmin=max(i-txy,1);rmax=min(i+txy,m);
                smin=max(j-txy,1);smax=min(j+txy,n);
                tmin=max(k-tz,1);tmax=min(k+tz,o);
                L=length(rmin:rmax)*length(smin:smax)*length(tmin:tmax);                
                rawi=reshape(raw(rmin:rmax,smin:smax,tmin:tmax,:),[L p]);
                x(1,:)=raw(i,j,k,:);
                D=100*sum(abs(bsxfun(@minus,rawi,x)),2)./sum(x);
                pos=D<thresh;
                S_NESMA(i,j,k,:)=mean(rawi(pos==1,:),1);
            end
        end
    end
end
end


