function [ExpDataF,ExpDataFvect] = filterMDmri(ExpData,threshIm)

ImageT2=ExpData.ImageT2;
ImageT1=ExpData.ImageT1;
ImageRefD=ExpData.ImageRefD;
ImageIsoD=ExpData.ImageIsoD;
ImageIsoDT2=ExpData.ImageIsoDT2;
ImageIsoDT1=ExpData.ImageIsoDT1;
ImageT1T2=ExpData.ImageT1T2;

clear rawF raw
for kkAll=1:length(threshIm)
    
    ImageT2org=squeeze(ImageT2(kkAll,:,:,:,:));
    ImageT1org=squeeze(ImageT1(kkAll,:,:,:,:));
    ImageRefDOrg=squeeze(ImageRefD(kkAll,:,:,:,:));
    ImageIsoDOrg=squeeze(ImageIsoD(kkAll,:,:,:,:));
    ImageIsoDT2org=squeeze(ImageIsoDT2(kkAll,:,:,:,:,:));
    ImageIsoDT1org=squeeze(ImageIsoDT1(kkAll,:,:,:,:,:));
    ImageT1T2org=squeeze(ImageT1T2(kkAll,:,:,:,:));

       
    ImageT2org1=permute(ImageT2org,[2 3 4 1]);
    ImageT1org1=permute(ImageT1org,[2 3 4 1]);
    
    clear ImageIsoDT2org1
    kk=1;
    for ii=1:size(ImageIsoDT2org,4)
        for jj=1:size(ImageIsoDT2org,5)
            ImageIsoDT2org1(:,:,:,kk)=ImageIsoDT2org(:,:,:,ii,jj);
            kk=kk+1;
        end
    end
    
    clear ImageIsoDT1org1
    kk=1;
    for ii=1:size(ImageIsoDT1org,4)
        for jj=1:size(ImageIsoDT1org,5)
            ImageIsoDT1org1(:,:,:,kk)=ImageIsoDT1org(:,:,:,ii,jj);
            kk=kk+1;
        end
    end
    
    
    ImageT1T2org1=permute(ImageT1T2org,[2 3 4 1]);
    
    % filter together
    
    raw{kkAll}=cat(4,ImageT2org1,ImageT1org1,ImageRefDOrg,ImageIsoDOrg,ImageIsoDT2org1,ImageIsoDT1org1,ImageT1T2org1);
    
    txy=5; tz=5; thresh=7;
    rawF{kkAll}=NESMA_filtering(raw{kkAll},txy,tz,thresh);
    
end

%% rearranging

indTmp=21;

clear ImageIsoDAll ImageRefHiAll ImageRefAll ImageT1All ImageRefDAll ImageT2All ImageIsoDT2All ImageT1T2All ImageIsoDT1All ImageRefDTIAll ImageDTIAll
for ii=1:length(threshIm)   

    ImageRefAll(ii,:,:,:) = squeeze(rawF{ii}(:,:,:,1));
    
    ImageRefDAll(ii,:,:,:) = squeeze(rawF{ii}(:,:,:,21+indTmp));
    
    tmp = permute(rawF{ii}(:,:,:,1:20),[4 1 2 3]);
    for oo1=1:size(tmp,1)
        ImageT2All(ii,oo1,:,:,:)=squeeze(tmp(oo1,:,:,:));
    end
    tmp = permute(rawF{ii}(:,:,:,21 : (21+indTmp-1)  ),[4 1 2 3]);
    for oo1=1:size(tmp,1)
        ImageT1All(ii,oo1,:,:,:)=squeeze(tmp(oo1,:,:,:));
    end
    tmp = rawF{ii}(:,:,:, (21+indTmp+1):(21+indTmp+16));
    for oo1=1:size(tmp,4)
        ImageIsoDAll(ii,:,:,:,oo1)=squeeze(tmp(:,:,:,oo1));
    end
    
    kk=1;
    for ii1=1:4
        for jj1=1:4
            ImageIsoDT2All(ii,:,:,:,ii1,jj1)=squeeze(rawF{ii}(:,:,:, 21+indTmp+16+kk ));
            kk=kk+1;
        end
    end
    
    kk=1;
    for ii1=1:4
        for jj1=1:4
            ImageIsoDT1All(ii,:,:,:,ii1,jj1)=squeeze(rawF{ii}(:,:,:, 21+indTmp+32+kk ));
            kk=kk+1;
        end
    end
    tmp = permute(rawF{ii}(:,:,:,(21+indTmp+49):(21+indTmp+64)),[4 1 2 3]) ;
    for oo1=1:size(tmp,1)
        ImageT1T2All(ii,oo1,:,:,:)=squeeze(tmp(oo1,:,:,:));
    end
    
    
end

ExpDataF.ImageRef=ImageRefAll;
ExpDataF.ImageT2=ImageT2All;
ExpDataF.ImageT1=ImageT1All;
ExpDataF.ImageRefD=ImageRefDAll;
ExpDataF.ImageIsoD=ImageIsoDAll;
ExpDataF.ImageIsoDT2=ImageIsoDT2All;
ExpDataF.ImageIsoDT1=ImageIsoDT1All;
ExpDataF.ImageT1T2=ImageT1T2All;




%% thresholding for mask
 
clear row_pixAll col_pixAll ImageRefAllMask NcurAll bwMaskCurAll
for ii=1:length(threshIm)
    for iiS=1:size(ImageRefAll,4)
        
        ImTemp=squeeze(ImageRefAll(ii,:,:,iiS));
        ImTemp(ImTemp<threshIm(ii))=0;
        ImTemp(ImTemp~=0)=1;
        bwMaskCurAll(ii,:,:,iiS)=ImTemp;
           
        [row_pixAll{ii,iiS},col_pixAll{ii,iiS}]=find(squeeze(bwMaskCurAll(ii,:,:,iiS)));
     
        NcurAll(ii,iiS)=length(row_pixAll{ii,iiS});
     
       
    end
end



%% vectorize voxels


clear ImageRefAllV ImageT1AllV ImageIsoDAllV ImageRefDAllV ImageT2AllV ImageIsoDT2AllV ImageIsoDT1AllV ImageT1T2AllV
for kkAll=1:length(threshIm)
    for iiS=1:size(squeeze(ImageRefAll(kkAll,:,:,:)),3)
        
        kk=1;
        
        for ii1=1:length(row_pixAll{kkAll,iiS})
            ImageRefAllV{kkAll,iiS}(kk)=squeeze(double(squeeze(ImageRefAll(kkAll,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS))));
            ImageT1AllV{kkAll,iiS}(:,kk)=squeeze(double(squeeze(ImageT1All(kkAll,:,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS))));
            ImageIsoDAllV{kkAll,iiS}(:,kk)=squeeze(double(squeeze(ImageIsoDAll(kkAll,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS,:))));
            ImageRefDAllV{kkAll,iiS}(:,kk)=squeeze(double(squeeze(ImageRefDAll(kkAll,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS))));
            ImageT2AllV{kkAll,iiS}(:,kk)=squeeze(double(squeeze(ImageT2All(kkAll,:,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS))));
            ImageIsoDT2AllV{kkAll,iiS}(:,:,kk)=squeeze(double(squeeze(ImageIsoDT2All(kkAll,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS,:,:))));
            ImageIsoDT1AllV{kkAll,iiS}(:,:,kk)=squeeze(double(squeeze(ImageIsoDT1All(kkAll,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS,:,:))));
            ImageT1T2AllV{kkAll,iiS}(:,kk)=squeeze(double(squeeze(ImageT1T2All(kkAll,:,row_pixAll{kkAll,iiS}(ii1),col_pixAll{kkAll,iiS}(ii1),iiS))));
            
            kk=kk+1;
        end
    end
end



ExpDataFvect.ImageRefAllV=ImageRefAllV;
ExpDataFvect.ImageT1AllV=ImageT1AllV;
ExpDataFvect.ImageIsoDAllV=ImageIsoDAllV;
ExpDataFvect.ImageRefDAllV=ImageRefDAllV;
ExpDataFvect.ImageT2AllV=ImageT2AllV;
ExpDataFvect.ImageIsoDT2AllV=ImageIsoDT2AllV;
ExpDataFvect.ImageIsoDT1AllV=ImageIsoDT1AllV;
ExpDataFvect.ImageT1T2AllV=ImageT1T2AllV;
ExpDataFvect.bwMaskCurAll=bwMaskCurAll;
ExpDataFvect.row_pixAll=row_pixAll;
ExpDataFvect.col_pixAll=col_pixAll;
ExpDataFvect.NcurAll=NcurAll;









