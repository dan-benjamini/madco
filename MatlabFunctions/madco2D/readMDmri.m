function [ExpPars, ExpData] = readMDmri(im_first_ind,pathname)

for ii=1:length(im_first_ind)
    im_indRef(ii)=im_first_ind(ii);
    im_ind_T2_1D(ii,:)=im_first_ind(ii):(im_first_ind(ii)+19);
    im_ind_D_1D_M3(ii)=im_ind_T2_1D(end)+1;im_ind_D_1D_M4(ii)=im_ind_D_1D_M3(ii)+1;im_ind_D_1D_M6(ii)=im_ind_D_1D_M4(ii)+1;
    im_ind_T1_1D(ii,:)=(im_ind_D_1D_M6(ii)+1):(im_ind_D_1D_M6(ii)+20);
    im_ind_T1ref_1D(ii)=im_indRef(ii);
    im_ind_DT1_M3(ii,:)=[im_ind_T1_1D(ii,end)+1,im_ind_T1_1D(ii,end)+4,im_ind_T1_1D(ii,end)+7,im_ind_T1_1D(ii,end)+10];
    im_ind_DT1_M4(ii,:)=im_ind_DT1_M3(ii,:)+1; im_ind_DT1_M6(ii,:)=im_ind_DT1_M4(ii,:)+1;
    im_ind_DT2_M3(ii,:)=[im_ind_DT1_M6(ii,end)+1,im_ind_DT1_M6(ii,end)+4,im_ind_DT1_M6(ii,end)+7,im_ind_DT1_M6(ii,end)+10];
    im_ind_DT2_M4(ii,:)=im_ind_DT2_M3(ii,:)+1;im_ind_DT2_M6(ii,:)=im_ind_DT2_M4(ii,:)+1;
    im_ind_T1T2(ii,:)=[im_ind_DT2_M6(ii,end)+1:im_ind_DT2_M6(ii,end)+16];
    im_indRefHi(ii)=im_ind_DT2_M6(ii,end)+1;
end


for ii=1:length(im_first_ind)


    [tmp,~]=read_2dseq_3D_TEpath(im_indRefHi(ii),pathname{ii});
    ImageRefHi(ii,:,:,:)=squeeze(tmp);

    % reference image
    [tmp,~]=read_2dseq_3D_TEpath(im_indRef(ii),pathname{ii});
    ImageRef(ii,:,:,:)=squeeze(tmp);

    % 1D-T1 data
    [tmp,TI(ii,:)]=read_2dseq_3D_TIpath(im_ind_T1_1D(ii,:),pathname{ii});
    ImageT1(ii,1:end-1,:,:,:)=tmp;
    [tmp,~]=read_2dseq_3D_TEpath(im_ind_T1ref_1D(ii),pathname{ii});
    ImageT1(ii,end,:,:,:)=tmp;

    % 1D-D data
    [ImageRefM3,ImageM3, bvalM3,TE_D_M3(ii) ]=read_2dseq_3D_DiffIsoAcqPath(im_ind_D_1D_M3(ii),pathname{ii});
    [ImageRefM4,ImageM4, bvalM4,TE_D_M4(ii) ]=read_2dseq_3D_DiffIsoAcqPath(im_ind_D_1D_M4(ii),pathname{ii});
    [ImageRefM6,ImageM6, bvalM6,TE_D_M6(ii) ]=read_2dseq_3D_DiffIsoAcqPath(im_ind_D_1D_M6(ii),pathname{ii});


    ImageRefD(ii,:,:,:)=(ImageRefM3+ImageRefM4+ImageRefM6)/3;

    nM3=length(bvalM3);nM4=length(bvalM4);nM6=length(bvalM6);
    n3=nM3-nM4;
    n4=nM4-nM6;
    n6=nM6;

    M3=geo_mean(ImageM3,5);
    M4=geo_mean(ImageM4,5);
    M6=geo_mean(ImageM6,5);

    ImageIso3=M3(:,:,:,1:n3);
    ImageIso4=1/5*(2*M3(:,:,:,(n3+1):(n3+n4))+3*M4(:,:,:,1:n4));
    ImageIso6=1/7*(10/5*M3(:,:,:,n3+n4+1:end) + 9/5*M4(:,:,:,n4+1:end) + 16/5*M6);


    ImageIsoDOrg=cat(4,ImageIso3,ImageIso4,ImageIso6);
    b_1Dorg=1e6*bvalM3;

    Dstart=1;

    ImageIsoD(ii,:,:,:,:)=ImageIsoDOrg(:,:,:,Dstart:end);
    b_1D(ii,:)=b_1Dorg(Dstart:end);

    % 1D-T2 data
    [ImageT2(ii,:,:,:,:),TE(ii,:)]=read_2dseq_3D_TEpath(im_ind_T2_1D(ii,:),pathname{ii});


    % D-T2

    for kk=1:length(im_ind_DT2_M3(ii,:))


        [ImageRefM3TE,ImageM3TE, bvalM3TE,TE_DT2_M3(ii,kk) ]=read_2dseq_3D_DiffIsoAcqPath(im_ind_DT2_M3(ii,kk),pathname{ii});
        [ImageRefM4TE,ImageM4TE, bvalM4TE,TE_DT2_M4(ii,kk) ]=read_2dseq_3D_DiffIsoAcqPath(im_ind_DT2_M4(ii,kk),pathname{ii});
        [ImageRefM6TE,ImageM6TE, bvalM6TE,TE_DT2_M6(ii,kk) ]=read_2dseq_3D_DiffIsoAcqPath(im_ind_DT2_M6(ii,kk),pathname{ii});


        nM3TE=length(bvalM3TE);nM4TE=length(bvalM4TE);nM6TE=length(bvalM6TE);
        n3TE=nM3TE-nM4TE;
        n4TE=nM4TE-nM6TE;
        n6TE=nM6TE;

        M3TE=geo_mean(ImageM3TE,5);
        M4TE=geo_mean(ImageM4TE,5);
        M6TE=geo_mean(ImageM6TE,5);

        ImageIso3TE=M3TE(:,:,:,1:n3TE);
        ImageIso4TE=1/5*(2*M3TE(:,:,:,(n3TE+1):(n3TE+n4TE))+3*M4TE(:,:,:,1:n4TE));
        ImageIso6TE=1/7*(10/5*M3TE(:,:,:,n3TE+n4TE+1:end) + 9/5*M4TE(:,:,:,n4TE+1:end) + 16/5*M6TE);

        ImageIsoDT2(ii,:,:,:,:,kk)=cat(4,ImageIso3TE,ImageIso4TE,ImageIso6TE);

        b_DT2(ii,:)=1e6*bvalM3TE;
    end



    % D-T1
    %         clear ImageIsoDT1
    for kk=1:length(im_ind_DT1_M3(ii,:))


        [ImageRefM3IR,ImageM3IR, bvalM3IR,TI_DT1_M3(ii,kk) ]=read_2dseq_3D_DiffIsoAcqIRpath(im_ind_DT1_M3(ii,kk),pathname{ii});
        [ImageRefM4IR,ImageM4IR, bvalM4IR,TI_DT1_M4(ii,kk) ]=read_2dseq_3D_DiffIsoAcqIRpath(im_ind_DT1_M4(ii,kk),pathname{ii});
        [ImageRefM6IR,ImageM6IR, bvalM6IR,TI_DT1_M6(ii,kk) ]=read_2dseq_3D_DiffIsoAcqIRpath(im_ind_DT1_M6(ii,kk),pathname{ii});



        nM3IR=length(bvalM3IR);nM4IR=length(bvalM4IR);nM6IR=length(bvalM6IR);
        n3IR=nM3IR-nM4IR;
        n4IR=nM4IR-nM6IR;
        n6IR=nM6IR;

        M3IR=geo_mean(ImageM3IR,5);
        M4IR=geo_mean(ImageM4IR,5);
        M6IR=geo_mean(ImageM6IR,5);

        ImageIso3IR=M3IR(:,:,:,1:n3IR);
        ImageIso4IR=1/5*(2*M3IR(:,:,:,(n3IR+1):(n3IR+n4IR))+3*M4IR(:,:,:,1:n4IR));
        ImageIso6IR=1/7*(10/5*M3IR(:,:,:,n3IR+n4IR+1:end) + 9/5*M4IR(:,:,:,n4IR+1:end) + 16/5*M6IR);

        ImageIsoDT1(ii,:,:,:,:,kk)=cat(4,ImageIso3IR,ImageIso4IR,ImageIso6IR);

        b_DT1(ii,:)=1e6*bvalM3IR;
    end


    % T1-T2
    [ImageT1T2(ii,:,:,:,:),TE_T1T2v(ii,:),TI_T1T2v(ii,:)]=read_2dseq_3D_TETIpath(im_ind_T1T2(ii,:),pathname{ii});
    %         TI_T1T2(ii,:)=unique(TI_T1T2v);TE_T1T2(ii,:)=unique(TE_T1T2v);%TE_T1T2(1)=[];


end

ExpPars.TI=TI;
ExpPars.TE=TE;
ExpPars.b_1D=b_1D;
ExpPars.b_DT2=b_DT2;
ExpPars.TE_DT2_M3=TE_DT2_M3;
ExpPars.b_DT1=b_DT1;
ExpPars.TI_DT1_M3=TI_DT1_M3;
ExpPars.TE_T1T2v=TE_T1T2v;
ExpPars.TI_T1T2v=TI_T1T2v;

ExpData.ImageRef=ImageRef;
ExpData.ImageT2=ImageT2;
ExpData.ImageT1=ImageT1;
ExpData.ImageRefD=ImageRefD;
ExpData.ImageIsoD=ImageIsoD;
ExpData.ImageIsoDT2=ImageIsoDT2;
ExpData.ImageIsoDT1=ImageIsoDT1;
ExpData.ImageT1T2=ImageT1T2;


