function MultiDprocFun(sampleNum,sliceNum,pathToSave,ExpDataFvect,ExpPars)

%Project Title: multidimensional MRI data processing %%%%%%%
%Author: Dan Benjamini
%Contact: dan.benjamini@nih.gov
% Last Modified by Dan Benjamini 13-Oct-2023 %%%%%%%%%%%%%%%

ImageRefAllV=ExpDataFvect.ImageRefAllV;
ImageT1AllV=ExpDataFvect.ImageT1AllV;
ImageIsoDAllV=ExpDataFvect.ImageIsoDAllV;
ImageRefDAllV=ExpDataFvect.ImageRefDAllV;
ImageT2AllV=ExpDataFvect.ImageT2AllV;
ImageIsoDT2AllV=ExpDataFvect.ImageIsoDT2AllV;
ImageIsoDT1AllV=ExpDataFvect.ImageIsoDT1AllV;
ImageT1T2AllV=ExpDataFvect.ImageT1T2AllV;
bwMaskCurAll=ExpDataFvect.bwMaskCurAll;
row_pixAll=ExpDataFvect.row_pixAll;
col_pixAll=ExpDataFvect.col_pixAll;
NcurAll=ExpDataFvect.NcurAll;
titles=ExpDataFvect.titles;

TI=ExpPars.TI;
TE=ExpPars.TE;
b_1D=ExpPars.b_1D;
b_DT2=ExpPars.b_DT2;
TE_DT2_M3=ExpPars.TE_DT2_M3;
b_DT1=ExpPars.b_DT1;
TI_DT1_M3=ExpPars.TI_DT1_M3;
TE_T1T2=ExpPars.TE_T1T2v;
TI_T1T2=ExpPars.TI_T1T2v;

if isempty(gcp('nocreate'))
    parpool(round(0.96*feature('numcores')));
end



kT2=@(TE,T2)exp(-TE./T2);
kT1=@(TIR,T1)exp(-TIR./T1);
kD=@(b,D)exp(-b.*D);


curImage=1:length(titles);


for kkAll=sampleNum %1:length(curImage)

    sliceNumAdj=sliceNum{kkAll};
    mkdir(fullfile(pathToSave,titles{kkAll}));
    save(fullfile(pathToSave,titles{kkAll},strcat('ws_MultiDproc_Sample',titles{kkAll},'_slice',num2str(sliceNumAdj(1)),'to',num2str(sliceNumAdj(end)),'.mat')));
    mSave = matfile(fullfile(pathToSave,titles{kkAll},strcat('ws_MultiDproc_Sample',titles{kkAll},'_slice',num2str(sliceNumAdj(1)),'to',num2str(sliceNumAdj(end)),'.mat')),'Writable',true);


    for iiS=sliceNumAdj %1:size(ImageRefAllV,2)
        
        
        
        
        ImT1cur=ImageT1AllV{kkAll,iiS};
        ImageIsoDcur=ImageIsoDAllV{kkAll,iiS};        
        ImRefCurD=ImageRefDAllV{kkAll,iiS};
        ImT2cur=ImageT2AllV{kkAll,iiS};
        ImageIsoDT2=ImageIsoDT2AllV{kkAll,iiS};
        ImageIsoDT1=ImageIsoDT1AllV{kkAll,iiS};
        ImageT1T2=ImageT1T2AllV{kkAll,iiS};
        
        
        Ncur=NcurAll(kkAll,iiS);
        if Ncur==0
            continue
        end
        
       
     
        % flip, subtract and normalize with Ref
        
        clear TIR TR min_ind
        
        ImT1curMod=NaN(size(ImT1cur,1)-2,size(ImT1cur,2));
        % ImT1curMod=NaN(size(ImT1cur,1)-1,size(ImT1cur,2),size(ImT1cur,3),length(curSlice));
        
        % TI(end)=[];
        for j=1:Ncur
            
            %             [a1_x,a2_x]=min(squeeze(ImT1cur(:,j)));
            [a1_x,a2_x]=min(squeeze(ImT1cur(1:end-1,j)));
            
            min_ind(j)=a2_x;
            IR_out=min_ind(j);
            TIR(j,:)=TI(~ismember(1:end,IR_out));
            
            %     tmp1=squeeze(ImT1cur(1:(end-1),row_pix(j),col_pix(j)));
            tmp1=squeeze(ImT1cur(1:end-1,j));
            tmp1(IR_out)=[];
            tmp2=tmp1;
            tmp2(1:IR_out-1,:)=-tmp2(1:IR_out-1,:);
            %
            %         tmp3=repmat(ImT1cur(end,row_pix(j),col_pix(j)),size(tmp2,1),1) - tmp2;
            tmp3=repmat(squeeze(ImT1cur(end,j)),size(tmp2,1),1) - tmp2;
            %         tmp3=repmat(squeeze(ImRefCur(row_pix{ii}(j),col_pix{ii}(j),ii)),size(tmp2,1),1) - tmp2;
            
            ImT1curMod(:,j)=tmp3/squeeze(ImT1cur(end,j))/2;
            %     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
            
        end
        
        % T1 1D
        NT1_1D=50;
        
        T1_1DA=logspace(log10(0.001),log10(10),NT1_1D);
        
        TIR_1D=TIR;
        
        clear E_T1
        for j=1:Ncur
            %         E_T1(:,j,ii)=squeeze(ImT1curMod(:,row_pix{ii}(j),col_pix{ii}(j),ii));
            E_T1(:,j)=squeeze(ImT1curMod(:,j));
            
        end
        
        alpha_T1= logspace( -1, 4, 50 );
        % alpha_T1= logspace( 1, 4, 30 );
        
        nT1_1D=NT1_1D;
        nTIR=size(TIR_1D,2);
        
        % 1D-T1
        clear E_mat_T1all
        parfor kk=1:Ncur
            for i=1:nTIR
                for j=1:NT1_1D
                    E_mat_T1all(kk,i,j)=kT1(TIR_1D(kk,i),T1_1DA(j));
                end
            end
        end
        
        
        
        
        %
        clear E_T1til sVecT1 sT1  E_mat_T1til
        for jj=1:Ncur
            [UT1, ST1, VT1]=svd(squeeze(E_mat_T1all(jj,:,:)));
            sT1(jj)=find(diag(ST1)<0.01,1)-1;
            % end
            % sT1=round(mean(sT1_1));
            
            UtilT1=UT1(:,1:sT1(jj));
            VtilT1=VT1(:,1:sT1(jj));
            StilT1=ST1(1:sT1(jj),1:sT1(jj));
            E_mat_T1til{jj}=StilT1*VtilT1';
            E_T1til{jj}(:)=UtilT1'*E_T1(:,jj);
            sVecT1{jj}(:)=diag(StilT1);
        end
        
        
        % 1D-T1
        
        tic
        clear res_T1_1Dsvd norm_T1_1Dsvd x_T1_1D_matsvd res_T1_1Dtmp norm_T1_1Dtmp x_T1_1D_mattmp  ...
            alphaT1Opt GCVfinalT1 alphaT1iterPre alphaT1iter
        clear res_T1_1DRefsvdIter norm_T1_1DRefsvdIter x_T1_1D_matRefsvdIter
        parfor kk=1:Ncur
            E_mat_T1tilCur=squeeze(E_mat_T1til{kk});
            sVecT1Cur=sVecT1{kk}';
            sT1cur=sT1(kk);
            E_sig_sub=E_T1til{kk}';
            GCViterPre=1000;
            alphaT1iter=sum(sVecT1Cur.^2)/sT1cur;
            alphaT1iterPre=alphaT1iter;
            
            criter=1;
            kkcount=1;
            while (criter>=0 && kkcount<50)
                % parfor kk=1:length(alpha_D)
                [res_T1_1DRefsvdIter, norm_T1_1DRefsvdIter, x_T1_1D_matRefsvdIter,~]= opt_cvx_1D_alt2new(E_mat_T1tilCur,E_sig_sub,alphaT1iter,nT1_1D);
                
                cVect=(E_sig_sub-E_mat_T1tilCur*x_T1_1D_matRefsvdIter')./alphaT1iter;
                
                mfree=sum((sVecT1Cur.^2)./(sVecT1Cur.^2 + alphaT1iter));
                mfreePrime=-sum(sVecT1Cur.^2./((sVecT1Cur.^2 + alphaT1iter).^2));
                
                
                tHat=sum(((sVecT1Cur.*cVect).^2)./(sVecT1Cur.^2 + alphaT1iter));
                
                % GCViter=sD.*(res_D_1DRefsvdIter.^2)./((sD-mfree).^2);
                
                GCViterUpdate=sT1cur.*alphaT1iter^2.*norm(cVect,2)^2./((sT1cur-mfree).^2);
                
                alphaT1iterUpdate=-(alphaT1iter.^2.*norm(cVect,2).^2 * mfreePrime)./(tHat.*(sT1cur-mfree));
                
                if kkcount==1
                    criter=1;
                else
                    criter=(GCViterUpdate-GCViterPre)/(alphaT1iter-alphaT1iterPre);
                end
                alphaT1iterPre=alphaT1iter;
                alphaT1iter=alphaT1iterUpdate;
                GCViterPre=GCViterUpdate;
                % GCViter=GCViterPre;
                % clear alphaDiterUpdate GCViter
                
                kkcount=kkcount+1;
                
            end
            alphaT1Opt(kk)=alphaT1iterPre;
            GCVfinalT1(kk)=GCViterPre;
            res_T1_1Dsvd(kk)=res_T1_1DRefsvdIter;
            norm_T1_1Dsvd(kk)=norm_T1_1DRefsvdIter;
            x_T1_1D_matsvd(:,kk)=x_T1_1D_matRefsvdIter;
        end
        
        toc
        
        
        for kk=1:Ncur
            Psi_T1{kkAll,iiS}(kk,:)=x_T1_1D_matsvd(:,kk)/sum(x_T1_1D_matsvd(:,kk));
        end
        
        
        %
        clear InvZeroCross
        for ii=1:Ncur
            % InvZeroCrossIm(row_pix(ii),col_pix(ii))=TI(find(~ismember(TI,TIR_1D(ii,:))));
            InvZeroCross(ii)=TI(find(~ismember(TI,TIR_1D(ii,:))));
            
        end
        
        fprintf(['Done: PsiT1, sample ',num2str(kkAll),'/',num2str(length(sampleNum)),' slice ',num2str(iiS),' \n'])
        mSave.Psi_T1(kkAll,iiS) = Psi_T1(kkAll,iiS);
        
        
        %D 1D
        ND_1D=50;
        
        D_1DA=logspace(log10(1e-13),log10(5e-9),ND_1D);
        
        E_mat_D=ones(length(b_1D),ND_1D);
        for i=1:length(b_1D)
            for j=1:ND_1D
                E_mat_D(i,j)=kD(b_1D(i),D_1DA(j));
            end
        end
        
        clear E_D
        for j=1:Ncur
            E_D(:,j)=squeeze(ImageIsoDcur(:,j))/squeeze(ImRefCurD(j));
        end
        
        
        
        nD_1D=ND_1D;
        
        [UD, SD, VD]=svd(E_mat_D);
        sD=find(diag(SD)<0.01,1)-1;
        UtilD=UD(:,1:sD);
        VtilD=VD(:,1:sD);
        StilD=SD(1:sD,1:sD);
        E_mat_Dtil=StilD*VtilD';
        
        clear E_Dtil
        for jj=1:Ncur
            E_Dtil(:,jj)=UtilD'*E_D(:,jj);
        end
        
        % 1D-D
        
        
        sVecD=diag(StilD);
        
        clear res_D_1Dsvd norm_D_1Dsvd x_D_1D_matsvd res_D_1Dtmp norm_D_1Dtmp x_D_1D_mattmp res_D_1DRefsvdIter norm_D_1DRefsvdIter x_D_1D_matRefsvdIter ...
            alphaDOpt GCVfinalD
        
        tic
        parfor kk=1:Ncur
            
            E_sig_sub=E_Dtil(:,kk);
            GCViterPre=1000;
            alphaDiter=sum(sVecD.^2)/sD/2;
            alphaDiterPre=alphaDiter/2;
            
            criter=1;
            kkcount=1;
            while (criter>=0 && kkcount<50)
                % parfor kk=1:length(alpha_D)
                [res_D_1DRefsvdIter, norm_D_1DRefsvdIter, x_D_1D_matRefsvdIter,~]= opt_cvx_1D_alt2new(E_mat_Dtil,E_sig_sub,alphaDiter,nD_1D);
                
                cVect=(E_sig_sub-E_mat_Dtil*x_D_1D_matRefsvdIter')./alphaDiter;
                
                mfree=sum((sVecD.^2)./(sVecD.^2 + alphaDiter));
                mfreePrime=-sum(sVecD.^2./((sVecD.^2 + alphaDiter).^2));
                
                
                tHat=sum(((sVecD.*cVect).^2)./(sVecD.^2 + alphaDiter));
                
                % GCViter=sD.*(res_D_1DRefsvdIter.^2)./((sD-mfree).^2);
                
                GCViterUpdate=sD.*alphaDiter^2.*norm(cVect,2)^2./((sD-mfree).^2);
                
                alphaDiterUpdate=-(alphaDiter.^2.*norm(cVect,2).^2*mfreePrime)./(tHat.*(sD-mfree));
                
                if kkcount==1
                    criter=1;
                else
                    criter=(GCViterUpdate-GCViterPre)/(alphaDiter-alphaDiterPre);
                end
                alphaDiterPre=alphaDiter;
                alphaDiter=alphaDiterUpdate;
                GCViterPre=GCViterUpdate;
                % GCViter=GCViterPre;
                % clear alphaDiterUpdate GCViter
                
                kkcount=kkcount+1;
                
            end
            
            alphaDOpt(kk)=alphaDiterPre;
            GCVfinalD(kk)=GCViterPre;
            res_D_1Dsvd(kk)=res_D_1DRefsvdIter;
            norm_D_1Dsvd(kk)=norm_D_1DRefsvdIter;
            x_D_1D_matsvd(:,kk)=x_D_1D_matRefsvdIter;
        end
        toc
        
        
        %
        for kk=1:Ncur
            Psi_D{kkAll,iiS}(kk,:)=x_D_1D_matsvd(:,kk)/sum(x_D_1D_matsvd(:,kk));
        end
        
        fprintf(['Done: PsiD, sample ',num2str(kkAll),'/',num2str(length(sampleNum)),' slice ',num2str(iiS),' \n'])

        mSave.Psi_D(kkAll,iiS) = Psi_D(kkAll,iiS);
        
        
        % T2 1D
        NT2_1D=50;
        T2_1DA=logspace(log10(1e-3),log10(5e-1),NT2_1D);
        TEexp=TE(2:end);
        
        clear E_T2
        for j=1:Ncur
            E_T2(:,j)=squeeze(ImT2cur(2:end,j))/squeeze(ImT2cur(1,j));
        end
        
        
        E_mat_T2=ones(length(TEexp),NT2_1D);
        for i=1:length(TEexp)
            for j=1:NT2_1D
                E_mat_T2(i,j)=kT2(TEexp(i),T2_1DA(j));
            end
        end
        
        
        nT2_1D=NT2_1D;
        
        
        [UT2, ST2, VT2]=svd(E_mat_T2);
        sT2=find(diag(ST2)<0.01,1)-1;
        UtilT2=UT2(:,1:sT2);
        VtilT2=VT2(:,1:sT2);
        StilT2=ST2(1:sT2,1:sT2);
        E_mat_T2til=StilT2*VtilT2';
        
        
        
        clear E_T2til
        for jj=1:Ncur
            E_T2til(:,jj)=UtilT2'*E_T2(:,jj);
        end
        
        
        % 1D-T2
        
        
        sVecT2=diag(StilT2);
        
        clear res_T2_1Dsvd norm_T2_1Dsvd x_T2_1D_matsvd res_T2_1Dtmp norm_T2_1Dtmp x_T2_1D_mattmp res_T2_1DRefsvdIter norm_T2_1DRefsvdIter x_T2_1D_matRefsvdIter ...
            alphaT2Opt GCVfinalT2
        
        tic
        parfor kk=1:Ncur
            
            E_sig_sub=E_T2til(:,kk);
            GCViterPre=1000;
            alphaT2iter=sum(sVecT2.^2)/sT2;
            alphaT2iterPre=alphaT2iter;
            
            criter=1;
            kkcount=1;
            while (criter>=0 && kkcount<50)
                % parfor kk=1:length(alpha_D)
                [res_T2_1DRefsvdIter, norm_T2_1DRefsvdIter, x_T2_1D_matRefsvdIter,~]= opt_cvx_1D_alt2new(E_mat_T2til,E_sig_sub,alphaT2iter,nT2_1D);
                
                cVect=(E_sig_sub-E_mat_T2til*x_T2_1D_matRefsvdIter')./alphaT2iter;
                
                mfree=sum((sVecT2.^2)./(sVecT2.^2 + alphaT2iter));
                mfreePrime=-sum(sVecT2.^2./((sVecT2.^2 + alphaT2iter).^2));
                
                
                tHat=sum(((sVecT2.*cVect).^2)./(sVecT2.^2 + alphaT2iter));
                
                % GCViter=sD.*(res_D_1DRefsvdIter.^2)./((sD-mfree).^2);
                
                GCViterUpdate=sT2.*alphaT2iter^2.*norm(cVect,2)^2./((sT2-mfree).^2);
                
                alphaT2iterUpdate=-(alphaT2iter.^2.*norm(cVect,2).^2*mfreePrime)./(tHat.*(sT2-mfree));
                
                if kkcount==1
                    criter=1;
                else
                    criter=(GCViterUpdate-GCViterPre)/(alphaT2iter-alphaT2iterPre);
                end
                alphaT2iterPre=alphaT2iter;
                alphaT2iter=alphaT2iterUpdate;
                GCViterPre=GCViterUpdate;
                % GCViter=GCViterPre;
                % clear alphaDiterUpdate GCViter
                
                kkcount=kkcount+1;
                
            end
            alphaT2Opt(kk)=alphaT2iterPre;
            GCVfinalT2(kk)=GCViterPre;
            res_T2_1Dsvd(kk)=res_T2_1DRefsvdIter;
            norm_T2_1Dsvd(kk)=norm_T2_1DRefsvdIter;
            x_T2_1D_matsvd(:,kk)=x_T2_1D_matRefsvdIter;
        end
        toc
        
        for kk=1:Ncur
            Psi_T2{kkAll,iiS}(kk,:)=x_T2_1D_matsvd(:,kk)/sum(x_T2_1D_matsvd(:,kk));
        end
        fprintf(['Done: PsiT2, sample ',num2str(kkAll),'/',num2str(length(sampleNum)),' slice ',num2str(iiS),' \n'])
        mSave.Psi_T2(kkAll,iiS) = Psi_T2(kkAll,iiS);
        
        
        %D-T2
        
        
        ImRefCurDT2=ImRefCurD;
        
        
        ImageIsoDT2cur=ImageIsoDT2./permute(repmat(ImRefCurDT2,[length(TE_DT2_M3) ,1,length(b_DT2)]),[1,3,2]);
        
        
        alpha2Dpre=logspace( -4, 4, 7 );
        
        alpha2D=alpha2Dpre(3);
        
        E_mat_D_2D=ones(length(b_DT2),ND_1D);
        for i=1:length(b_DT2)
            for j=1:ND_1D
                E_mat_D_2D(i,j)=kD(b_DT2(i),D_1DA(j));
            end
        end
        
        
        E_mat_T2_2D=ones(length(TE_DT2_M3),length(T2_1DA));
        for i=1:length(TE_DT2_M3)
            for j=1:length(T2_1DA)
                
                E_mat_T2_2D(i,j)=kT2(TE_DT2_M3(i),T2_1DA(j));
                
                
            end
        end
        
        
        %
        
        clear res_T2D norm_x_T2D x_matCon_T2D res_T2Ddiff norm_x_T2Ddiff x_matCon_T2Ddiff
        clear res_T2Dtmp1 norm_x_T2Dtmp1 x_matCon_T2Dtmp1 res_T2Dtmp2 norm_x_T2Dtmp2 x_matCon_T2Dtmp2
        tic
        parfor kk=1:Ncur
            
            E_sig_mat=squeeze(ImageIsoDT2cur(:,:,kk));
            [res_T2Dtmp1(:,kk), norm_x_T2Dtmp1(:,kk), x_matCon_T2Dtmp1(:,kk,:,:)]= opt_cvx2dSVD_MADCO_Brain_a(E_mat_D_2D,E_mat_T2_2D,E_sig_mat,alpha2D,Psi_D{kkAll,iiS}(kk,:),Psi_T2{kkAll,iiS}(kk,:));
            %         [res_T2Dtmp2(:,kk), norm_x_T2Dtmp2(:,kk), x_matCon_T2Dtmp2(:,kk,:,:)]= opt_cvx2dSVD_MADCO_Brain_a(E_mat_D_2D,E_mat_T2_2D,E_sig_mat,alpha2D,Psi_D{ii}(kk,:),Psi_T2diff{ii}(kk,:));
            
            %         parsave('outputLoop',[ii kk]);
            
        end
        res_T2D=res_T2Dtmp1;
        norm_x_T2D=norm_x_T2Dtmp1;
        x_matCon_T2D=x_matCon_T2Dtmp1;
        
        %     res_T2Ddiff{ii}=res_T2Dtmp2;
        %     norm_x_T2Ddiff{ii}=norm_x_T2Dtmp2;
        %     x_matCon_T2Ddiff{ii}=x_matCon_T2Dtmp2;
        
        toc
        
        
        for ii=1:Ncur
            Psi_T2D{kkAll,iiS}(:,:,ii)=squeeze(x_matCon_T2D(1,ii,:,:))/sum(sum(squeeze(x_matCon_T2D(1,ii,:,:)),1),2);
        end
        
        Psi_T2D{kkAll,iiS}=permute(Psi_T2D{kkAll,iiS},[2 1 3]);
        
        fprintf(['Done: PsiT2D, sample ',num2str(kkAll),'/12, slice ',num2str(iiS),' \n'])
        mSave.Psi_T2D(kkAll,iiS) = Psi_T2D(kkAll,iiS);
        
        
        
        
        
        % %
        %D-T1
       
        ImRefCurDT1=ImRefCurD;
        
        ImageIsoDT1cur=ImageIsoDT1;
        
        [~, indDfor2D]=intersect(b_1D,b_DT1);% the corresponding b-values indices from the 1D diffusion
        ImRefCurDT1_TI=ImageIsoDcur(indDfor2D,:);
        
        
        
        clear  ImDT1curMod TIR2D
        
        
        indT1D_1=logical(round(InvZeroCross,3)<round(TI_DT1_M3(1),3)); %use all TIs
        indT1D_2=logical(round(InvZeroCross,3)==round(TI_DT1_M3(1),3)); %discard 1st TI
        indT1D_3=logical(round(InvZeroCross,3)<=round(TI_DT1_M3(2),3) & round(InvZeroCross,3)>round(TI_DT1_M3(1),3)); % flip sign of 1st TI, discard 2nd TI
        indT1D_4=logical(round(InvZeroCross,3)<round(TI_DT1_M3(3),3) & round(InvZeroCross,3)>round(TI_DT1_M3(2),3)); %use all TIs, flip sign of 1st and 2nd TIs
        indT1D_5=logical(round(InvZeroCross,3)==round(TI_DT1_M3(3),3)); %discard 3rd TI, flip sign
        indT1D_6=logical(round(InvZeroCross,3)<round(TI_DT1_M3(4),3) & round(InvZeroCross,3)>round(TI_DT1_M3(3),3)); %use all TIs, flip sign of 1-3 TIs
        indT1D_7=logical(round(InvZeroCross,3)>round(TI_DT1_M3(4),3)); %keep all, flip sign
        %     indT1D_6=logical(round(InvZeroCross{jj},3)==round(TI_DT1_M3(4),3)); %discard 4rd TI, flip sign
        indT1D_8=logical((indT1D_1+indT1D_2+indT1D_3+indT1D_4+indT1D_5+indT1D_6)==0); %all the rest
        %     indT1D_7=logical((indT1D_1+indT1D_2+indT1D_3+indT1D_4+indT1D_5)==0); %all the rest
        
        
        %     length(find(indT1D_1==0))-length(find(indT1D_2==1))-length(find(indT1D_3==1))-length(find(indT1D_4==1))-length(find(indT1D_5==1))
        %
        %     find((indT1D_1+indT1D_2+indT1D_3+indT1D_4+indT1D_5+indT1D_6+indT1D_7)==0)
        
        
        for ii=1:Ncur
            
            if indT1D_1(ii)
                TIR2D{ii}=TI_DT1_M3(1:end-1);
                clear ImDT1curModtemp
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,1:(end-1),ii));
                    tmp2a=tmp1;
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    
                    
                    ImDT1curModtemp(kk,:)=tmp3/ImRefCurD(ii)/2;
                    %                     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            
            
            
            if indT1D_2(ii)
                TIR2D{ii}=TI_DT1_M3(2:end);
                clear ImDT1curModtemp
                
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,:,ii));
                    tmp1(1)=[];
                    tmp2a=tmp1;
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            
            if indT1D_3(ii)
                TIR2D{ii}=TI_DT1_M3([1,3:end]);
                clear ImDT1curModtemp
                
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,:,ii));
                    tmp1(2)=[];
                    tmp2a=tmp1;
                    tmp2a(1)=-tmp1(1);
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            
            if indT1D_4(ii)
                TIR2D{ii}=TI_DT1_M3(1:end);
                clear ImDT1curModtemp
                
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,1:(end),ii));
                    tmp2a=tmp1;
                    tmp2a(1:2)=-tmp1(1:2);
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            if indT1D_5(ii)
                TIR2D{ii}=TI_DT1_M3([1:2,4]);
                clear ImDT1curModtemp
                
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,:,ii));
                    tmp1(3)=[];
                    tmp2a=tmp1;
                    tmp2a(1:2)=-tmp1(1:2);
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            if indT1D_6(ii)
                TIR2D{ii}=TI_DT1_M3(1:end);
                clear ImDT1curModtemp
                
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,1:(end),ii));
                    tmp2a=tmp1;
                    tmp2a(1:3)=-tmp1(1:3);
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                     ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            if indT1D_7(ii)
                TIR2D{ii}=TI_DT1_M3(1:end);
                clear ImDT1curModtemp
                
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,:,ii));
                    tmp2a=-tmp1;
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                      ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
            if indT1D_8(ii)
                TIR2D{ii}=TI_DT1_M3(1:end);
                clear ImDT1curModtemp
                for kk=1:size(ImageIsoDT1cur,1)
                    
                    tmp1=squeeze(ImageIsoDT1cur(kk,:,ii));
                    tmp2a=tmp1;
                    
                    tmp3=repmat(ImRefCurDT1_TI(kk,ii),1,size(tmp2a,2)) - tmp2a;
                    ImDT1curModtemp(kk,:)=tmp3/squeeze(ImRefCurD(ii))/2;
                    %                      ImT1curMod2(:,row_pix(j),col_pix(j))=tmp3/ImT1cur(end,row_pix(j),col_pix(j))/2;
                end
                ImDT1curMod{ii}=ImDT1curModtemp;
            end
            
        end
        
        
        
        
        alpha2Dpre=logspace( -4, 4, 7 );
        
        
        alpha2D=alpha2Dpre(3);
        
        E_mat_D_2D=ones(length(b_DT1),ND_1D);
        for i=1:length(b_DT1)
            for j=1:ND_1D
                E_mat_D_2D(i,j)=kD(b_DT1(i),D_1DA(j));
            end
        end
        
        
        clear res_T1D norm_x_T1D x_matCon_T1D res_T1Dtmp norm_x_T1Dtmp x_matCon_T1Dtmp
        tic
        parfor kk=1:Ncur
            TIR2Dlocaliter=TIR2D{kk};
            E_mat_T1_2D=ones(size(TIR2Dlocaliter,2),length(T1_1DA));
            
            for i=1:size(TIR2Dlocaliter,2)
                for j=1:length(T1_1DA)
                    E_mat_T1_2D(i,j)=kT1(TIR2Dlocaliter(i),T1_1DA(j));
                end
            end
            
            E_sig_mat=(ImDT1curMod{kk});
            [res_T1Dtmp(:,kk), norm_x_T1Dtmp(:,kk), x_matCon_T1Dtmp(:,kk,:,:)]= opt_cvx2dSVD_MADCO_Brain_a(E_mat_D_2D,E_mat_T1_2D,E_sig_mat,alpha2D,Psi_D{kkAll,iiS}(kk,:),Psi_T1{kkAll,iiS}(kk,:));
            
            
        end
        
        res_T1D=res_T1Dtmp;
        norm_x_T1D=norm_x_T1Dtmp;
        x_matCon_T1D=x_matCon_T1Dtmp;
        
        toc
        
        
        
        for ii=1:Ncur
            Psi_T1D{kkAll,iiS}(:,:,ii)=squeeze(x_matCon_T1D(1,ii,:,:))/sum(sum(squeeze(x_matCon_T1D(1,ii,:,:)),1),2);
        end
        
        
        Psi_T1D{kkAll,iiS}=permute(Psi_T1D{kkAll,iiS},[2 1 3]);
        
        fprintf(['Done: PsiT1D, sample ',num2str(kkAll),'/',num2str(length(sampleNum)),' slice ',num2str(iiS),' \n'])
        mSave.Psi_T1D(kkAll,iiS) = Psi_T1D(kkAll,iiS);
        
     
        
        %2D T1-T2 vector notation
        
        
        
        ImageT1T2Cur=ImageT1T2;
        % ImageT1T2matCur=permute(ImageT1T2matCur,[1 2 4 3 5]);
        [~, TEind2D]=intersect(round(1000*TE),round(1000*TE_T1T2));
        ImT1T2curRef=ImT2cur(TEind2D,:);
        %
        
        clear  ImT1T2curMod TIR2Da TE2Da
        
        for ii=1:Ncur
            kk1=1;
            for kk=1:length(TI_T1T2)
                if TI_T1T2(kk)<InvZeroCross(ii)
                    
                    TIR2Da{ii}(kk1)=TI_T1T2(kk);
                    TE2Da{ii}(kk1)=TE_T1T2(kk);
                    
                    tmp1=-ImageT1T2Cur(kk,ii);
                    [~,TEind2Dcur]=intersect(round(1000*TE),round(1000*TE_T1T2(kk)));
                    tmp2=ImT2cur(TEind2Dcur,ii) - tmp1;
                    ImT1T2curMod{ii}(kk1)=tmp2/squeeze(ImT2cur(1,ii))/2;
                    kk1=kk1+1;
                end
                
                if TI_T1T2(kk)>InvZeroCross(ii)
                    
                    TIR2Da{ii}(kk1)=TI_T1T2(kk);
                    TE2Da{ii}(kk1)=TE_T1T2(kk);
                    
                    tmp1=ImageT1T2Cur(kk,ii);
                    [~,TEind2Dcur]=intersect(round(1000*TE),round(1000*TE_T1T2(kk)));
                    tmp2=ImT2cur(TEind2Dcur,ii) - tmp1;
                    ImT1T2curMod{ii}(kk1)=tmp2/squeeze(ImT2cur(1,ii))/2;
                    kk1=kk1+1;
                end
                
            end
        end
        
        
        
        
        
        % 2D T1-T2 - with full D grid
        alpha2Dpre=logspace( -4, 4, 7 );
        
        
        alpha2D=alpha2Dpre(3);
        
        
        %
        
        clear res_T1T2 norm_x_T1T2 x_matCon_T1T2 res_T1T2tmp norm_x_T1T2tmp x_matCon_T1T2tmp
        
        tic
        parfor kk=1:Ncur
            TIR2Dlocaliter=TIR2Da{kk};
            TE2Dlocaliter=TE2Da{kk};
            
            E_mat_T1T2=ones(size(TIR2Dlocaliter,2),length(T1_1DA)*length(T2_1DA));
            
            for i1=1:length(TIR2Dlocaliter)
                
                kk2=1;
                for j2=1:length(T1_1DA)
                    for j3=1:length(T2_1DA)
                        E_mat_T1T2(i1,kk2)=kT1(TIR2Dlocaliter(i1),T1_1DA(j2)).*kT2(TE2Dlocaliter(i1),T2_1DA(j3));
                        kk2=kk2+1;
                    end
                    
                end
                
            end
            
            
            
            E_sigCur=ImT1T2curMod{kk}';
            [res_T1T2tmp(kk), norm_x_T1T2tmp(kk), x_matCon_T1T2tmp(:,kk)]= opt_cvxMADCO_2D(E_mat_T1T2,E_sigCur,alpha2D,Psi_T1{kkAll,iiS}(kk,:),Psi_T2{kkAll,iiS}(kk,:))
            
            
        end
        
        res_T1T2=res_T1T2tmp;
        norm_x_T1T2=norm_x_T1T2tmp;
        x_matCon_T1T2=x_matCon_T1T2tmp;
        
        toc
        
        for kk=1:Ncur
            
            tmp=reshape(x_matCon_T1T2(:,kk),length(T2_1DA),length(T1_1DA)); % !!P(T2,T1)!!
            Psi_T1T2{kkAll,iiS}(:,:,kk) = tmp/sum(tmp(:));
            
        end
        
        
        
        Psi_T1T2{kkAll,iiS}=permute(Psi_T1T2{kkAll,iiS},[2 1 3]); % !!P(T1,T2)!!
        
        fprintf(['Done: PsiT1T2, sample ',num2str(kkAll),'/',num2str(length(sampleNum)),' slice ',num2str(iiS),' \n'])
        mSave.Psi_T1T2(kkAll,iiS) = Psi_T1T2(kkAll,iiS);

        
    end
end


