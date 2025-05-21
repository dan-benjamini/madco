function [ImageRef,Image, bvalOreintAvg,TE ]=read_2dseq_3D_DiffIsoAcqPath(im_ind, pathname)
% 
% close all
% clear all
% clc
% % 
% im_ind=20;


% pathname=uigetdir('/Users/danbenj/Documents/RESEARCH/Experiments/','Pick Data Directory');
 
curdir=strcat(pathname,'/',num2str(im_ind(1)),'/pdata/1/2dseq');
    
    version_file = [pathname, '/',num2str(im_ind(1)), '/acqp'];
    version_hand=fopen(version_file);
    version=fscanf(version_hand,'%c');
    version_place = findstr(version, '##$BYTORDA=');
    tmp=textscan(version(version_place+11:version_place+17), '%s');
    version=tmp{1}{1};
    if findstr(version,'little')
        reco=fopen(curdir,'r','ieee-le');        % little endian
    else
        reco=fopen(curdir,'r');                  % big endian
    end
    
[b_cell, ~]=find_multi_method(pathname,im_ind(1),'PVM_DwEffBval');
bval(:,1)=cellfun(@str2num,[b_cell{:}]);
[bPers_cell, ~]=find_multi_method(pathname,im_ind(1),'PVM_DwBvalEach');
bPresval(:,1)=cellfun(@str2num,[bPers_cell{:}]);


[TE(1), ~]=find_multi_method_mod(pathname,im_ind(1),'EchoTime');
[Dir_cell, ~]=find_multi_method_mod(pathname,im_ind(1),'PVM_DwDir');
[BmatVectCell, ~]=find_multi_method_mod(pathname,im_ind(1),'PVM_DwBMat');

[dims_cell,num_dim]=find_multi_reco(pathname,im_ind(1),'RECO_size',1);
[slope_cell, num_grad]=find_multi_reco(pathname,im_ind(1),'RECO_map_slope',1);


BmatVect=cellfun(@str2num,[BmatVectCell{:}]);
Bmat=reshape(BmatVect,[3,3,length(BmatVect)/9]);

dims=cellfun(@str2num,[dims_cell{:}]);
slope=cellfun(@str2num,[slope_cell{:}]);
Dir=cellfun(@str2num,[Dir_cell{:}]);
num_dir=length(Dir)/3;

A0=0;
if (length(bval)/num_dir)>length(bPresval)
    A0=1;
end


A=fread(reco,'int16');

A1=permute(reshape(A,[dims(2),dims(1),dims(3),num_grad]),[2 1 3 4]);

if A0    
    ImageRef=squeeze(A1(:,:,:,1))/slope(1);
    for ii=2:length(slope)
        ImagePre(:,:,:,ii-1)=squeeze(A1(:,:,:,ii))/slope(ii);
    end
else
    ImageRef=zeros(size(squeeze(A1(:,:,:,1))));
    for ii=1:length(slope)
        ImagePre(:,:,:,ii)=squeeze(A1(:,:,:,ii))/slope(ii);
    end
end

bvalPerNum=length(bPresval);
bvalOreintAvgMat=zeros(bvalPerNum,num_dir);

for i=1:num_dir
    bvalOreintAvgMat(:,i)=bval((bvalPerNum*(i-1)+2):(i*bvalPerNum)+1);
    Image(:,:,:,:,i)=ImagePre(:,:,:,(bvalPerNum*(i-1)+1):i*bvalPerNum);

end

bvalOreintAvg=mean(bvalOreintAvgMat,2);
TE=TE/1000;




