function [Image,TI]=read_2dseq_3D_TIpath(im_ind,pathname)

% close all
% clear all
% clc
% % 
% im_ind=49:68;
% sliceNum=16;
% pathname=uigetdir('/Users/danbenj/Documents/RESEARCH/Experiments/','Pick Data Directory');


for i=1:length(im_ind)

 
curdir=strcat(pathname,'/',num2str(im_ind(i)),'/pdata/1/2dseq');
    
    version_file = [pathname, '/',num2str(im_ind(i)), '/acqp'];
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
    
[TI(i), ~]=find_multi_method_mod(pathname,im_ind(i),'PVM_SelIrInvTime');

[TR(i), ~]=find_multi_method_mod(pathname,im_ind(i),'PVM_RepetitionTime');

[dims_cell,num_dim]=find_multi_reco(pathname,im_ind(i),'RECO_size',1);
[slope_cell, num_grad]=find_multi_reco(pathname,im_ind(i),'RECO_map_slope',1);

dims=cellfun(@str2num,[dims_cell{:}]);
slope=cellfun(@str2num,[slope_cell{:}]);


A=fread(reco,'int16');

A1(i,:,:,:)=permute(reshape(A,[dims(2),dims(1),dims(3)]),[2 1 3]);



if i==1
    Image=zeros(length(im_ind),dims(1),dims(2),dims(3));
end



 
Image(i,:,:,:)=squeeze(A1(i,:,:,:))/slope;
    
    
fclose(reco);
fclose(version_hand);
    clear A
    clear A1
    clear reco

end

TI=TI/1000;
TR=TR/1000;


 
