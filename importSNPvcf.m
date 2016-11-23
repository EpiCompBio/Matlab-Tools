function data_converted = importSNPvcf(filename,headerline)

if nargin < 2
    headerline = 36;
end

temp1 = importdata(filename,'\n');
temp1 = temp1(headerline+1:end,1);
for i = 1:length(temp1)
    temp2(i,:)=strsplit(temp1{i,1},'\t');
end
data_imported = temp2';
clear temp*;

noVarFields = 9;
VarInfo_fieldnames = strrep(data_imported(1:noVarFields,1),'#','');
for i = 1:noVarFields
    data_converted.VarInfo.(VarInfo_fieldnames{i}) = data_imported(i,2:end);
end

noOfSamples = size(data_imported(noVarFields+1:end,1),1);
idx_replicate = false(size(data_imported(noVarFields+1:end,1),1),1);
for i = 1:noOfSamples
    temp1 = strsplit(data_imported{noVarFields+i,1},'_');
    if length(temp1) > 2
        temp2{1,1} = strcat(temp1{1,1},'_',temp1{1,2});
        temp2{1,2} = strcat(temp1{1,3},'_',temp1{1,4});
        idx_replicate(i) = 1;
    else
        temp2 = temp1;
    end        
        data_converted.SampleInfo(i,:) = temp2;
        clear temp*;
end

for i = 1:noOfSamples
    for j = 1:size(data_imported,2)-1
        temp = strsplit(data_imported{i+noVarFields,j+1},':');
        data_converted.Data.GT(i,j) = temp(1);
        data_converted.Data.DS(i,j) = str2double(temp{2});
        data_converted.Data.GP(i,j) = temp(3);
    end
end


% delete replicates
data_converted.SampleInfo(idx_replicate,:) = [];
data_converted.Data.GT(idx_replicate,:) = [];
data_converted.Data.DS(idx_replicate,:) = [];
data_converted.Data.GP(idx_replicate,:) = [];

end