function SNPdosagedata = SNPdata2table(snpdata_converted)

SNPdosagedata = array2table(snpdata_converted.Data.DS);
SNPdosagedata.Properties.RowNames = snpdata_converted.SampleInfo(:,1);

for i = 1:size(snpdata_converted.Data.DS,2)
    SNPdosagedata.Properties.VariableNames{:,i} = ['Ch',snpdata_converted.VarInfo.CHROM{i},'_Pos',snpdata_converted.VarInfo.POS{i},'_',snpdata_converted.VarInfo.REF{i},'_',snpdata_converted.VarInfo.ALT{i}];
    SNPdosagedata.Properties.VariableUnits{:,i} = snpdata_converted.VarInfo.INFO{i};
end

end