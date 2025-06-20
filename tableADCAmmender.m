% simple function to ammend our previous ADC dataset with new SNP statuses 

addpath('/Users/pennerc/Documents/Documents - BGS-STU23-138/MATLAB/GPNMB_Global_Analysis/GPNMB-Global-Analysis');

data2ammend=readtable('/Users/pennerc/Documents/Documents - BGS-STU23-138/updatedADCData2Ammend.xlsx');

origData=readtable('/Users/pennerc/Downloads/InqueryDatasets/FinalizedSets/ADC_MMSE_TotalPatients.xlsx');

        
        
        for dd=1:height(data2ammend)
        ptAtPlay=data2ammend.INDDID(dd);

        mainLoc=find(origData.INDDID==ptAtPlay);
         oldSNP=origData.rs199347(mainLoc);
            if ~cellfun(@isempty,oldSNP)
                error('something wrong, old SNP present')

            end


        
        newSNP=[data2ammend.rs199347_a1{dd}, data2ammend.rs199347_a2{dd}];

            if sum(newSNP=='A')==1
            snp2Ammend='CT';
            elseif sum(newSNP=='A')==2
                snp2Ammend='TT';
            elseif sum(newSNP=='B')==2
                snp2Ammend='CC';
            end
            origData.rs199347(mainLoc)={snp2Ammend};
            unique(origData.ADCCDx1(mainLoc))

        end



        writetable(origData,'/Users/pennerc/Downloads/InqueryDatasets/FinalizedSets/ADC_MMSE_TotalPatients_withUpdate.xlsx')