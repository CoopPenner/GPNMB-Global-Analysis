%quick test of a few hits

origDataRaw= readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 1/20210120_olinkdata.csv');

repDataRaw=readtable('/Users/pennerc/Desktop/2_Cognitive Decline Project/data sets/olink data 2/20220726_olinkreplication.csv');

finnGenData=readtable('/Volumes/PC60/PDGWASStuff/GWASCompleteMatchFinnGen500kMBEvalFinal.csv');

ukBioBankData=readtable('/Volumes/PC60/PDGWASStuff/GWASCompleteMatchUKBioBankData500KMBWindow_Final.csv');



GWASData2use=finnGenData;
olinkData2use=repDataRaw;

easyAccessData=olinkData2use(:,12:end); % I know I will make a stupid math mistake unless I knockout the non genes for indexing lol 
origNames=olinkData2use(:,12:end).Properties.VariableNames;
 %origplate=olinkData2use.PlateID;
 %origBatch=olinkData2use.Batch;
 %daysPost=olinkData2use.DaysBetween;
dx=olinkData2use.Group;
age=olinkData2use.AGE__YRS;
%tubeNum=olinkData2use.TubeNumber;
namesAtPlay=unique(GWASData2use.matched_gene);
sexAtPlay=olinkData2use.Sex;
INDDID=olinkData2use.INDDID;
sampleDate=olinkData2use.SampleDate;
% Identify groups
isHC = strcmp(dx, 'NC');
isPD = strcmp(dx, 'PD');


%% finding what of our gwas hits are in the olink data to a reasonable
% extent 

genePres=false(1,length(namesAtPlay));
        for dd=1:length(namesAtPlay)
            if sum(strcmp(origNames,namesAtPlay{dd} )   )>0
                origLoc = find(strcmp(origNames, namesAtPlay{dd}));
                arrayAtPlay=table2array(easyAccessData(:,origLoc));
    
                    if iscell(arrayAtPlay(1))
                    nanValues=sum(contains(arrayAtPlay,'NA')); %accounting for inconsitent notation for nans between replicates
                    else
                        nanValues=sum(isnan(arrayAtPlay));
                    end
    
                        if sum(nanValues)< 50
                    genePres(dd)=  sum(strcmp(origNames,namesAtPlay{dd} )   )>0; 
                        end
    
            end
        end

%% separate nan filter catch irritated by this

nanPres=false(1,width(easyAccessData));
        for dd=1:width(easyAccessData)

                arrayAtPlay=table2array(easyAccessData(:,dd));
    
                    if iscell(arrayAtPlay(1))
                    nanValues=sum(contains(arrayAtPlay,'NA')); %accounting for inconsitent notation for nans between replicates
                    else
                        nanValues=sum(isnan(arrayAtPlay));
                    end
    
                    nanPres(dd)=  sum(nanValues)  ; 
   
            end
        

            normFactor=mean(table2array(easyAccessData(:,~nanPres))'  );


usableGenes=namesAtPlay(genePres);


%%        quick and dirty glme to look at diff between PD and HC

filter2use=strcmp(dx,'NA') | cellfun(@isempty, dx)| strcmp(dx,'AD');

[qValsTrue,PValsTrue] = quickAnovaForGWASAnalysis(usableGenes,easyAccessData,~filter2use,genePres,origNames,dx,age,sexAtPlay,INDDID,normFactor);

% plotting significant hits
        
sigIdx=find(qValsTrue<.05 ) ;

colorMap = struct('HC', [0.13 0.55 0.13], 'PD', [0.6 0.2 0.8]   ); 

  for dd= 1:length(sigIdx)
        
             origLoc = find(strcmp(origNames, usableGenes{sigIdx(dd)}  )  );
             val2use=table2array(easyAccessData(:,origLoc)); val2use(val2use<0)=0;
             protAtPlay=usableGenes{sigIdx(dd)};


            % Make plot
            a=figure; hold on;
            
            % HC Bar
            b = bar(1, nanmedian(val2use(isHC & ~filter2use )), 'FaceColor', colorMap.HC, 'FaceAlpha', 0.3, 'BarWidth', 1.2);
            scatter(rand(1, sum(isHC& ~filter2use  )) * 0.6 + 0.7, val2use(isHC & ~filter2use  ), ...
                'o', 'MarkerFaceColor', colorMap.HC, 'MarkerEdgeColor', colorMap.HC);
            
            % WT Bar
            b = bar(3, nanmedian(val2use(isPD & ~filter2use  )  ), 'FaceColor', colorMap.PD, 'FaceAlpha', 0.3, 'BarWidth', 1.2);
            scatter(rand(3, sum(isPD & ~filter2use  )) * 0.6 + 2.7, val2use(isPD& ~filter2use  ), ...
                'o', 'MarkerFaceColor', colorMap.PD, 'MarkerEdgeColor', colorMap.PD);
            
            % Format
            xlim([0 4])
            set(gca, 'XTick', [1 3], 'XTickLabel', {'HC', 'PD'},'FontSize',15)
            ylabel(' Olink Protein Level','FontSize',20)
            box on
                 
            title([protAtPlay, ' Protein Level Comparison PD vs NC q=', num2str(qValsTrue(sigIdx(dd)))])

           

            saveas(a, ['/Volumes/PC60/PDGWASStuff/',protAtPlay,'.jpg'])
            
  end

%%

refDataTable = readtable('/Volumes/PC60/InqueryDatasets/FinalizedSets/PD_Broad_Overview_data_including_SNP.xlsx');

[resultsSlope,resultsTot] = runProteinCognitiveModels(usableGenes, easyAccessData, genePres, origNames, refDataTable, 1, 'slope',INDDID,normFactor,sampleDate)



%% Now finally a permutation test to see the likelihood of this many sig hits in a given group

rng(42);
permNum=1000;
filter2use=strcmp(dx,'NA') | cellfun(@isempty, dx)| strcmp(dx,'AD');


qCollect=nan(1,permNum);
bonfCollect=qCollect;
for qq=1:permNum

randVec=randperm(sum(~nanPres)); randVec=randVec(1:sum(genePres));

usableGenes=origNames(~nanPres); usableGenes=usableGenes(randVec);
[qValsPerm,pValsPerm] = quickAnovaForGWASAnalysis(usableGenes,easyAccessData,~filter2use,genePres,origNames,dx,age,sexAtPlay,daysPost,INDDID,normFactor);

qCollect(qq)=sum(  qValsPerm <  .05);
bonfCollect(qq)= sum(pValsPerm < (.05/length(pValsPerm))) ;


end


figure

histogram(bonfCollect)


%%

% 
% 
% intraSamplePlate1_1a=173; % index for RBM301-01 pooled plasma
% intraSamplePlate1_1b=174; % index for RBM301-01 
% 
% intraSamplePlate1_2a = 87; %index for 143922
% intraSamplePlate1_2b = 88; %index for 143923
% 
% inter_samplePlate1_Plate3_1a = 175; %index for RBM302-01 
% inter_samplePlate1_Plate3_1b  = 176; %index for RBM302-02 
% 
% 
% inter_samplePlate1_Plate3_2a = 77; %index for 142965
% inter_samplePlate1_Plate3_2b = 78; %index for 142966
% 
% intraSamplePlate3a=182;
% intraSamplePlate3b=183;
% 
% intraSamplePlate2_1a=179;
% intraSamplePlate2_1b=180;
% 
% allSamp=origDataRaw.SampleID;
% 


% % LACTB=origDataRaw.LACTB2;
% SCARB2=origDataRaw.SCARB2;
% MSR1=origDataRaw.MSR1;
% % CTSD=origDataRaw.CTSD;
% % CSTB=origDataRaw.CSTB;
% DAB2=origDataRaw.DAB2;
% FBP1=origDataRaw.FBP1;
% FCGR2A=origDataRaw.FCGR2A;
% FES=origDataRaw.FES;
% HLAE=origDataRaw.HLA_E;
% %ITGB1=origDataRaw.ITGB1 These only present in one GWAS
% % LRP1=origDataRaw.LRP1
% % ANXA1=origDataRaw.ANXA1;
% % LTA4H=origDataRaw.LTA4H;
% MERTK=origDataRaw.MERTK;
% MIF=origDataRaw.MIF;
% NCF2=origDataRaw.NCF2;
% % NRP1=origDataRaw.NRP1;
% % PLAG27
% %TYMP
% SLAMF8=origDataRaw.SLAMF8;
% TGFBI=origDataRaw.TGFBI;
% % FBP1 LAP3 LSP1 TNFAIP8
% TESTR=origDataRaw.USP