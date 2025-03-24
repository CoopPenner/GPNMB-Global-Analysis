function [] = tdpBrainStagerGPNMB(pt2use, pathTable, disName )


%initializing all variables
snpStat=pathTable.rs199347;
overGPNMB= contains(snpStat,'TT'); %the major allele
het=contains(snpStat, 'CT');
underGPNMB=contains(snpStat,'CC'); %the minor allele
emptySNP=cellfun(@isempty,snpStat);

ageAtDeath = pathTable.AgeatDeath;
ageAtOnset = pathTable.GlobalAgeOnset; 
disDur = ageAtDeath - ageAtOnset;
bioSex = pathTable.Sex;
ptID=pathTable.INDDID;





% from the paper.... 
% 1 agranular motor cortex brainstem motor nuclei (V VII and X-XII) and SC
% We have almost no samples for Brainstem so I am excluding that



%2 PF gyrus and MF gyrus reticular formation precerebellar nuclei and red
%nucleus 
%not sure exactly where to nab precerebellar nuclei... brainstem,
%cerebellum spinal cord? PF and MF gyrus will be neocortical and MF red
%nucleus blends with the reticular formation rostrally... idk... should
%this be Pons, midbrain? it's at the level of the SN. I will take both Pons
%and Medulla for stage 2 so it will be Neocortical, Pons, Medulla.

%3 is caudate, putamen, Substantia Nigra, Mamillary Bodies

%4 is Amygdala Hippocampus EC and DG...




%I'd also like to do something a bit simpler, based on path densities it
%looks like there are two separate streams, dorsal from spine and caudal
%from mc
%have a spinal and cerebral progression scale
%Spinal: 1) SC 2) Medulla 3) Pons (simple)
%cerebral: 1) MC 2) Neocortical Cingulate 3) SN GP CP and MB 4) DG Amyg
%Hippo EC Temporal Sulcus

%% pt 1 the Brettschneider scheme

cortical1={'MC','SC'}; cortical2={'Cing','Neocortical','Pons','Med'}; cortical3={'SN','GP','CP','MB'}; cortical4={'CS','Amyg','DG','EC', 'TS'};

pathPosBlank=nan(1,height(pathTable));

pathCut=1;
[cort1Path]=TDPStager(cortical1,pathTable,pathCut,pt2use);
cort1PathPos=pathPosBlank;cort1PathPos(cort1Path>1)=true; % greater than an average of 1

[cort2Path]=TDPStager(cortical2,pathTable,pathCut,pt2use);
cort2PathPos=pathPosBlank;cort2PathPos(cort1Path>1 & cort2Path>1)=true;


[cort3Path]=TDPStager(cortical3,pathTable,pathCut,pt2use);
cort3PathPos=pathPosBlank;cort3PathPos(cort1Path>1 & cort2Path>1 & cort3Path>1)=true;


[cort4Path]=TDPStager(cortical4,pathTable,pathCut,pt2use);
cort4PathPos=pathPosBlank;cort4PathPos(cort1Path>1 & cort2Path>1 & cort3Path>1 & cort4Path>1)=true;


%now ranking path burden in patients with sufficent samples and outputting
%a simple stack bar graph 


path2use=nansum([cort1PathPos;cort2PathPos;cort3PathPos;cort4PathPos]);
remVals=~(pt2use) | emptySNP; %this includes pt2use ie non patients will be nanned out


figure



overGPNMBRat=[    sum(path2use(overGPNMB & ~remVals)==0  )/ sum((overGPNMB & ~remVals)  )*100,...
    sum(path2use(overGPNMB & ~remVals)==1  )/ sum((overGPNMB & ~remVals)  )*100,...
sum(path2use(overGPNMB & ~remVals)==2  )/ sum((overGPNMB & ~remVals) )*100,...
sum(path2use(overGPNMB & ~remVals)==3  )/ sum((overGPNMB & ~remVals) )*100,...
sum(path2use(overGPNMB & ~remVals)==4  )/ sum((overGPNMB & ~remVals) )*100];


    



hetGPNMBRat=[    sum(path2use(het & ~remVals)==0  )/ sum((het & ~remVals)  )*100,...
    sum(path2use(het & ~remVals)==1  )/ sum((het & ~remVals)  )*100,...
sum(path2use(het & ~remVals)==2  )/ sum((het & ~remVals) )*100,...
sum(path2use(het & ~remVals)==3  )/ sum((het & ~remVals) )*100,...
sum(path2use(het & ~remVals)==4  )/ sum((het & ~remVals) )*100];






underGPNMBRat=[    sum(path2use(underGPNMB & ~remVals)==0  )/ sum((underGPNMB & ~remVals)  )*100,...
    sum(path2use(underGPNMB & ~remVals)==1  )/ sum((underGPNMB & ~remVals)  )*100,...
sum(path2use(underGPNMB & ~remVals)==2  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==3  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==4  )/ sum((underGPNMB & ~remVals) )*100];



b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('unstaged','stage 1','stage 2','stage 3','stage 4', 'FontSize', 15)

ylabel('percent of total cases', 'FontSize', 15)

title(['Original BrettScheider  TDP43 Spread Scores in',  disName, ' Patients', ' p>.05' ], 'FontSize', 20)

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];



[p, n] = ranksum(path2use(  (overGPNMB' & ~remVals' &  path2use>0 )), path2use( underGPNMB' & ~remVals' & path2use>0  )   )  % Mann-Whitney U



%  [n,p]=prop_test([sum(path2use(overGPNMB & ~remVals )==4 ),sum(path2use(underGPNMB & ~remVals )==4   ) ] ...
%     ,[sum(overGPNMB & ~remVals), sum(underGPNMB& ~remVals)  ],true)
% 
% 
% 
% figure
% 
% path2use=cort4PathPos;
% 
% 
% b=bar([1,3,5], [nansum(path2use(overGPNMB & ~remVals )  )/ sum(overGPNMB(~remVals)), ...
%    nansum(path2use(het & ~remVals )  )/ sum(het(~remVals)),...
%    nansum(path2use(underGPNMB & ~remVals )  )/ sum(underGPNMB(~remVals))])
% 
% 
% 














%% pt 2 My  purely cortical pathology 

cortical1={'MC', 'Neocortical'};  cortical2={'SN','GP','CP','Cing','MB'}; cortical3={'CS','Amyg','DG','EC'};

pathPosBlank=nan(1,height(pathTable));

pathCut=1;
[cort1Path]=TDPStager(cortical1,pathTable,pathCut,pt2use);
cort1PathPos=pathPosBlank;cort1PathPos(cort1Path>1)=true;

[cort2Path]=TDPStager(cortical2,pathTable,pathCut,pt2use);
cort2PathPos=pathPosBlank;cort2PathPos(cort1Path>1 & cort2Path>1)=true;


[cort3Path]=TDPStager(cortical3,pathTable,pathCut,pt2use);
cort3PathPos=pathPosBlank;cort3PathPos(cort1Path>1 & cort2Path>1 & cort3Path>1)=true;


% [cort4Path]=TDPStager(cortical4,pathTable,pathCut,pt2use);
% cort4PathPos=pathPosBlank;cort4PathPos(cort1Path>1 & cort2Path>1 & cort3Path>1 & cort4Path>1)=true;


%now ranking path burden in patients with sufficent samples and outputting
%a simple stack bar graph 


path2use=nansum([cort1PathPos;cort2PathPos;cort3PathPos]);







remVals=~(pt2use) | emptySNP; %this includes pt2use ie non patients will be nanned out


figure



overGPNMBRat=[    sum(path2use(overGPNMB & ~remVals)==0  )/ sum((overGPNMB & ~remVals)  )*100,...
    sum(path2use(overGPNMB & ~remVals)==1  )/ sum((overGPNMB & ~remVals)  )*100,...
sum(path2use(overGPNMB & ~remVals)==2  )/ sum((overGPNMB & ~remVals) )*100,...
sum(path2use(overGPNMB & ~remVals)==3  )/ sum((overGPNMB & ~remVals) )*100];
% sum(path2use(overGPNMB & ~remVals)==4  )/ sum((overGPNMB & ~remVals) )*100];


    



hetGPNMBRat=[    sum(path2use(het & ~remVals)==0  )/ sum((het & ~remVals)  )*100,...
    sum(path2use(het & ~remVals)==1  )/ sum((het & ~remVals)  )*100,...
sum(path2use(het & ~remVals)==2  )/ sum((het & ~remVals) )*100,...
sum(path2use(het & ~remVals)==3  )/ sum((het & ~remVals) )*100];







underGPNMBRat=[    sum(path2use(underGPNMB & ~remVals)==0  )/ sum((underGPNMB & ~remVals)  )*100,...
    sum(path2use(underGPNMB & ~remVals)==1  )/ sum((underGPNMB & ~remVals)  )*100,...
sum(path2use(underGPNMB & ~remVals)==2  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==3  )/ sum((underGPNMB & ~remVals) )*100];



b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('unstaged','stage 1','stage 2','stage 3', 'FontSize', 15)

ylabel('percent of total cases', 'FontSize', 15)

[p, ~] = ranksum(path2use((overGPNMB' ) & ~remVals' & path2use>0), path2use( underGPNMB' & ~remVals'  & path2use>0)   );  % Mann-Whitney U


title(['Cortical Based TDP43 Spread Scores in',  disName, ' Patients', ' p=', num2str(p) ], 'FontSize', 20)

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];



 % [n,p]=prop_test([sum(path2use(overGPNMB & ~remVals )==1 ),sum(path2use(underGPNMB & ~remVals )==1   ) ] ...
 %    ,[sum(overGPNMB & ~remVals), sum(underGPNMB& ~remVals)  ],true)

% figure
% 
% % path2use=cort4PathPos;
% 
% 
% b=bar([1,3,5], [nansum(path2use(overGPNMB & ~remVals )  )/ sum(overGPNMB(~remVals)), ...
%    nansum(path2use(het & ~remVals )  )/ sum(het(~remVals)),...
%    nansum(path2use(underGPNMB & ~remVals )  )/ sum(underGPNMB(~remVals))])

%% part3 My pure spine based spread model

spine1={'SC'}; spine2={'Med'}; spine3={'Pons'};

pathCut=1;

pathPosBlank=nan(1,height(pathTable));

[spine1Path]=TDPStager(spine1,pathTable,pathCut,pt2use);
spine1PathPos= pathPosBlank; spine1PathPos(spine1Path>1)=true;



[spine2Path]=TDPStager(spine2,pathTable,pathCut,pt2use);
spine2PathPos= pathPosBlank; spine2PathPos(spine1Path>1 & spine2Path>1 )=true;


[spine3Path]=TDPStager(spine3,pathTable,pathCut,pt2use);
spine3PathPos= pathPosBlank; spine3PathPos(spine1Path>1 & spine2Path>1 & spine3Path>1 )=true;



path2use=nansum([spine1PathPos;spine2PathPos;spine3PathPos]);

remVals=~(pt2use) | emptySNP; %this includes pt2use ie non patients will be nanned out


figure



overGPNMBRat=[sum(path2use(overGPNMB & ~remVals)==0  )/ sum((overGPNMB & ~remVals)  )*100,...
    sum(path2use(overGPNMB & ~remVals)==1  )/ sum((overGPNMB & ~remVals)  )*100,...
sum(path2use(overGPNMB & ~remVals)==2  )/ sum((overGPNMB & ~remVals) )*100,...
sum(path2use(overGPNMB & ~remVals)==3  )/ sum((overGPNMB & ~remVals) )*100];

    


hetGPNMBRat=[sum(path2use(het & ~remVals)==0  )/ sum((het & ~remVals)  )*100,...
    sum(path2use(het & ~remVals)==1  )/ sum((het & ~remVals)  )*100,...
sum(path2use(het & ~remVals)==2  )/ sum((het & ~remVals) )*100,...
sum(path2use(het & ~remVals)==3  )/ sum((het & ~remVals)  )*100];




underGPNMBRat=[sum(path2use(underGPNMB & ~remVals)==0  )/ sum((underGPNMB & ~remVals) )*100,...
    sum(path2use(underGPNMB & ~remVals)==1  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==2  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==3  )/ sum((underGPNMB & ~remVals)  )*100];


b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('unstaged','stage 1','stage 2','stage 3', 'FontSize', 15)

ylabel('percent of total cases', 'FontSize', 15)

title(['Spine Based TDP43 Spread Scores in',  disName, ' Patients', ' p=' ], 'FontSize', 20)

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];




figure
%nansum([spine1PathPos;spine2PathPos;spine3PathPos]);

remVals=~(pt2use) | emptySNP; %this includes pt2use ie non patients will be nanned out

b=bar([1,3,5], [nansum(path2use(overGPNMB & ~remVals )  )/ sum(overGPNMB(~remVals)), ...
   nansum(path2use(het & ~remVals )  )/ sum(het(~remVals)),...
   nansum(path2use(underGPNMB & ~remVals )  )/ sum(underGPNMB(~remVals))])


 % [n,p]=prop_test([sum(path2use(overGPNMB & ~remVals )>=3 ),sum(path2use(underGPNMB & ~remVals )>=3  ) ] ...
 %    ,[sum(overGPNMB & ~remVals), sum(underGPNMB& ~remVals)  ],true)


[p, ~] = ranksum(path2use((overGPNMB' ) & ~remVals' & path2use>0), path2use( underGPNMB' & ~remVals'  & path2use>0)   )  % Mann-Whitney U







% 
% 
% 
% ageAtDeath = pathTable.AgeatDeath;
% ageAtOnset = pathTable.GlobalAgeOnset; 
% disDur = ageAtDeath - ageAtOnset;
% bioSex = pathTable.Sex;
% ptID=pathTable.INDDID;
% 
% 
% % Create a table with all collected data
% modelData = table(ptID(~remVals), ageAtDeath(~remVals), bioSex(~remVals), disDur(~remVals), snpStat(~remVals), path2use(~remVals), ...
%     'VariableNames', {'ptID', 'ageAtDeath', 'Sex', 'diseaseDuration', 'SNP', 'pathScore'});
% 
% % Convert categorical variables
% modelData.Sex = categorical(modelData.Sex);
% modelData.ptID = categorical(string(modelData.ptID));
% modelData.SNP = categorical(modelData.SNP);
% 
% % Fit GLME with BrainRegion as a cofactor
% glme = fitglme(modelData, ...
%     'pathScore ~ 1 + Sex + ageAtDeath + diseaseDuration + SNP  + (1|ptID)', ...
%     'Distribution', 'Normal', 'Link', 'Identity');
% 
% disp(glme);
% 
% title(['Spine Based TDP43 Spread Scores in',  disName, ' Patients' ], 'FontSize', 20)
% 
% 
% glme.Coefficients
% 
% % b=bar([1,3,5], [sum(path2use(overGPNMB & ~remVals )  )/ sum(overGPNMB(~remVals)), ...
% %    sum(path2use(het & ~remVals )  )/ sum(het(~remVals)),...
% %    sum(path2use(underGPNMB & ~remVals )  )/ sum(underGPNMB(~remVals))])
% 
% 


% groupLabels = ones(size(path2use)); % Placeholder for group labels
% groupLabels(overGPNMB & ~remVals) = 1;  % Over GPNMB
% groupLabels(underGPNMB & ~remVals) = 2; % Under GPNMB
% groupLabels(het & ~remVals) = 3;        % Heterozygous carriers
% 
% p = kruskalwallis(path2use(~remVals), groupLabels(~remVals));




    
% 
% hetGPNMBRat=[sum(val2Test(~remVals & het & pt2use )==4  )/ sum(((~remVals & het & pt2use   )))*100,...
% sum(val2Test(~remVals & het & pt2use )==3  )/ sum(((~remVals & het & pt2use )))*100,...
% sum(val2Test(~remVals & het & pt2use )==2  )/ sum(((~remVals & het & pt2use )))*100,...
% sum(val2Test(~remVals & het & pt2use )==1  )/ sum(((~remVals & het & pt2use )))*100,...
% sum(val2Test(~remVals & het & pt2use )==0  )/ sum(((~remVals & het & pt2use )))*100]; 
% 
% 
% 
% b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;
% 
% 
% 
% legend('3+','2+','1+','rare','none', 'FontSize', 15)
% 
% ylabel('percent of total cases', 'FontSize', 15)
% 
% title([brainAreaAtPlay,' ', pathName, ' Burden in ', disName, ' Patients' ], 'FontSize', 20)
% 
% ylim([0,120])
% 
% a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};
% 
% a.YTick(a.YTick>100)=[];






end