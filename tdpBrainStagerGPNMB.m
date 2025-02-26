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




%I will talk to neuropathologists here about my scheme... but I am going to
%have a spinal and cerebral progression scale
%Spinal: 1) SC 2) Medulla 3) Pons (simple)
%cerebral: 1) MC 2) Neocortical/Cingulate 3) SN GP CP and MB 4) DG Amyg
%Hippo EC



cortical1={'MC'}; cortical2={'Cing','Neocortical'}; cortical3={'SN','GP','CP'}; cortical4={'CS','Amyg','DG','EC'};


pathCut=1;
[cort1Path,cort1PathPos]=TDPStager(cortical1,pathTable,pathCut,pt2use);

[cort2Path,cort2PathPos]=TDPStager(cortical2,pathTable,pathCut,pt2use);
cort2PathPos(cort1PathPos==0)=0;

[cort3Path,cort3PathPos]=TDPStager(cortical3,pathTable,pathCut,pt2use);
cort3PathPos(cort2PathPos==0)=0;


[cort4Path,cort4PathPos]=TDPStager(cortical4,pathTable,pathCut,pt2use);
cort4pathpos(cort3PathPos==0)=0;


%now ranking path burden in patients with sufficent samples and outputting
%a simple stack bar graph 




figure

path2use=cort4PathPos;

remVals=isnan(path2use) | emptySNP; %this includes pt2use ie non patients will be nanned out

b=bar([1,3,5], [sum(path2use(overGPNMB & ~remVals )  )/ sum(overGPNMB(~remVals)), ...
   sum(path2use(het & ~remVals )  )/ sum(het(~remVals)),...
   sum(path2use(underGPNMB & ~remVals )  )/ sum(underGPNMB(~remVals))])



spine1={'SC'}; spine2={'Med'}; spine3={'Pons'};

pathCut=1

[spine1Path,spine1PathPos]=TDPStager(spine1,pathTable,pathCut,pt2use);

[spine2Path,spine2PathPos]=TDPStager(spine2,pathTable,pathCut,pt2use);
spine2PathPos(spine1PathPos==0)=0;


[spine3Path,spine3PathPos]=TDPStager(spine3,pathTable,pathCut,pt2use);
spine3PathPos(spine2PathPos==0)=0;



path2use=spine1PathPos+spine2PathPos+spine3PathPos;

remVals=isnan(path2use) | emptySNP |path2use==0;

figure
b=bar([1,3,5], [nanmean(path2use(overGPNMB & ~remVals )  ), ...
   nanmean(path2use(het & ~remVals )  ),...
   nanmean(path2use(underGPNMB & ~remVals )  )])

underVal=path2use(underGPNMB & ~remVals ); overVal=path2use(overGPNMB & ~remVals );


figure



overGPNMBRat=[sum(path2use(overGPNMB & ~remVals)==1  )/ sum((overGPNMB & ~remVals)  )*100,...
sum(path2use(overGPNMB & ~remVals)==2  )/ sum((overGPNMB & ~remVals) )*100,...
sum(path2use(overGPNMB & ~remVals)==3  )/ sum((overGPNMB & ~remVals) )*100]

    


hetGPNMBRat=[sum(path2use(het & ~remVals)==1  )/ sum((het & ~remVals)  )*100,...
sum(path2use(het & ~remVals)==2  )/ sum((het & ~remVals) )*100,...
sum(path2use(het & ~remVals)==3  )/ sum((het & ~remVals)  )*100];




underGPNMBRat=[sum(path2use(underGPNMB & ~remVals)==1  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==2  )/ sum((underGPNMB & ~remVals) )*100,...
sum(path2use(underGPNMB & ~remVals)==3  )/ sum((underGPNMB & ~remVals)  )*100]


b=bar([1,3,5], [overGPNMBRat;hetGPNMBRat;underGPNMBRat], 'stacked') ;



legend('stage 1','stage 2','stage 3', 'FontSize', 15)

ylabel('percent of total cases', 'FontSize', 15)

title(['Spine Based TDP43 Spread Scores in',  disName, ' Patients', ' p=' ], 'FontSize', 20)

ylim([0,120])

a=gca; a.XTickLabel={'AA overProduction', 'AG', 'GG underProduction'};

a.YTick(a.YTick>100)=[];











ageAtDeath = pathTable.AgeatDeath;
ageAtOnset = pathTable.GlobalAgeOnset; 
disDur = ageAtDeath - ageAtOnset;
bioSex = pathTable.Sex;
ptID=pathTable.INDDID;


% Create a table with all collected data
modelData = table(ptID(~remVals), ageAtDeath(~remVals), bioSex(~remVals), disDur(~remVals), snpStat(~remVals), path2use(~remVals), ...
    'VariableNames', {'ptID', 'ageAtDeath', 'Sex', 'diseaseDuration', 'SNP', 'pathScore'});

% Convert categorical variables
modelData.Sex = categorical(modelData.Sex);
modelData.ptID = categorical(string(modelData.ptID));
modelData.SNP = categorical(modelData.SNP);

% Fit GLME with BrainRegion as a cofactor
glme = fitglme(modelData, ...
    'pathScore ~ 1 + Sex + ageAtDeath + diseaseDuration + SNP  + (1|ptID)', ...
    'Distribution', 'Normal', 'Link', 'Identity');

disp(glme);

title(['Spine Based TDP43 Spread Scores in',  disName, ' Patients' ], 'FontSize', 20)


glme.Coefficients

% b=bar([1,3,5], [sum(path2use(overGPNMB & ~remVals )  )/ sum(overGPNMB(~remVals)), ...
%    sum(path2use(het & ~remVals )  )/ sum(het(~remVals)),...
%    sum(path2use(underGPNMB & ~remVals )  )/ sum(underGPNMB(~remVals))])


 [n,p]=prop_test([sum(path2use(overGPNMB & ~remVals )==2 ),sum(path2use(underGPNMB & ~remVals )==2   ) ] ...
    ,[sum(overGPNMB & ~remVals), sum(underGPNMB& ~remVals)  ],true)


[p, ~] = ranksum(path2use(  (overGPNMB & ~remVals)), path2use( (underGPNMB) & ~remVals));  % Mann-Whitney U



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