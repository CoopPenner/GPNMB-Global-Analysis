a%ALS comprehensive GPNMB eval


clinDataTot=readtable('/Volumes/PC60/InqueryDatasets/ALS_Comprehensive_inquery_output.xlsx');


GPNMBSNP=clinDataTot.rs199347;
ID=clinDataTot.INDDID;
ageAtOnset=(clinDataTot.GlobalAgeOnset);
ageAtDeath=(clinDataTot.AgeatDeath);
disDur=ageAtDeath-ageAtOnset;
Sex=clinDataTot.Sex;
milService=clinDataTot.VeteranService;
tobaccoYears=clinDataTot.SmokingPackYears;
snpStat=clinDataTot.rs199347;
pathID=clinDataTot.INDDID;
SOD1=contains(clinDataTot.Mutation_Summary,'SOD1');
geneCarrier=~cellfun(@isempty,clinDataTot.Mutation_Summary) & ~SOD1;
trachDate=clinDataTot.TracheostomyDate;
hospDate=clinDataTot.HospiceEnrollmentDate;
pegDate=clinDataTot.PEGDate;
onsetDate=clinDataTot.ALSSymptomOnsetDate;
diagDate=clinDataTot.DiagnosisDate;
onsetSite=clinDataTot.ALSSymptomOnsetSite;



overGPNMB= contains(GPNMBSNP,'TT'); %the major allele
het=contains(GPNMBSNP, 'CT');
underGPNMB=contains(GPNMBSNP,'CC'); %the minor allele


UMN=contains(onsetSite, 'UMN');
LMN=contains(onsetSite,'LMN');
cogOnset=contains(onsetSite,'Cognitive');
bulbarOnset=contains(onsetSite,'Bulbar');
cervicalOnset=contains(onsetSite,'Cervical');









val2Plot=  disDur ;
remVals=disDur |cellfun(@isempty,GPNMBSNP) ;
%really clear and interesting trend for age at onset and disease duration



snpVar=categorical(GPNMBSNP(~remVals))';
ageVar=ageAtOnset(~remVals)';
sexVar=categorical(Sex(~remVals))';
geneCarrierVar= categorical(geneCarrier(~remVals));
SOD1Var=categorical(SOD1(~remVals));
onsetSiteVar=categorical(onsetSite(~remVals));


varNames=["Sex","cogSlope", "rs199347","ageAtTest","edLevel",'dx','ID','startScore','endScore'];

glmeTable = table(sexVar ,slopeVar,snpVar,ageVar,edVar,diagVar,IDVar, startVar, endVar, 'Variablenames',varNames);

glmeTable.rs199347 = reordercats(glmeTable.rs199347, {'TT', 'CC', 'CT'});



glme = fitglme(glmeTable,...
		'cogSlope ~ 1  + rs199347 + Sex + startScore + (1|ageAtTest)   ',...
		'Distribution','normal','Link','identity','FitMethod','Laplace',...
		'DummyVarCoding','reference')


glmeHeldOut=fitglme(glmeTable,...
		'cogSlope ~ 1 + Sex + startScore + (1|ageAtTest)   ',...
		'Distribution','normal','Link','identity','FitMethod','Laplace',...
		'DummyVarCoding','reference');


compare(glme,glmeHeldOut)






figure
subplot(1,2,1)
ALSSpecGPNMBPlotr(GPNMBSNP,val2Plot, true(1,length(disDur)), 'Disease Duration', 'ALS',remVals)


val2Plot=  ageAtOnset ;
remVals=disDur |cellfun(@isempty,GPNMBSNP) ;
subplot(1,2,2)
ALSSpecGPNMBPlotr(GPNMBSNP,val2Plot, true(1,length(ageAtOnset)), 'Age At Onset', 'ALS',remVals)












       pureTime=datenum(relDates); %convert it to pure time to sort
        pureTimeSample=datenum(sampleDate2use);
        [sorted,n]=sort(pureTime); %sorting on absolute value of time
        testYear=testYear(n); testMonth=testMonth(n); relScores=relScores(n); %then just outputting the reorder


        
        %now make every test date relative to the sample date

            testYear=testYear-sampleYear;
            testYear=testYear*12;
            testMonth=testMonth-sampleMonth;
            
            relativeTime=testYear+testMonth;

            preSamp=relativeTime<0;
            relativeTime(preSamp)=[]; relScores(preSamp)=[]; %remove any test dates that came before the sample was taken

        if ~isempty(relScores)

         [uniqueA, i,j] = unique(relativeTime,'first');
         idxToRemove = find(not(ismember(1:numel(relativeTime),i)))   ;
         relativeTime(idxToRemove)=[]; relScores(idxToRemove)=[]; %sometimes for some reason there are duplicate times
        relativeTime((isnan(relScores)))=[];
        relScores(isnan(relScores))=[]; 
        relativeTime=relativeTime/12; %everything contextualized in years, since that is the standard by which slopes are presented


        scoreVec=diff(relScores);  sigChange=find( (relScores-relScores(1) ) <=-2 ,1); %here we identify when scores sink two points below first measure.      
         P=polyfit(relativeTime,relScores,1);
        slopeCollect(tt)=P(1);
     
if isnan(P(1)) && condAtPlay==2
d=1
end

        if ~isempty(sigChange)      
    year2DecCollect(tt)= relativeTime(sigChange)  ;
        else
     year2DecCollect(tt)= nan ;
        end
    startScore(tt)=relScores(1); %collect starting point
    ageCollect(tt)= data2use.AGE__YRS(tt); %collect age at first check

    if strcmp(data2use.Sex(tt),'Male')
    sexCollect(tt)= 1; %coding male as 1
    elseif strcmp(data2use.Sex(tt),'Female')
    sexCollect(tt)= 2; %coding female as 2
    end

        end
            end

            if condAtPlay==1
        origSlope=slopeCollect;
        origID=countID;
        origDec=year2DecCollect;
        origDement=startScore;
        idxToRemoveOrig=idxToRemoveTot;
        origAge=ageCollect;
        origSex=sexCollect;
            elseif condAtPlay==2
        repSlope=slopeCollect;
        repID=countID;
        repDec=year2DecCollect;
        repDement=startScore;
        idxToRemoveRep=idxToRemoveTot;
        repAge=ageCollect;
        repSex=sexCollect;
            end
end















%really clear and interesting trend for age at onset and disease duration
figure
subplot(1,2,1)
ALSSpecGPNMBPlotr(GPNMBSNP,disDur, true(1,length(disDur)), 'Disease Duration', 'ALS',LMN, SOD1)
subplot(1,2,2)
ALSSpecGPNMBPlotr(GPNMBSNP,ageAtOnset, true(1,length(disDur)), 'Age At Onset', 'ALS',LMN,SOD1)

