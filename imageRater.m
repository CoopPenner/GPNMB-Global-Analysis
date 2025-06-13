ratings=readtable('/Users/pennerc/Desktop/ImageRatingsACPTS.xlsx')


% Assume `ratings` is your table
labels = ratings.Genotype;
data = ratings{:, 2:end}; % extract numeric part (Var2–Var7)
avgData=[];
dataPoint=[];
LabelsAvg={'KO PBS', 'WT PBS', '2uM aBeta WT', '2uM aBeta KO', ' 3uM aBeta WT', '3uM aBeta KO'};

countr=1;
for dd=1:2:11

avgData(countr)=nanmean( [data(dd,:), data(dd+1,:)  ]);

dataPoint=[dataPoint; nanmean([data(dd,:); data(dd+1,:)]   ) ];

countr=countr+1;
end

% Boxplot with grouping
figure;
boxplot(data', {labels}, 'factorseparator', 1);
xlabel('Group');
ylabel('Rating');
title('Boxplot of Ratings by Condition and Variable');
xtickangle(45);





figure

bar(avgData)
a=gca;
a.XTickLabel=LabelsAvg;
hold on

    for dd=1:size(dataPoint,1)
    
    randVec=  (rand(1,5)/2.5) + dd;

    scatter(randVec, dataPoint(dd,:))
    
    end

    ylim([0,4])


[n,p]=ttest(  [mean(data(5,1:5);data(6,1:5)  ]  ,[mean(data(7,1:5);data(8,1:5)  ])
