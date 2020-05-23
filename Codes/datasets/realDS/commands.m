%%% Reading commands

%% Yeast data
yeast = readtable('yeast.dat');

Yeast = table2array(yeast(:,2:9));
figure; gscatter(Yeast(:,1),Yeast(:,2),ones(size(Yeast,1),1),'b'); grid on;
figure; gscatter(Yeast(:,1),Yeast(:,2)./100,ones(size(Yeast,1),1),'b'); grid on;

YeastMod = Yeast(:,[1:5 7:8]);

sum(strcmp(yeast{:,'Var10'},'CYT'))
strcmp(yeast{:,'Var10'},'CYT')
yeast_CYT = table2array(yeast(strcmp(yeast{:,'Var10'},'CYT'),1:8));

sum(strcmp(yeast{:,'Var10'},'NUC'))
strcmp(yeast{:,'Var10'},'NUC')
yeast_NUC = table2array(yeast(strcmp(yeast{:,'Var10'},'NUC'),1:8));

sum(strcmp(yeast{:,'Var10'},'MIT'))
strcmp(yeast{:,'Var10'},'MIT')
yeast_MIT = table2array(yeast(strcmp(yeast{:,'Var10'},'MIT'),1:8));


%% Mamography data
mamography = readtable('mammography.csv');

sum(strcmp(mamography{:,7},{'''-1'''}))
sum(strcmp(mamography{:,7},{'''1'''}))

mamography_1 = table2array(mamography(strcmp(mamography{:,7},'''-1'''),1:6));
mamography_2 = table2array(mamography(strcmp(mamography{:,7},'''1'''),1:6));


