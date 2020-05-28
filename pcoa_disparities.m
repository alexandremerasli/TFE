distance_matrix = dlmread('dist_EEG_mMDS.csv');
a = mink(distance_matrix,2);
distance_matrix = distance_matrix-repmat(min(a(2,:)),26,26);
distance_matrix = distance_matrix - diag(diag(distance_matrix));
%disparitiesN = dlmread('dist_EEG_mMDS.csv');
N = 26;
narrative = [ 0  2  4  5  7  9 11 13 16 17 18 22 25]+1;
stimulus = [ 1  3  6  8 10 12 14 15 19 20 21 23 24]+1;
subjects = 1:N;
condition = [0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 1 0 0 0 1 1 1 0 1 1 0]';
K=2;
plotFigure = true;
N_init = 100;
%% cMDS
%close all;
N_init = 1;
labelmMDS = zeros(N_init,N);
corr_matrix = zeros(N,N);
% plotFigure = true;

for nb_init=1:N_init
    disp(nb_init)
    [points,strainC] = cmdscale(distance_matrix,2);
    %[points,strainC] = cmdscale(disparitiesN,2);
    %[points,stressM,disparitiesM] = mdscale(disparitiesN,2,"Criterion","metricstress");
    disparitiesC = distance_matrix;
    if plotFigure
        figure();
        scatter(points(narrative,1),points(narrative,2),[],'blue','filled');
        hold on;
        scatter(points(stimulus,1),points(stimulus,2),[],'red','filled');
        
        for i=1:N
            text(points(i,1)+max(points,[],'all')/20,points(i,2),int2str(i-1));
        end
    end
    [~,ord] = sortrows([disparitiesC(:) distance_matrix(:)]);
    computedDistancesC = pdist(points);
    computedDistancesC = squareform(computedDistancesC);
    
    if plotFigure
        figure();
        plot(distance_matrix,computedDistancesC,'bo', ...
        distance_matrix(ord),disparitiesC(ord),'r.-');
        xlabel('Dissimilarities'); ylabel('Distances');
        legend({'Distances' 'Disparities'},'Location','NW');
    end
    corrC = corrcoef(computedDistancesC(:),disparitiesC(:));
    corrC = corrC(2,1);
    
    label = kmeans(points,K,'Replicates',20)-1;
    %label = spectralcluster(points,K)-1;
    NA = subjects(label==1);
    SSA = subjects(label==2);
    accuracy_cMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
    labelcMDS(nb_init,:) = label;
    %disp(nb_init)
    corr_matrix = corr_matrix + (label==label')/N_init;
end