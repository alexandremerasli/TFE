distance_matrix = dlmread('dist_EEG.csv');
%distance_matrix = dlmread('dist_EDA.csv');
%distance_matrix = dlmread('dist_IBI.csv');
N = 26;
narrative = [ 0  2  4  5  7  9 11 13 16 17 18 22 25]+1;
stimulus = [ 1  3  6  8 10 12 14 15 19 20 21 23 24]+1;
subjects = 1:N;
condition = [0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 1 0 0 0 1 1 1 0 1 1 0]';
K=2;
plotFigure = false;
N_init = 10;
algo = 'K-Means';
algo = 'Spectral Clustering';
algo = 'Hierarchical Clustering';
algo = 'K-Medoids';
%% cMDS
close all;
N_init = 1;
labelmMDS = zeros(N_init,N);
corr_matrix = zeros(N,N);
plotFigure = true;

for nb_init=1:N_init
    disp(nb_init)
    [points,strainC] = cmdscale(distance_matrix,2);
    %[points,strainC] = cmdscale(disparitiesN,2);
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
    
    % Apply clustering
    if strcmp(algo,'K-Means')
        label = kmeans(points,K,'Replicates',100)-1;
    elseif strcmp(algo,'Spectral Clustering')
        label = spectralclusterTest(points,K)-1;
    elseif strcmp(algo,'K-Medoids')
        label = kmedoids(points,K,'Replicates',100)-1;
    elseif strcmp(algo,'Hierarchical Clustering')
        Z = linkage(points,'ward');
        label = cluster(Z,'Maxclust',2)-1;
    end
    
    % SC precomputed
    %affinity_matrix = exp(-1 * distance_matrix.*distance_matrix);
    %label = spectralcluster(affinity_matrix,K,'Distance','precomputed')-1;
    
    NA = subjects(label==1);
    SSA = subjects(label==2);
    accuracy_cMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
    labelcMDS(nb_init,:) = label;
    %disp(nb_init)
    corr_matrix = corr_matrix + (label==label')/N_init;
end

for i=1:N_init
    tmpAcc = sum(labelcMDS(1,:)==labelcMDS(i,:))/N;
    if tmpAcc<0.5
        labelcMDS(i,:) = 1-labelcMDS(i,:);
    end
end

probaLabel = 1/N_init*sum(labelcMDS,1);
label = (probaLabel>0.5)';
accuracy_cMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
%silhouette_cMDS = mean(silhouette(points,label))
silhouettecMDS = silhouette_score(points,label)
dbIndex_cMDS = evalclusters(points,'kmeans','DaviesBouldin','klist',[2:2]).CriterionValues

disp(probaLabel)
wellClassified = subjects(condition==label)-1
misClassified = subjects(condition~=label)-1

%label = kmeans(points,K)-1;
%NA = subjects(label==1);
%SSA = subjects(label==2); 
%accuracy_cMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition)); 
%% mMDS
close all;
N_init = 300;
labelmMDS = zeros(N_init,N);
corr_matrix = zeros(N,N);
% plotFigure = true;
plotFigure = false;

for nb_init=1:N_init
    %s=rng(2);
    %init = rand(N,2);
    [points,stressM,disparitiesM] = mdscale(distance_matrix,2,"Criterion","metricstress","Start","random");
    %[points,stressM,disparitiesM] = mdscale(distance_matrix,2,"Criterion","metricstress");
    
    if plotFigure
        figure();
        scatter(points(narrative,1),points(narrative,2),[],'blue','filled');
        hold on;
        scatter(points(stimulus,1),points(stimulus,2),[],'red','filled');
        
        for i=1:N
            text(points(i,1)+max(points,[],'all')/20,points(i,2),int2str(i-1));
        end
    end
    [~,ord] = sortrows([disparitiesM(:) distance_matrix(:)]);
    computedDistancesM = pdist(points);
    computedDistancesM = squareform(computedDistancesM);
    
    if plotFigure
        figure();
        plot(distance_matrix,computedDistancesM,'bo', ...
        distance_matrix(ord),disparitiesM(ord),'r.-');
        xlabel('Dissimilarities'); ylabel('Distances');
        legend({'Distances' 'Disparities'},'Location','NW');
    end
    corrM = corrcoef(computedDistancesM(:),disparitiesM(:));
    corrM = corrM(2,1);
    
    % Apply clustering
    if strcmp(algo,'K-Means')
        label = kmeans(points,K,'Replicates',100)-1;
    elseif strcmp(algo,'Spectral Clustering')
        label = spectralclusterTest(points,K)-1;
    elseif strcmp(algo,'K-Medoids')
        label = kmedoids(points,K,'Replicates',100)-1;
    elseif strcmp(algo,'Hierarchical Clustering')
        Z = linkage(points,'ward');
        label = cluster(Z,'Maxclust',2)-1;
    end
    
    % SC precomputed
    %affinity_matrix = exp(-1 * distance_matrix.*distance_matrix);
    %label = spectralcluster(affinity_matrix,K,'Distance','precomputed')-1;
    
    NA = subjects(label==1);
    SSA = subjects(label==2);
    accuracy_mMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
    labelmMDS(nb_init,:) = label;
    %disp(nb_init)
    corr_matrix = corr_matrix + (label==label')/N_init;
end

for i=1:N_init
    tmpAcc = sum(labelmMDS(1,:)==labelmMDS(i,:))/N;
    if tmpAcc<0.5
        labelmMDS(i,:) = 1-labelmMDS(i,:);
    end
end

probaLabel = 1/N_init*sum(labelmMDS,1);
label = (probaLabel>0.5)';
accuracy_mMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
%silhouette_mMDS = mean(silhouette(points,label))
silhouettemMDS = silhouette_score(points,label)
disp(probaLabel)
wellClassifiedmMDS = subjects(condition==label)-1
misClassifiedmMDS = subjects(condition~=label)-1
%% nMDS
N_init = 300;
labelmMDS = zeros(N_init,N);
corr_matrix = zeros(N,N);

for nb_init=1:N_init
    %s=rng(2);
    %init = rand(N,2);
    [points,stressM,disparitiesN] = mdscale(distance_matrix,2,"Criterion","stress","Start","random");
    %[points,stressM,disparitiesN] = mdscale(distance_matrix,2,"Criterion","stress");
    
    if plotFigure
        figure();
        scatter(points(narrative,1),points(narrative,2),[],'blue','filled');
        hold on;
        scatter(points(stimulus,1),points(stimulus,2),[],'red','filled');
        
        for i=1:N
            text(points(i,1)+max(points,[],'all')/20,points(i,2),int2str(i-1));
        end
    end
    [~,ord] = sortrows([disparitiesN(:) distance_matrix(:)]);
    computedDistancesN = pdist(points);
    computedDistancesN = squareform(computedDistancesN);
    
    if plotFigure
        figure();
        plot(distance_matrix,computedDistancesN,'bo', ...
        distance_matrix(ord),disparitiesN(ord),'r.-');
        xlabel('Dissimilarities'); ylabel('Distances');
        legend({'Distances' 'Disparities'},'Location','NW');
    end
    corrN = corrcoef(computedDistancesN(:),disparitiesN(:));
    corrN = corrN(2,1);
    
    % Apply clustering
    if strcmp(algo,'K-Means')
        label = kmeans(points,K,'Replicates',100)-1;
    elseif strcmp(algo,'Spectral Clustering')
        label = spectralclusterTest(points,K)-1;
    elseif strcmp(algo,'K-Medoids')
        label = kmedoids(points,K,'Replicates',100)-1;
    elseif strcmp(algo,'Hierarchical Clustering')
        Z = linkage(points,'ward');
        label = cluster(Z,'Maxclust',2)-1;
    end
    
    % SC precomputed
    %affinity_matrix = exp(-1 * distance_matrix.*distance_matrix);
    %label = spectralcluster(affinity_matrix,K,'Distance','precomputed')-1;
    SSA = subjects(label==2);
    accuracy_nMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
    labelnMDS(nb_init,:) = label;
    
    corr_matrix = corr_matrix + (label==label')/N_init;
end

for i=1:N_init
    tmpAcc = sum(labelnMDS(1,:)==labelnMDS(i,:))/N;
    if tmpAcc<0.5
        labelnMDS(i,:) = 1-labelnMDS(i,:);
    end
end

probaLabel = 1/N_init*sum(labelnMDS,1);
label = (probaLabel>0.5)';
accuracy_nMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
%silhouette_nMDS = mean(silhouette(points,label))
silhouettenMDS = silhouette_score(points,label)
disp(probaLabel)
wellClassifiednMDS = subjects(condition==label)-1
misClassifiednMDS = subjects(condition~=label)-1
%label = kmeans(points,K)-1;
%NA = subjects(label==1);
%SSA = subjects(label==2);
%accuracy_nMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition)); 

%% Run Python
% system('python MDS_normalized.py');
