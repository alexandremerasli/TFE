distance_matrix = dlmread('dist_EEG_nMDS.csv');
a = mink(distance_matrix,2);
distance_matrix = distance_matrix-repmat(min(a(2,:)),26,26);
distance_matrix = distance_matrix - diag(diag(distance_matrix));
%weights = dlmread('ISC steepness.csv');
%[a,idx] = sort(weights(:,1));
%weights = weights(idx,2);
weights = dlmread('weights_EEG.csv');
N = 26;
narrative = [ 0  2  4  5  7  9 11 13 16 17 18 22 25]+1;
stimulus = [ 1  3  6  8 10 12 14 15 19 20 21 23 24]+1;
subjects = 0:N-1;
condition = [0 1 0 1 0 0 1 0 1 0 1 0 1 0 1 1 0 0 0 1 1 1 0 1 1 0]';
K=2;
plotFigure = false;
N_init = 100;
%% nMDS
[points,stressN,disparitiesN] = mdscale(distance_matrix,2,"Start","cmdscale","Criterion","metricstress");
%[points,stressN,disparitiesN] = mdscale(distance_matrix,2,"Start","cmdscale","Criterion","strain","Weights",weights);
%writematrix(disparitiesN,'disparities.csv');
[points,stressN,disparitiesN] = mdscale(distance_matrix,2,"Start","cmdscale");

    figure();
    scatter(points(narrative,1),points(narrative,2),[],'blue','filled');
    hold on;
    scatter(points(stimulus,1),points(stimulus,2),[],'red','filled');

    for i=1:N
        text(points(i,1)+max(points,[],'all')/20,points(i,2),int2str(i-1));
    end

[~,ord] = sortrows([disparitiesN(:) distance_matrix(:)]);
computedDistancesN = pdist(points);
computedDistancesN = squareform(computedDistancesN);

    figure();
    plot(distance_matrix,computedDistancesN,'bo', ...
    distance_matrix(ord),disparitiesN(ord),'r.-');
    xlabel('Dissimilarities'); ylabel('Distances/Disparities')
    legend({'Distances' 'Disparities'},'Location','NW');
    corrN = corrcoef(computedDistancesN(:),disparitiesN(:));
    corrN = corrN(2,1);

label = kmeans(points,K,"Replicates",1000)-1;
disp(label');
NA = subjects(label==1);
SSA = subjects(label==2);
accuracy_nMDS = max(1-1/N*sum(label==condition),1/N*sum(label==condition));
misClassified = subjects(condition~=label)