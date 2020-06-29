function [silhouette] = silhouette_score(X,label)

N = max(size(X));
subjects = 1:N;
C0 = subjects(label==0);
C1 = subjects(label==1);

D = dist(X');
a = zeros(N,1);
b = zeros(N,1);
s = zeros(N,1);
for i=1:N
    if ismember(i,C0)
        C = C0;
        nonC = C1;
    else
        C = C1;
        nonC = C0;
    end
    a(i,1) = 1/(length(C)-1)*sum(D(i,C)-D(i,i));
    b(i,1) = 1/length(nonC)*sum(D(i,nonC));
    s(i,1) = (b(i,1)-a(i,1))/max([a(i,1),b(i,1)]);
end

%silhouette = max(s);
silhouette = mean(s);

end