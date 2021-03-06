# -*- coding: utf-8 -*-

# https://github.com/Behrouz-Babaki/COP-Kmeans

import random
import numpy as np
from sklearn.utils import check_random_state
from sklearn.metrics import pairwise_distances

class COPKmeans():    
    
    def __init__(self, n_clusters,must_link=[],cannot_link=[],
               initialization='kmpp',
               max_iter=300, tol=1e-4,n_init=10):
        self.n_clusters = n_clusters
        self.ml = must_link
        self.cl = cannot_link
        self.init = initialization
        self.max_iter = max_iter
        self.tol = tol
        self.labels_ = None
        self.n_init = n_init
        
    def fit(self,dataset):
        n_init = self.n_init
        random_state = check_random_state(None)

        seeds = random_state.randint(np.iinfo(np.int32).max,size=n_init)
        self.inertia_ = None
        
        for seed in seeds:
            # run a COP k-means once
            labels, inertia, centers = self.copkmeans_single(dataset)
            # determine if these results are the best so far
            if self.inertia_ is None or inertia < self.inertia_:
                self.labels_ = labels.copy()
                self.centers_ = centers.copy()
                self.inertia_ = inertia.copy()
                
        return self
    
    def compute_inertia(self,dataset,labels):
        subjects = np.arange(0,len(labels))
        foundNA = subjects[labels==0]
        cluster0 = dataset[labels==0]
        center0 = np.mean(cluster0,axis=0)
        foundSSA = subjects[labels==1]
        cluster1 = dataset[labels==1]
        center1 = np.mean(cluster1,axis=0)

        distArray0 = pairwise_distances(cluster0,np.tile(center0,(cluster0.shape[0],1)))[:,0]
        distArray1 = pairwise_distances(cluster1,np.tile(center1,(cluster1.shape[0],1)))[:,0]
    
        return np.sum(distArray0) + np.sum(distArray1)
        
    def copkmeans_single(self,dataset):
        ml = self.ml
        cl = self.cl
        initialization = self.init
        max_iter = self.max_iter
        tol = self.tol
        k = self.n_clusters
        
        ml, cl = self.transitive_closure(ml, cl, len(dataset))
        ml_info = self.get_ml_info(ml, dataset)
        tol = self.tolerance(tol, dataset)

        centers = self.initialize_centers(dataset, k, initialization)

        for _ in range(max_iter):
            clusters_ = [-1] * len(dataset)
            for i, d in enumerate(dataset):
                indices, _ = self.closest_clusters(centers, d)
                counter = 0
                if clusters_[i] == -1:
                    found_cluster = False
                    while (not found_cluster) and counter < len(indices):
                        index = indices[counter]
                        if not self.violate_constraints(i, index, clusters_, ml, cl):
                            found_cluster = True
                            clusters_[i] = index
                            for j in ml[i]:
                                clusters_[j] = index
                        counter += 1

                    if not found_cluster:
                        return None, None

            clusters_, centers_ = self.compute_centers(clusters_, dataset, k, ml_info)
            shift = sum(self.l2_distance(centers[i], centers_[i]) for i in range(k))
            if shift <= tol:
                break

            centers = centers_

        labels = np.array(clusters_)
        inertia = self.compute_inertia(dataset,labels)
        centers = np.array(centers_)
        return labels, inertia, centers

    def l2_distance(self,point1, point2):
        return sum([(float(i)-float(j))**2 for (i, j) in zip(point1, point2)])

    # taken from scikit-learn (https://goo.gl/1RYPP5)
    def tolerance(self,tol, dataset):
        n = len(dataset)
        dim = len(dataset[0])
        averages = [sum(dataset[i][d] for i in range(n))/float(n) for d in range(dim)]
        variances = [sum((dataset[i][d]-averages[d])**2 for i in range(n))/float(n) for d in range(dim)]
        return tol * sum(variances) / dim

    def closest_clusters(self,centers, datapoint):
        distances = [self.l2_distance(center, datapoint) for
                     center in centers]
        return sorted(range(len(distances)), key=lambda x: distances[x]), distances

    def initialize_centers(self,dataset, k, method):
        if method == 'random':
            ids = list(range(len(dataset)))
            random.shuffle(ids)
            return [dataset[i] for i in ids[:k]]

        elif method == 'kmpp':
            chances = [1] * len(dataset)
            centers = []

            for _ in range(k):
                chances = [x/sum(chances) for x in chances]
                r = random.random()
                acc = 0.0
                for index, chance in enumerate(chances):
                    if acc + chance >= r:
                        break
                    acc += chance
                centers.append(dataset[index])

                for index, point in enumerate(dataset):
                    cids, distances = self.closest_clusters(centers, point)
                    chances[index] = distances[cids[0]]

            return centers

    def violate_constraints(self,data_index, cluster_index, clusters, ml, cl):
        for i in ml[data_index]:
            if clusters[i] != -1 and clusters[i] != cluster_index:
                return True

        for i in cl[data_index]:
            if clusters[i] == cluster_index:
                return True

        return False

    def compute_centers(self,clusters, dataset, k, ml_info):
        cluster_ids = set(clusters)
        k_new = len(cluster_ids)
        id_map = dict(zip(cluster_ids, range(k_new)))
        clusters = [id_map[x] for x in clusters]

        dim = len(dataset[0])
        centers = [[0.0] * dim for i in range(k)]

        counts = [0] * k_new
        for j, c in enumerate(clusters):
            for i in range(dim):
                centers[c][i] += dataset[j][i]
            counts[c] += 1

        for j in range(k_new):
            for i in range(dim):
                centers[j][i] = centers[j][i]/float(counts[j])

        if k_new < k:
            ml_groups, ml_scores, ml_centroids = ml_info
            current_scores = [sum(self.l2_distance(centers[clusters[i]], dataset[i])
                                  for i in group)
                              for group in ml_groups]
            group_ids = sorted(range(len(ml_groups)),
                               key=lambda x: current_scores[x] - ml_scores[x],
                               reverse=True)

            for j in range(k-k_new):
                gid = group_ids[j]
                cid = k_new + j
                centers[cid] = ml_centroids[gid]
                for i in ml_groups[gid]:
                    clusters[i] = cid

        return clusters, centers

    def get_ml_info(self,ml, dataset):
        flags = [True] * len(dataset)
        groups = []
        for i in range(len(dataset)):
            if not flags[i]: continue
            group = list(ml[i] | {i})
            groups.append(group)
            for j in group:
                flags[j] = False

        dim = len(dataset[0])
        scores = [0.0] * len(groups)
        centroids = [[0.0] * dim for i in range(len(groups))]

        for j, group in enumerate(groups):
            for d in range(dim):
                for i in group:
                    centroids[j][d] += dataset[i][d]
                centroids[j][d] /= float(len(group))

        scores = [sum(self.l2_distance(centroids[j], dataset[i])
                      for i in groups[j])
                  for j in range(len(groups))]

        return groups, scores, centroids

    def transitive_closure(self,ml, cl, n):
        ml_graph = dict()
        cl_graph = dict()
        for i in range(n):
            ml_graph[i] = set()
            cl_graph[i] = set()

        def add_both(d, i, j):
            d[i].add(j)
            d[j].add(i)

        for (i, j) in ml:
            add_both(ml_graph, i, j)

        def dfs(i, graph, visited, component):
            visited[i] = True
            for j in graph[i]:
                if not visited[j]:
                    dfs(j, graph, visited, component)
            component.append(i)

        visited = [False] * n
        for i in range(n):
            if not visited[i]:
                component = []
                dfs(i, ml_graph, visited, component)
                for x1 in component:
                    for x2 in component:
                        if x1 != x2:
                            ml_graph[x1].add(x2)
        for (i, j) in cl:
            add_both(cl_graph, i, j)
            for y in ml_graph[j]:
                add_both(cl_graph, i, y)
            for x in ml_graph[i]:
                add_both(cl_graph, x, j)
                for y in ml_graph[j]:
                    add_both(cl_graph, x, y)

        for i in ml_graph:
            for j in ml_graph[i]:
                if j != i and j in cl_graph[i]:
                    raise Exception('inconsistent constraints between %d and %d' %(i, j))

        return ml_graph, cl_graph