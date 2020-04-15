import numpy as np

class KMedoids():    
    
    def __init__(self, n_clusters, tmax=100,printInfo=False,init=None):
        self.n_clusters = n_clusters
        self.tmax = tmax
        self.printInfo = printInfo
        self.init = init
        
    def fit(self,D):
        k = self.n_clusters
        tmax = self.tmax

        # determine dimensions of distance matrix D
        m, n = D.shape
        N = n

        # randomly initialize an array of k medoid indices
        M = self.safeInit(D,k,n)

        # create a copy of the array of medoid indices
        Mnew = np.copy(M)
        # initialize a dictionary to represent clusters
        C = {}
        for t in range(tmax):
            # determine clusters, i.e. arrays of data indices
            J = np.argmin(D[:,M], axis=1)
            for kappa in range(k):
                C[kappa] = np.where(J==kappa)[0]

            # update cluster medoids
            for kappa in range(k):
                J = np.mean(D[np.ix_(C[kappa],C[kappa])],axis=1)
                j = np.argmin(J)
                Mnew[kappa] = C[kappa][j]
            np.sort(Mnew)

            # check for convergence
            if np.array_equal(M, Mnew):
                if (self.printInfo):
                    print("procedure converged after %i iterations" %t)
                break
                M = np.copy(Mnew)

        if (t==tmax-1): # procedure did not converge
            if (self.printInfo):
                print("procedure did not converge")
            # final update of cluster memberships
            J = np.argmin(D[:,M], axis=1)
            for kappa in range(k):
                C[kappa] = np.where(J==kappa)[0]

        NA_medoids = C[0] # Maybe SSA
        SSA_medoids = C[1] # Maybe NA
        label = np.zeros(N)
        for subj in range(N):
            if subj in NA_medoids:
                label[subj] = 0
            else:
                label[subj] = 1
        self.labels_ = label
        return self

    def safeInit(self,D,k,n):
        if (self.init == None):
            M = np.sort(np.random.choice(n, k))
            C = {}
            J = np.argmin(D[:,M], axis=1)
            for kappa in range(k):
                C[kappa] = np.where(J==kappa)[0]

            # update cluster medoids
            for kappa in range(k):
                if (len(D[np.ix_(C[kappa],C[kappa])]) == 0):
                    if (self.printInfo):
                        print("Computing new init")
                    M = self.safeInit(D,k,n)
        else:
            # initialize k medoid indices thanks to HC results
            M = np.zeros((k,1),dtype=int)
            J = np.mean(D[np.ix_(self.init[0],self.init[0])],axis=1)
            j = np.argmin(J)
            M[0] = self.init[0][j]

            J = np.mean(D[np.ix_(self.init[1],self.init[1])],axis=1)
            j = np.argmin(J)
            M[1] = self.init[1][j]
            
        return M