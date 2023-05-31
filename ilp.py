import numpy as np
import sys

class GomoryCut:
    def __init__(self):
        self.A = None
        self.b = None
        self.c = None
        self.precision = 1e-6
        self.pdigits = 6

    def read_file(self,input_file):
        with open(input_file) as f:
            lines = f.readlines()
            n,m = lines[0].split(" ")
            n,m = int(n), int(m)
            self.b = np.array(list(map(int, lines[1].split(" "))))
            self.c = np.array(list(map(int, lines[2].split(" "))))
            self.A = np.zeros((m,n))
            for i in range(m):
                self.A[i] = np.array(list(map(int, lines[i+3].split(" "))))

        self.b = np.array(self.b).reshape(-1,1)
        self.c = np.array(-self.c).reshape(-1,1)

    def primal_simplex(self):
        while(np.min(self.tableau[0,1:]) < -self.precision):
    
            j  = 1
            while(j < self.cols and self.tableau[0,j] >= -self.precision):
                j+=1
            
            l = 0
            min_ratio = 0
            for i in range(1, self.rows):
                if(self.tableau[i,j] > self.precision):
                    ratio = self.tableau[i,0]/self.tableau[i,j]
                    if(l == 0 or ratio < min_ratio or self.basis[l-1] > self.basis[i-1] and ratio == min_ratio):
                        l = i
                        min_ratio = ratio

            if(l == 0):
                return False
    
            self.tableau[l] /= self.tableau[l,j]
            
            for i in range(self.rows):
                if(i != l):
                    self.tableau[i] -= self.tableau[l]*self.tableau[i,j]
            
            self.basis[l-1] = j

       
        return True

    def dual_simplex(self):
        while(np.min(self.tableau[1:,0]) < -self.precision):
           
            l = 1
            while(l < self.rows and self.tableau[l,0] >= -self.precision):
                l += 1
            
            j = 0
            min_ratio = 0
            for i in range(1, self.cols):
                if(self.tableau[l,i] < -self.precision):
                    ratio =  - self.tableau[0,i]/self.tableau[l,i]
                    if(j == 0 or ratio < min_ratio):
                        j = i
                        min_ratio = ratio
            
            if(j == 0):
                #print("Dual cost +inf")
                return False

            self.tableau[l] /= self.tableau[l,j]
            
            for i in range(self.rows):
                if(i != l):
                    self.tableau[i] -= self.tableau[l]*self.tableau[i,j]
            
            self.basis[l-1] = j

        return True

    def solve(self):
        self.m = self.A.shape[0]
        self.n = self.A.shape[1]

        self.rows = self.m+1
        self.cols = self.m+self.n+1

        self.tableau = np.zeros((self.m + 1, self.m + self.n + 1))
        self.tableau[0,0] = 0
        self.tableau[0,1:] = np.concatenate((self.c.T, np.zeros((1,self.m))), axis=1)
        self.tableau[1:,0] = self.b.flatten()
        self.tableau[1:,1:self.n + 1] = self.A
        self.tableau[1:,self.n+1:] = np.eye(self.m)

        self.basis = np.arange(self.n+1, self.n+self.m+1)

        while(1):
            t = self.primal_simplex()
            if t == False:
                return None
            
            k = 0
            for i in range(1, self.rows):
                if(self.tableau[i,0]%1 <= self.precision or 1-self.tableau[i,0]%1 <= self.precision):
                    continue
                else:
                    k = i
                    break

            if(k == 0):
                solution = np.zeros(self.n)
                for i in range(1, self.rows):
                    if(self.basis[i-1] <= self.n):
                        solution[self.basis[i-1]-1] = self.tableau[i,0]
                return list(map(round, solution.tolist()))

            self.tableau = np.concatenate((self.tableau, np.zeros(self.rows).reshape(-1,1)), axis=1)
            self.cols+=1
            self.tableau = np.concatenate((self.tableau, np.zeros(self.cols).reshape(1,-1)), axis=0)
            self.rows+=1

            for i in range(0, self.cols-1):
                if(self.tableau[k,i]%1 <= self.precision or 1-self.tableau[k,i]%1 <= self.precision):
                    continue
                else:
                    self.tableau[self.rows-1, i] = - (self.tableau[k,i]%1)
            self.tableau[self.rows-1, self.cols-1] = 1
                    
            self.basis = np.concatenate((self.basis,np.array([self.cols-1]).reshape(1)))

            t = self.dual_simplex()
            if t == False:
                return None

        return None    
        

def gomory(filename):
    gc = GomoryCut()
    gc.read_file(filename)
    t = gc.solve()

    return t

# print(gomory("test.txt"))