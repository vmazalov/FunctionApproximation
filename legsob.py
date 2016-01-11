from sympy import Rational
from operator import add
from math import sqrt, exp

def addl(L1,L2):
    if len(L1)>len(L2):
        L1,L2 = L2,L1
    return [L1[i]+L2[i] for i in range(len(L1))] + L2[len(L1):]
    
def cmull(c,L):
    return map(lambda x: c*x, L)
    
def shiftl(L):
    return [0] + L

def differentiate(L):
    dL = []
    for i in range(1,len(L)):
        dL.append(L[i]*i)
    return dL
        
def integrate01(L):
    return sum([L[i]*Rational(1,i+1) for i in range(len(L))])
    
def evaluate(L,x):
    t = 1
    s = 0
    for i in range(len(L)):
        s += L[i] * t
        t *= x
    return s

class LegendreSobolevPolynomials(object):
    n = 17
    mu = [1,1]
    
    oldn = 0
    oldmu = []

    def __setattr__(self,name,value):
        if name == 'n':
            self.oldn = self.n
        if name == 'mu':
            self.oldmu = self.mu
        object.__setattr__(self,name,value)
        
    def __getattribute__(self,name):
        if name not in ['n','mu','oldn','oldmu']:
            if self.n != self.oldn or self.mu != self.oldmu:
                self.oldn = self.n
                self.oldmu = self.mu
                self.update()
        return object.__getattribute__(self,name)

    def update(self):
        self.LSP = self.LS()
        self.dLSP = map(differentiate, self.LSP)
        self.d2LSP = map(differentiate,self.dLSP)
        self.LSPdiff = map(lambda i: addl(self.LSP[i], cmull(-self.mu[1],self.d2LSP[i])), range(len(self.LSP)))
        self.LSPNorms = map(self.LSnorm, self.LSP)
        self.LSPN = [map(lambda x: x/sqrt(self.LSPNorms[i]), self.LSP[i]) for i in range(len(self.LSP))]
        self.LSPdiffN = [map(lambda x: x/sqrt(self.LSPNorms[i]), self.LSPdiff[i]) for i in range(len(self.LSPdiff))]
        self.e0 = map(lambda L: evaluate(L,0), self.dLSP)
        self.e0N = map(lambda i: self.e0[i]/sqrt(self.LSPNorms[i]), range(len(self.e0)))
        self.e1 = map(lambda L: evaluate(L,1), self.dLSP)   
        self.e1N = map(lambda i: self.e1[i]/sqrt(self.LSPNorms[i]), range(len(self.e1)))
        
    def LS(self):
        """ Computes Legendre-Sobolev polynomials up to degree `n'
            These polynomials are orthogonal w.r.t. the Sobolev inner product
            <f,g> = int(f(x)g(x)+mu[1]f'(x)g'(x)+...+mu[s]f^(s)(x)g^(s)(x),x=0..1)
            mu[0] must be equal to 1
            
            Reference: W.Gautschi, M.Zhang, 
                       Computing orthogonal polynomials in Sobolev spaces, 
                       Numerische Mathematik, Volume 71, Issue 2, Pages 159 - 183, 1995.
            
            Function 'leg_sob_polys' is stored in 'leg_sob.m', to be used in 'expansions_proc.mw'
            
            Implemented in Python by: Oleg Golubitsky (December, 2008)
            email: oleg.golubitsky@gmail.com
        """
        n = self.n
        mu = self.mu
        s = len(mu) - 1
           
        def Tau(j,k):
            if j<0 or j>=k:
                return 0
            else:
                return tau[j,k]
                
        def delta(j,k):
            if j==k:
                return 1
            else:
                return 0
                
        def nu(k):
            return Rational(1,k+1)
            
        sigma = {}
        for l in range(2*n):
            sigma[0,l] = nu(l)
        beta = {}
        beta[0,0] = sigma[0,1]/sigma[0,0];
    
        tau = {}
        tau[0,1] = 1
        for k in range(1, 2*n-2):
            for j in range(k+1):
                tau[j,k+1] = Tau(j-1,k)+delta(j,k)
                
        mus = {}
    
        for i in range(1,s+1):
            for u in range(i+1):
                for v in range(i+1):
                    for k in range(2*n-1):
                        if u == 0: 
                            if v == 0: 
                                mus[0,k,i,u,v] = nu(k)*mu[i];
                            else:
                                mus[0,k,i,u,v] = sum([Tau(j,k)*mus[0,j,i,0,v-1] for j in range(v-1,k)])
                        else:
                            mus[0,k,i,u,v] = 0
    
        for k in range(1,n):
            for l in range(k,2*n-k):
                sigma[k,l] = sigma[k-1,l+1] + sum([i*(mus[k-1,l,i,i-1,i]-mus[k-1,l,i,i,i-1]) for i in range(1,s+1)]) - sum([beta[j,k-1]*sigma[k-1-j,l] for j in range(k)])
            beta[0,k] = sigma[k,k+1]/sigma[k,k]-sigma[k-1,k]/sigma[k-1,k-1];
            for j in range(k-1,0,-1):
                beta[k-j,k] = sigma[j,k+1]/sigma[j,j]-sigma[j-1,k]/sigma[j-1,j-1]-sum([beta[l-j,l]*sigma[l,k]/sigma[l,l] for l in range(j,k)])
            beta[k,k] = sigma[0,k+1]/sigma[0,0]-sum([beta[l,l]*sigma[l,k]/sigma[l,l] for l in range(k)])
            if k < n - 1:
                for i in range(1,s+1):
                    for u in range(i+1):
                        for v in range(i+1):
                            for l in range(k,2*n-k-1):
                                mus[k,l,i,u,v] = mus[k-1,l+1,i,u,v] + u*mus[k-1,l,i,max(u-1,0),v] - v*mus[k-1,l,i,u,max(v-1,0)] - sum([beta[j,k-1]*mus[k-1-j,l,i,u,v] for j in range(k-u)])
        pi = [0 for i in range(n+1)]
        pi[0] = [1]
            
        for k in range(n):            
            L = [cmull(beta[j,k],pi[k-j]) for j in range(k+1)]
            al = reduce(addl, L)
            al = cmull(-1,al)
            pi[k+1] = addl(shiftl(pi[k]), al)
        return pi

    def LSnorm(self,L):
        c = 1
        d = len(L) - 1
        s = 0
        for i in range(d,max(-1,d-len(self.mu)),-1):
            sL = [0 for j in range(i)] + L
            s += integrate01(sL) * c * self.mu[d-i]
            L = differentiate(L)
            c *= i
        return s

    def expansion(self,X,T):
        M = moments(X,T,self.n)
        part1 = momentsToExpansion(M, self.LSPdiffN)
        part2 = map(lambda i: X[-1] * self.e1N[i] - X[0] * self.e0N[i], range(len(self.e0N)))
        return addl(part1,cmull(self.mu[1],part2))

    
def arclen(P):
    """ Computes the list of arc lengths for points in P.
        For example: arclen([[1,1],[4,5],[7,9]]) = [0.0,5.0,10.0]
    """
    dim = len(P[0])
    s = 0.0
    L = [s]
    for i in range(len(P)-1):
        s += sqrt(sum([(P[i][j] - P[i+1][j])**2 for j in range(dim)]))
        L.append(s)
    return L

def moments(X,T,n):
    """ Computes the list of the first (n+1) moments for X(T)
        scaled to [0,1]
    """
    TP = [float(t) for t in T]
    M = []
    dT = 1.0 / (T[-1] - T[0])
    t = 1.0
    for k in range(n+1):
        s = (X[-2] + X[-1]) * TP[-1] - (X[0] + X[1]) * TP[0]
        for i in range(1,len(T)-1):
            s += (X[i-1] - X[i+1]) * TP[i]
        t *= dT
        M.append(s * t / (2.0*k+2))
        for i in range(len(T)):
            TP[i] *= T[i]
    return M

def momentsToExpansion(M,P):
    E = []
    for i in range(len(M)):
        s = 0.0
        for j in range(i+1):
            s += M[j] * P[i][j]
        E.append(s)
    return E

def scaleTranslate(Ex, Ey):
    E = reduce(add,[[Ex[i],Ey[i]] for i in range(1,len(Ex))])
    w = sum(map(lambda x: x*x, E))
    if w > 0:
        E = cmull(1/w, E)
    return E

