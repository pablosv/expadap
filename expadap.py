"""
This file contains the core classes to numerically study exploratory
adaptation. These are environment, cell, and simulation. There are other
classes that are auxiliary, such as chemicals, interactions, and xxx. These
should perhaps go into a tools file. Such a tools file could also contain
analysis functions, as well as storage functions/classes.
"""

import random
import numpy as np
import scipy.sparse as sp

class environment: 
  """
  Contains the trait vectors as well as the fitness function
  """
  def __init__(self, stress = 'tanh', M0 = 0.1, mu = 0.1, eps = 0.1, N = 100,
               c = 0.2, alpha = 100., y0 = 0, g = 0.1):
    self.stress = stress # type of environmental fitness
    self.M0 = M0 # stress amplitude
    self.mu = mu # stress sensitivity
    self.eps = eps # stress range
    self.N = N # size of the cell network
    self.c = c # sparseness of trait vectors b
    self.alpha = alpha # trait variance factor
    self.g = g # network gain
    self.y0 = y0 # target phenotype
    self.b = self.traits() # traits

  def __getitem__(self, index):
    """
    Default return are traits
    """
    return np.array((self.b.todense()))[0][index]

  def __setitem__(self, index, value):
    """
    Default input are traits
    """
    self.b[0,index] = value

  def traits(self):
    """
    Generates a trait vector b with sparseness c. Its not-null elements are
    sampled from a normal distribution with variance var that scales with N.
    """
    self.b = sp.rand(1, self.N, density = self.c, format='csr')
    notnull = len(self.b.data)
    var = self.alpha / (self.c * self.N * self.g**2)
    self.b.data[:] = np.random.normal(0, var, notnull)
    return self.b

  def M(self, y):
    """
    Calculate stress for the corresponding function
    """
    cases = {'piecewise-constant': self.M_pwc, 'piecewise-linear': self.M_pwl,
             'quartic': self.M_q4, 'quadratic': self.M_q2, 'tanh': self.M_th }                 
    return cases[self.stress](y)

  # environmental fitness functions
  def M_pwc(self, y):
    """
    Piecewise constant stress
    """
    y = y - self.y0
    M = np.piecewise(y, [np.abs(y)>self.eps, np.abs(y)<=self.eps],
                        [self.M0, 0.])
    return M

  def M_pwl(self, y):
    """
    Piecewise linear and flat stress
    """
    y = y - self.y0
    M = np.piecewise(y, [np.abs(y)>self.eps, np.abs(y)<=self.eps],
                        [lambda y: np.abs(y)-self.eps, 0.])
    return M
  
  def M_th(self, y):
    """
    Hyperbolic tangent stress
    """
    y = y - self.y0
    M = (self.M0/2.) * (1+np.tanh((np.abs(y)-self.eps)/self.mu))
    return M

  def M_q4(self, y):
    """
    Quartic stress
    """
    y = y - self.y0
    M = np.piecewise(y,[np.abs(y)>self.eps*2**.25, np.abs(y)<=self.eps*2**.25],
                       [self.M0, lambda y: (self.M0/2)*(y/self.eps)**4])
    return M

  def M_q2(self, y):
    """
    Quadratic stress
    """
    y = y - self.y0
    M = np.piecewise(y,[np.abs(y)>self.eps*2**.5, np.abs(y)<=self.eps*2**.5],
                       [self.M0, lambda y: (self.M0/2)*(y/self.eps)**2])
    return M

class cell: 
  """
  Contains the state variables of the cell. It's a combination of the
  'chemicals' class, which contains the chemical state of the cell (x, phi)
  and the 'interactions' class, which contains the interactions (Jij, T) and
  performs basic operations on them (eigenvalues, ...)
  """
  def __init__(self, nettype = 'sf-sf', sat ='cubical', N = 100, g = 0.1,
               s = 1., beta = 1., gamma = 1.2, a = 1., p = .01, rho =0.1):
    self.nettype = nettype # network type
    self.sat = sat # type of saturating function
    self.N = N # size of the cell network
    self.g = g # network gain
    self.s = s # saturation sensitivity
    self.beta = beta # scale of exponential distribution
    self.gamma = gamma # scale of scale-free distribution
    self.p = p # scale of binomial distribution
    self.rho = rho # density of sparse matrix

    graph  = interactions(self) # interaction network
    chem = chemicals(self) # chemical composition, x

    self.T = graph.T
    self.W = graph.W
    self.J = graph.J
    self.x = chem.x
    self.y = 0 # phenotypes, y

  def __getitem__(self, index):
    """
    Default return are chemical concentrations
    """
    return self.x[index]

  def __setitem__(self, index, value):
    """
    Default input is the chemical concentration
    """
    self.x[index] = value

class simulation:
  """
  Run the iterative loop
  """
  def __init__(self, env, cell, T = 10, dt = .1):
    self.T = T # total time steps
    self.dt = dt # storage steps
    self.env = env
    self.cell = cell
    if (env.N,env.g) != (cell.N, cell.g):
      print 'Cell and environment not compatible'
      quit()

  def update(self):
    # update chemicals
    self.cell[:] = self.cell[:]+(np.sin(self.cell.W.dot(self.cell.x))-self.cell[:])*self.dt
              #+ np.rand.rand()*self.parameters.dt**0.5

    # update phenotype
    self.cell.y = np.dot(self.env[:], self.cell[:])

    # update interactions
    self.cell.J = np.sqrt(self.env.M(self.cell.y))*np.random.normal(0, 1, self.cell.J.shape[0])
    self.cell.W = sp.csr_matrix((self.cell.J,self.cell.T.nonzero()),
                 self.cell.T.shape)


  def run(self):
    for t in np.arange(0,self.T,self.dt):
      if t==0: print('save parameters and initial state')
      if t%self.dt==0: print('save state')
      self.update()
      print cell.x[0]




# TOOLS TO BE USED, MAYBE ANOTHER FILE?
class interactions:
  """
  Taking the parameters as input this class generates an interaction matrix
  and calculates the topology 
  """
  def __init__(self, cell):
    self.T = self.generate(cell)
    self.J = np.random.normal(0, 1, np.count_nonzero(self.T))
    self.W = sp.csr_matrix((self.J,self.T.nonzero()),
              self.T.shape)

  def __getitem__(self, tup):
    """
    By default the interactions class returns the Jij
    """
    i, j = tup
    return self.J[i,j]

  def __setitem__(self, tup, item):
    """
    And by default the interaction class gets Jij as input
    """
    i, j = tup
    self.J[i,j] = item

  def generate(self, cell):
    """
    This function reads the type of network desired from the parameters and
    calls the corresponding function to generate it.
    """
    cases = {'sf-sf': self.net_sfsf, 'sf-exp': self.net_sfexp,
             'exp-sf': self.net_expsf, 'sf-bin': self.net_sfbin,
             'bin-sf': self.net_binsf, 'erdos-reny': self.net_er}    
    return cases[cell.nettype](cell)

  # Functions that calculate network properties
  def eigenvalues(self):
    """
    Calculate the eigenvalues of the interaction matrix
    """
    u = np.eigenvalues(self.J)

    return u

  def property1(self):
    """
    This function calculates a property of the network
    """  

    return 'Nothing for now'

  # Functions that generate interaction networks
  def net_sfsf(self, cell):
    """
    Generate interaction networks scale-free -- scale-free
    """
    W = sp.rand(cell.N,cell.N,density=0.1) #np.random.pareto(.1,(2,2))
    return W

  def net_sfexp(self, cell):
    """
    Generate interaction networks scale-free -- exponential
    """ 
    # create graph class
    # for i<100 attempts:
    #   G.in  = expo(length)
    #   G.out = expo(length)
    #   if G.test() == True: break
    # if i=99: print 'No network found in 100 attempts' 
    # return

    W = np.random.exponential(.1,(2,2))
    return 1

  def net_expsf(self, cell):
    """
    Generate interaction networks scale-free -- exponential
    """
    W = 1
    return 1

  def net_sfbin(self, cell):
    """
    Generate interaction networks scale-free -- binomial
    """
    din  = np.round(np.random.pareto(cell.gamma, cell.N)).astype(int)
    dout = np.random.binomial(cell.N, cell.p, cell.N)
    G = graph(din = din, dout = dout)
    
    print 'Generating degree sample'
    while din.sum()!=dout.sum() or G.test()==False:
      din  = np.round(np.random.pareto(cell.gamma, cell.N)).astype(int)
      dout = np.random.binomial(cell.N, cell.p, cell.N)
      # print din.sum()-dout.sum()
      G.__init__(din = din, dout = dout)

    print 'Constructing topology'
    G.construct()
    # perform test bw bds and G.T


    return G.T

  def net_binsf(self, cell):
    """
    Generate interaction networks binomial -- scale-free
    """
    din = np.random.binomial(cell.N, cell.p, cell.N)
    dout = np.round(np.random.pareto(cell.gamma, cell.N)).astype(int)
    G = graph(din = din, dout = dout)
    
    print 'Generating degree sample'
    while din.sum()!=dout.sum() or G.test()==False:
      din = np.random.binomial(cell.N, cell.p, cell.N)
      dout = np.round(np.random.pareto(cell.gamma, cell.N)).astype(int)
      # print din.sum()-dout.sum()
      G.__init__(din = din, dout = dout)

    print 'Constructing topology'
    G.construct()

    return G.T

  def net_er(self, cell):
    """
    Generate interaction networks binomial -- binomial (i.e., erdos-reny)
    """
    din  = np.random.binomial(cell.N, cell.p, cell.N)
    dout = np.random.binomial(cell.N, cell.p, cell.N)
    G = graph(din = din, dout = dout)
    
    print 'Generating degree sample'
    while din.sum()!=dout.sum() or G.test()==False:
      din  = np.random.binomial(cell.N, cell.p, cell.N)
      dout = np.random.binomial(cell.N, cell.p, cell.N)
      G.__init__(din = din, dout = dout)

    print 'Constructing topology'
    G.construct()

    return G.T

class graph:
    """
    This is the graph class. It constructs directed graphs given the in and out
    degree distributions. It also provides the weight of the graph. 
    """
    def __init__(self, din=np.array([2,2,1,1,1]), dout=np.array([2,1,3,1,0])):
      self.din  = din  # in connections
      self.dout = dout # out connections
      self.labels = 1+np.arange(len(din)) # nodes' labels
      self.D = self.normal(np.array((self.din, self.dout, self.labels))) # BDS
      self.T = np.zeros((len(din),len(din))) # sparse matrix with zeros.
      self.w = 1. # weight of T

    def construct(self):
      """
      We follow the algorithm in Kim et al, using their numbering of the steps
      """
      self.T = np.zeros((len(self.din),len(self.din))) # (0.0) reset topology
      D = self.D # (0.1) working BDS
      if self.test() == False: return 'Non-graphable'# (0.2) test
      nout = np.count_nonzero(D[1,:]) # (0.3) number of non-null out-nodes
      
      for k in np.arange(nout):
        D = self.normal(D) # (7) set D to normal order
        i = np.argmax(D[1,:]>0) # (1.0) working index
        a = D[2,i] # (1.1) working label
        stubs = D[1,i] # (1.2) number of stubs
        chi = np.append(a, D[2,np.where(D[0,:]==0)[0]]) # (2) forbidden labels   
        
        for s in np.arange(stubs):
          if np.count_nonzero(D[1,:])==1 or np.count_nonzero(D[0,:])==1: # last node
            bs = D[2,np.where(D[0,:]==1)[0]]
            self.T[a-1,bs-1] = 1 # wire
            # print 'yay'
            break
          if s>0: D = self.normal(D) # (6.0)  set D to normal order
          i = np.argmax(D[2,:]==a) # (6.1) index of a
          A = self.allowed(D, i, a, chi) # (3) graphicable labels
          # print str(A)
          # print str(D)
          b = np.random.choice(A) # (4.1) choose an in-node by its label
          self.T[a-1,b-1] = 1 # (4.2) make link between labeled nodes
          j = np.argmax(D[2,:]==b) # (4.3) index of b 
          D[0,j], D[1,i] = D[0,j] - 1, D[1,i] - 1 # (4.4) calculate residual
          #self.w *= np.size(np.array(A)) # (4.5) path weight
          chi = np.append(chi,b) # (5) b is now forbidden
      return None

    def test(self):
      """
      Fulkerson-Ryser test of graphability for a given in/out degree sequences
      """
      if (self.din>len(self.din)-1).any(): print 'E1'; return False # din < N
      if (self.dout>len(self.din)-1).any(): print 'E2'; return False # dout < N
      if self.din.sum()!=self.dout.sum(): print 'E3'; return False # sums
      for k in np.arange(0,len(self.din)-1):
        if self.LpKp(k, self.D)>0: print 'E4'; return False # Lk > Rk
      return True

    def allowed(self, D, i, a, chi):
      """
      For the current D, a  and chi find the set of in-nodes (j) that preserve
      graphicality. This is step (3) in Kim et al.
      """
      L = np.setdiff1d(D[2,:],chi)[:D[1,i]] # (3.1) compute L
      R = np.setdiff1d(D[2,:],np.append(chi,L)) # (3.2) compute R
      Dp = np.zeros(np.shape(D)) # allocate D'
      Dp[:] = D[:] # (3.3.1) create D'
      Dp[0,np.in1d(D[2,:],L[:-1])] = Dp[0,np.in1d(D[2,:],L[:-1])] -1 # (3.3.2)
      Dp[1,i] = 1 # (3.3.3) set out nodes
      Dp = self.normal_in(Dp) # (3.4) rearrange Dp
      i = np.argmax(Dp[2,:]==a) # recalculate i
      kini = 0 if i!=0 else 1 # (3.5.1) set initial k
      k0, N = -1, np.shape(Dp)[1] # (3.5.2) set reference k0, and size
      for k in range(kini,N): # (3.5)
        if self.LpKp(k, Dp) == 0: k0 = k; break
      if k0==-1 or k0==N-1: return np.setdiff1d(D[2,:],chi) # (3.5) return A
      # print('E3');
      qp = k0+1+np.argmax(np.in1d(Dp[2,:],R)[k0+1:]==True) # (3.6.1) q'
      q = np.argmax(D[2,:]==Dp[2,qp]) # (3.6.2) q
      return np.setdiff1d(D[2,:q],chi) # (3.7) return A

    def normal(self, D):
      """
      Put the bi-degree sequence D in normal order. Since 'sort' uses
      increasing order we need to reverse the array
      """
      Dr = D[:,::-1]
      D = Dr[:,np.lexsort((Dr[1,:],Dr[0,:]))][:,::-1]
      return D

    def normal_in(self, Dp):
      """
      Put the bi-degree sequence D in normal order only considering the in
      nodes. This is faster than self.normal
      """
      Dp = Dp[:,np.argsort(-Dp[0,:])]
      return Dp

    def LpKp(self, k, D):
      """
      Test the L'(k) = R'(k) condition (equations 4 and 5).
      """
      Lp = D[0,:k+1].sum()
      Rp = np.minimum(np.ones(k+1)*k,D[1,:k+1]).sum()+np.minimum(np.ones(len(D[1,k+1:]))*(k+1),D[1,k+1:]).sum()
      return Lp - Rp

class chemicals:
  """
  This class contains the chemical concentrations as well as their saturating
  functions
  """
  def __init__(self, cell):
    self.x = np.random.rand(cell.N)

  def __getitem__(self, index):
    """
    By default this class returns x
    """
    return self.x[index]

  def __setitem__(self,index, value):
    """
    And by default this class sets x as input
    """
    self.x[index] = value

  def phi(self, cell):
    """
    This function reads the type of saturation and calculates it for x
    """
    cases = {'sigmoidal': self.sat_sig, 'cubical': self.sat_cub}                 
    return cases[cell.sat](self, cell)

  def sat_sig(self, cell):
    """
    Sigmoidal saturation function
    """
    a = 1+1
    return a

  def sat_cub(self, cell):
    """
    Cubical saturation function
    """
    return self.x-(cell.s/(3*2))*self.x**3
