import random
import numpy as np

class parameters:
  """
  Parameters of the simulation, including the 
  """
  def __init__(self, nettype = 'sf-sf', sat ='cubical', stress = 'piecewise',
               N = 100, T = 1000, dt = 100, g = 0.1, s = 1.0, c = 0.1,
               alpha = 0.1, Nt = 1, y0 = 0, M0 = 0.1, mu = 0.1, eps = 0.1):
    self.nettype = nettype # network type
    self.sat = sat # type of saturating function
    self.env = env # type of environment
    self.N = N # size of the cell network
    self.T = T # total time steps
    self.dt = dt # storage steps
    self.g = g # network gain
    self.s = s # saturation sensitivity
    self.c = c # sparseness of trait vectors b
    self.alpha = alpha # 
    self.Nt = Nt # number of traits
    self.y0 = y0 # target phenotype
    self.M0 = M0 # stress amplitude
    self.mu = mu # stress sensitivity
    self.eps = eps # stress range

class environment: 
  """
  Contains the trait vectors as well as the fitness landscape
  """
  def __init__(self, parameters):
    self.b  = np.random(parameters.N, sigma) # traits

  def __getitem__(self, tup):
    """
    Default return are traits
    """
    i, j = tup
    return self.b[i,j]

  def __setitem__(self, tup, value):
    """
    Default input are traits
    """
    i, j = tup
    self.b[i,j] = value

  def M(self, parameters):
    """
    Calculate stress
    """
    cases = {'piecewise': M_pw, 'quartic': M_q4, 'quadratic': M_q2,
             'sigmoidal': M_sig }                 
    return cases[parameters.stress](self, parameters)

  def M_pw(self, parameters):
    """
    Piecewise linear stress
    """
    return 




class cell: 
  """
  Contains all the state variables of the cell. It's a combination of the
  'chemicals' class, which contains the chemical state of the cell (x, phi)
  and the 'interactions' class, which contains the interactions (Jij, T) and
  performs basic operations on them (eigenvalues, ...)
  """
  def __init__(self, parameters):
    self.net  = interactions(parameters) # interaction network, J
    self.chem = chemicals(parameters) # chemical composition, x
    self.y = np.zeros(parameters.M) # phenotypes, y

  def __getitem__(self, index):
    """
    Default return are chemical concentrations
    """
    return self.chem[index]

  def __setitem__(self, index, value):
    """
    Default input is the chemical concentration
    """
    self.chem[index] = value


class simulation: # first run parameters, then cell, then environment, then simulate
  """
  Run the iterative loop
  """
  def __init__(self, parameters, cell, env):
    self.parameters = parameters

  def evolve(self):
    # update chemicals
    cell[:]= cell[:]+(cell.net[:,:]*cell.phi()-cell[:])*self.parameters.dt +
             np.rand.rand()*self.parameters.dt**0.5

    # update phenotype
    cell.y = env[:,:]*cell[:]

    # update connectivities with Ito rule
    self.net[:,:] = np.sqrt(env.M(cell.y)*parameters.D)*np.rand.rand(size(J))

  def run(self):
    for t in range(self.param.T):
      self.flock.move()
      if t==0: print('save initial flock with parameters: not implemented')
      if t%self.param.dt==0: print('store flock positions and angles: not implemented') #




# TOOLS TO BE USED, MAYBE ANOTHER FILE?
class chemicals:
  """
  This class contains the chemical concentrations as well as their saturating
  functions
  """
  def __init__(self, parameters):
    self.x = np.random(parameters.N)

  def __getitem__(self, index):
    """
    By default this class returns x
    """
    return self.x[i]

  def __setitem__(self,index, value):
    """
    And by default this class sets x as input
    """
    self.x[index] = value

  def phi(self):
    """
    This function reads the type of saturation and calculates it for x
    """
    cases = {'sigmoidal': sat_sig, 'cubical': sat_cub}                 
    return cases[parameters.saturation](self, parameters)

  def sat_sig(self, parameters):
    """
    Sigmoidal saturation function
    """
    a = 1+1
    return a

  def sat_cub(self, parameters):
    """
    Cubical saturation function
    """
    return self.x-(parameters.sat/(3*2))*self.x**3


class interactions:
  """
  Taking the parameters as input this class generates an interaction matrix
  and calculates the topology 
  """
  def __init__(self,parameters):
    self.W0 = self.generate(parameters)
    self.T  = scipy.sign(self.W0)
    self.J  = self.W0

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

  def generate(self, parameters):
    """
    This function reads the type of network desired from the parameters and
    calls the corresponding function to generate it.
    """
    cases = {'sf-sf': net_sfsf, 'sf-exp': net_sfexp, 'exp-sf': net_expsf,
             'sf-exp': net_sfexp, 'sf-bin': net_sfbin, 'bin-sf': net_binsf}                 
    return cases[parameters.nettype](parameters)

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
  def net_sfsf(parameters):
    """
    Generate interaction networks scale-free -- scale-free
    """
    W = np.random.random(parameters.N,parameters.N)
    return W

  def net_sfexp(parameters):
    """
    Generate interaction networks scale-free -- exponential
    """
    W = np.random.random(parameters.N,parameters.N)

    return W



class potential:
  """
  Taking the parameters as input this class generates an interaction matrix
  and calculates the topology 
  """
  def __init__(self,parameters):
    self.W = self.generate_W(parameters)
    self.T = topology of T
    self.J = strengths of J

  def generate_W(self,parameters):
    switcher = {'sf-sf': net_sfsf, 'sf-exp': net_sfsf, 'exp-sf': net_expsf,
                'sf-exp': net_sfexp, 'sf-bin': net_sfbin, 'bin-sf': net_binsf}                 
    return switcher[parameters.nettype](parameters)

  def net_sfsf(parameters):
    W = np.random.random(parameters.N,parameters.N)

    return W



## from old softwareeee!!!!


class analysis:
  """
  Perform analysis of the simulation, for which we pass... stored flocks?
  """
  def __init__(self, flock):
    self.flock = flock

class bucket_grid:
  """
  Fixed radius nearest neighbour search using a grid of buckets of size r.
  Width and height are needed to wrap around (periodic boundary conditions).
  """
  def __init__(self, points, param):
    self.param = param
    self.n = self.param.width / self.param.r # # of horizontal buckets
    self.m = self.param.height / self.param.r # # of vertical buckets 
    self.buckets = {} # dictionary of buckets -> points

    # create dictionary for given points
    for p in points:
      self.buckets.setdefault(self.get_index(p), []).append(p)
 
  def get_index(self, p):
    # returns bucket coordinates of point
    return ( int(floor(p[0]/self.param.width*self.n)),
             int(floor(p[1]/self.param.height*self.m)) )

  def neighbours(self, p):
    # position of the central bucket
    i, j = self.get_index(p)
    # this is the number of adjacent buckets we need to check
    cx = int(ceil(float(self.param.r)/self.param.width*self.n))
    cy = int(ceil(float(self.param.r)/self.param.height*self.m))
    neighbours = []
    # check all neighbouring buckets
    for a in range(-cx, 1+cx):
      for b in range(-cy, 1+cy):
        # add points
        neighbours += filter(
          lambda q: dist(p, q) < self.param.r,
          self.buckets.setdefault( ( (i+a)%self.n, (j+b)%self.m ), [])
          )
    return neighbours

class flock:
  """
  Many of those things together.
  """
  # define a class param with all parameters
  def __init__(self, param):
    self.param = param

    # create birds
    self.birds = [ bird(self.param, # can change to param?
                        [ random.random()*self.param.width, random.random()*self.param.height ], 
                        2*pi*random.random() 
                       ) for i in range(self.param.N) ]
    # self.birds_old = copy.deepcopy(self.birds) 

  def move(self):
    # create the buckets
    grid = bucket_grid(self.birds, self.param)

    # update the angles
    for b in self.birds:
      sin_tot = 0.
      cos_tot = 0.
      # loop over neighbours
      neighbours = grid.neighbours(b) 
      for n in neighbours:
          sin_tot += np.sin(n.phi)
          cos_tot += np.cos(n.phi)
      # update ONE SHOULD UPDATE ALL ELEMENTS AT THE SAME TIME, I need birds_o
      b.phi = atan2(sin_tot, cos_tot) + self.param.n/2.*(1-2.*random.random())

    # move them
    for b in self.birds:
      b.move()

class plotter:
  """
  Make movie of the flock, also plot correlations and so on...
  """
  # def frame
  # def movi