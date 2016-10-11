import random
import numpy as np
from math import pi, atan2, fabs, floor, ceil, sqrt

# distance
def dist(p1, p2):
  return sqrt(fabs(p1[0]-p2[0])**2 + fabs(p1[1]-p2[1])**2)

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

class bird:
  """
  Flies.  
  """
  def __init__(self, param, pos, phi):
    self.param = param
    self.pos   = pos
    self.phi   = phi

  def __getitem__(self, index):
    """
    Return position
    """
    return self.pos[index]

  def move(self):
    # update position given speed and angle
    self.pos = map(sum, zip(self.pos, [ self.param.speed*np.cos(self.phi),
                                        self.param.speed*np.sin(self.phi) ]))
    # implement bc
    self.pos = [ self.pos[0]%self.param.width,
                 self.pos[1]%self.param.height ] 

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

class param:
  """
  Parameters of the 
  """
  def __init__(self, width = 512, height = 512, N = 250, T=100, r = 10,
    n = 0.5, speed = 10, dt = 10):
    self.width = width # box width
    self.height = height # box height
    self.N = N # number of birds
    self.T = T # total time steps
    self.r = r # interaction radius
    self.n = n # noise strength
    self.speed = speed # linear speed of birds
    self.dt = dt # time interval for storage

class simulate:
  """
  Run the iterative loop
  """
  def __init__(self, param):
    self.param = param
    self.flock = flock(self.param)

  def run(self):
    for t in range(self.param.T):
      self.flock.move()
      if t==0: print('save initial flock with parameters: not implemented')
      if t%self.param.dt==0: print('store flock positions and angles: not implemented') #

class analysis:
  """
  Perform analysis of the simulation, for which we pass... stored flocks?
  """
  def __init__(self, flock):
    self.flock = flock

class plotter:
  """
  Make movie of the flock, also plot correlations and so on...
  """
  # def frame
  # def movie