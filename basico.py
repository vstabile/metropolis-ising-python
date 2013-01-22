from random import random, randint
from math import exp

class Ising(object):
  def __init__(self, L = 64, J = 1, estado = 'random'):
    self.L = L        # Sitios por lado
    self.N = L*L      # Total de sitios
    self.J = J        # Constante de interacao
    self.__inicializa_rede(estado)

  @property
  def T(self):        # Temperatura
    return self._T

  @T.setter
  def T(self, value): 
    self._T = float(value)
    self.__pre_calcula_acceptance_ratio()

  def __inicializa_rede(self, estado):
    if estado == 'up':
      self.rede = [[1 for j in xrange(self.L)] for i in xrange(self.L)]
    elif estado == 'down':
      self.rede = [[-1 for j in xrange(self.L)] for i in xrange(self.L)]
    elif estado == 'random':
      self.rede = [[1 for j in xrange(self.L)] for i in xrange(self.L)]
      for i in xrange(self.L):
        for j in xrange(self.L):
          if random() < 0.5: self.rede[i][j] = -1

  def __pre_calcula_acceptance_ratio(self):
    self.__acceptance_ratio = {
      2: exp(-4*self.J/self.T),
      4: exp(-8*self.J/self.T)
    }

  def itera(self):
    for n in xrange(self.N):
      i, j = randint(0, self.L - 1), randint(0, self.L - 1)
      soma = self.rede[(i-1)%self.L][j] + self.rede[i][(j-1)%self.L] \
           + self.rede[(i+1)%self.L][j] + self.rede[i][(j+1)%self.L]
      indice = self.rede[i][j]*soma
      delta_energia = 2*self.J*indice
      if delta_energia <= 0 or random() < self.__acceptance_ratio[indice]:
        self.rede[i][j] *= -1

  def magnetizacao(self):
    rede_flat = [spin for linha in self.rede for spin in linha]
    return float(sum(rede_flat))/self.N

  def energia(self):
    soma = 0
    for i in xrange(self.L):
      for j in xrange(self.L):
        soma += self.rede[i][j]*(self.rede[(i-1)%self.L][j] + self.rede[i][(j-1)%self.L])
    return -self.J*float(soma)/self.N