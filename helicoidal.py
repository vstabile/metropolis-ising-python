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
      self.rede = [1]*self.N
    elif estado == 'down':
      self.rede = [-1]*self.N
    elif estado == 'random':
      self.rede = [1]*self.N
      for i in xrange(self.N):
        if random() < 0.5: self.rede[i] = -1

  def __pre_calcula_acceptance_ratio(self):
    self.__acceptance_ratio = {
      2: exp(-4*self.J/self.T),
      4: exp(-8*self.J/self.T)
    }

  def itera(self):
    for n in xrange(self.N):
      i = randint(0, self.N - 1)
      soma = self.rede[(i-1)%self.N] + self.rede[(i-self.L)%self.N] + \
             self.rede[(i+1)%self.N] + self.rede[(i+self.L)%self.N]
      indice = self.rede[i]*soma
      delta_energia = 2*self.J*indice
      if delta_energia <= 0 or random() < self.__acceptance_ratio[indice]:
        self.rede[i] *= -1

  def magnetizacao(self):
    return float(sum(self.rede))/self.N

  def energia(self):
    soma = 0
    for i in xrange(self.N):
      soma += self.rede[i]*(self.rede[(i-1)%self.N] + self.rede[(i-self.L)%self.N])
    return -self.J*float(soma)/self.N