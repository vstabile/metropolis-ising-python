from random import random, randint, seed
from math import exp, ceil
from multiprocessing import cpu_count, Pool

def unwrap_method(arg, **kwarg):
   return Ising._itera_dominio(*arg, **kwarg)

class Ising(object):
  def __init__(self, L = 64, J = 1, estado = 'random'):
    self.L = L        # Sitios por lado
    self.N = L*L      # Total de sitios
    self.J = J        # Constante de interacao
    global pool
    pool = Pool(processes = cpu_count())
    self.nucleos = cpu_count() # Numero de nucleos de processamento
    self.__inicializa_rede(estado)
    self.__decompoe_dominios()

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
    self.__dominios = pool.map(unwrap_method, zip([self]*len(self.__dominios), self.__dominios))
    self.__desloca_dominios()

  def _itera_dominio(self, dominio):
    seed()
    n = len(dominio) # Sitios no dominio
    ultimo_sitio_interno = n - self.L  - 1
    for k in xrange(n):
      i = randint(self.L, ultimo_sitio_interno)
      soma = dominio[i-1] + dominio[i+1] + dominio[i-self.L] + dominio[i+self.L]
      indice = dominio[i]*soma
      delta_energia = 2*self.J*indice
      if delta_energia <= 0 or random() < self.__acceptance_ratio[indice]:
        dominio[i] *= -1
    return dominio

  def magnetizacao(self):
    return float(sum([sum(dominio) for dominio in self.__dominios]))/self.N

  def energia(self):
    self.__recompoe_dominios()
    soma = 0
    for i in xrange(self.N):
      soma += self.rede[i]*(self.rede[(i-1)%self.N] + self.rede[(i-self.L)%self.N])
    return -self.J*float(soma)/self.N

  def __decompoe_dominios(self):
    n = int(ceil(float(self.L)/self.nucleos))*self.L # Sitios por dominio
    self.__dominios = [self.rede[x:x+n] for x in xrange(0, self.N, n)]

  def __desloca_dominios(self):
    numero_dominios = len(self.__dominios)
    for x in xrange(0, numero_dominios - 1):
      self.__dominios[x] += self.__dominios[(x+1)%numero_dominios][0:self.L]
    for x in xrange(0, numero_dominios - 1):
      del self.__dominios[x][0:self.L]

  def __recompoe_dominios(self):
    self.rede = []
    for dominio in self.__dominios:
      self.rede += dominio