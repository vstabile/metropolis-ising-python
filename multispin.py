from random import random, randint, seed
from math import exp, ceil, log
from multiprocessing import cpu_count, Pool
from sys import maxint
from gmpy import popcount

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
    self.bits = maxint.bit_length()
    self.__inicializa_rede(estado)
    self.__decompoe_dominios()

  @property # Getter para temperatura
  def T(self):
    return self._T

  @T.setter # Setter para temperatura
  def T(self, value):
    self._T = float(value)
    self.__pre_calcula_q()

  def __inicializa_rede(self, estado):
    if estado == 'random':
        self.rede = [randint(0, maxint) for x in xrange(self.N)]
    elif estado == 'up':
        self.rede = [maxint for x in xrange(0, self.N)]
    elif estado == 'down':
        self.rede = [0 for x in xrange(0, self.N)]

  def __pre_calcula_q(self):
    self.__q = {
      0: 1 - 2*exp(-8*self.J/self.T),
      1: 1 - 2*(exp(-4*self.J/self.T) - exp(-8*self.J/self.T))/(1 - exp(-8*self.J/self.T))
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
      # Calcula termos da expressao mestre
      a1 = dominio[i] ^ dominio[i+1]
      a2 = dominio[i] ^ dominio[i-1]
      a3 = dominio[i] ^ dominio[i+self.L]
      a4 = dominio[i] ^ dominio[i-self.L]
      R1 = a1 | a2 | a3 | a4
      R2 = ((a1 | a2) & (a3 | a4)) | ((a1 & a2) | (a3 & a4))
      # Calcula palavras aleatorias
      r0 = self.__palavra_aleatoria(self.__q[0])
      r1 = self.__palavra_aleatoria(self.__q[1])
      # Decide se inverte o spin
      dominio[i] = dominio[i] ^ (R2 | (R1 & r1) | r0)
    return dominio

  def magnetizacao(self):
    magnetizacoes = [0 for bit in xrange(0, self.bits)]
    i = 0
    for mask in [1 << x for x in xrange(0, self.bits)]:
      magnetizacoes[i] += float(sum([sum([(sitio & mask and 1)*2 - 1 for sitio in dominio]) for dominio in self.__dominios]))
      i += 1
    return sum(magnetizacoes)/(self.N*self.bits)

  def energia(self):
    redes = self.__redes()
    energia = 0
    for rede in redes:
      soma = 0
      for i in xrange(self.N):
        soma += rede[i]*(rede[(i-1)%self.N] + rede[(i-self.L)%self.N])
      energia -= self.J*float(soma)/self.N
    return energia/(self.N*self.bits)

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

  def __palavra_aleatoria(self, q):
    if random() < q:
      return 0
    else:
      return randint(0, maxint)

  def __redes(self):
    self.__recompoe_dominios()
    redes = [[] for bit in xrange(0, self.bits)]
    for inteiro in self.rede:
      i = 0
      for mask in [1 << x for x in xrange(0, self.bits)]:
        if inteiro & mask:
          redes[i].append(1)
        else:
          redes[i].append(-1)
        i += 1
    return redes