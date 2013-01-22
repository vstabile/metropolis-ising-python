import time

S = 10

for l in xrange(10):
  L = 2**l # Sitios por lado
  N = L**2 # Total de sitios

  # Inicializa sistema
  ising = Ising(L, 1, 'down')
  ising.T = 2

  inicio = time.time()

  for s in xrange(S):
    ising.itera()
    ising.magnetizacao()

  fim = time.time()

  print (fim - inicio)/(S*N)
