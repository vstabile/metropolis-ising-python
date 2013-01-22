Otimizações Computacionais na Simulação Monte Carlo do Modelo de Ising
======

Implementações em Python do algoritmo de Monte Carlo usando condições de contorno helicoidais, computação paralela e código multispin.
    
Uso da classe Ising
------------

Inicializa um sistema com `512` sítios por lado da rede, constante de interação igual a `1` e estado inicial com spins aleatórios:

    ising = Ising(512, 1, 'random')

Define a temperatura na qual o sistema será simulado:

    ising.T = 2.2

Realiza N processos de Markov, tal que N é o número de sítios do sistema:

    ising.itera()

Mede a magnetização do sistem:

    ising.magnetizacao()

Mede a energia do sistema

    ising.energia()