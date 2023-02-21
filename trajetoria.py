import math
import matplotlib.pyplot as plt
import numpy as np
import random  
import sympy as sp
import pandas as pd

### Propriedades Físicas
g0 = 9.81
R0 = 6378.3888

alturas = []
densidades = []
pressoes = []
gravidades = []

for altura in np.arange (0,3,0.001):
    densidade = 1.255*math.exp(-altura/8.882)
    pressao = ((1.0113*(10**5))*math.exp(-altura/7.502))/10**6
    g = g0 * (R0/(R0 + altura))**2

    gravidades.append(g) 
    alturas.append(altura)
    pressoes.append(pressao)
    densidades.append(densidade)

g_medio = sum(gravidades)/len(gravidades)
p_medio = sum(densidades)/len(densidades)

plt.figure(figsize = (20,5))
plt.subplot(1, 3, 1)
plt.plot(alturas, pressoes)
plt.title('Pressão Atmosférica')
plt.ylabel('Pressâo em MPa')
plt.xlabel('Altura em Km')
plt.grid()

plt.subplot(1, 3, 2)
plt.plot(alturas, densidades)
plt.title('Densidade do Ar')
plt.ylabel('densidade em Kg/m³')
plt.xlabel('Altura em Km')
plt.grid()

plt.subplot(1, 3, 3)
plt.plot(alturas, gravidades)
plt.title('Variação da Gravidade')
plt.ylabel('gravidade em m/s²')
plt.xlabel('Altura em Km')
plt.grid()
plt.show()

### Dados do Foguete
massa_foguete = 30 #Kg
massa_motor_vivo = 13.4480 #Kg
massa_motor_morto = 7.9436 #Kg
massa_perdida = massa_motor_vivo - massa_motor_morto
massa_final = massa_foguete - massa_perdida
Isp = 129.8
mp = 5.017
mt = (massa_foguete + massa_final) /2
t1 = 2.697
p = 1.255
Cd = 0.5
area_transversal_cilindrica = (math.pi * (0.16**2))/4

### Métodos de análise de trajetória

## Equação de Tsiolkovsky para um estágio

ve = Isp * g_medio
v = ve * np.log(massa_foguete/massa_final)

#Equação horária da velocidade
t1 = 2.697
a = v / t1
h = (a * t1**2)/2
t2 = v/g_medio
H = h + v*t2 - (g_medio*t2**2)/2

v_tempos = []
tempos =[]

for t in np.arange(0,t1,0.01):
    v_tempo = a * t
    v_tempos.append(v_tempo)
    tempos.append(t)
    t_ultimo = t
for t in np.arange(0.01,t2,0.01):
    v_tempo_queda = v_tempo - g_medio * t
    v_tempos.append(v_tempo_queda)
    tempos.append(t + t_ultimo)

alturas = []
altura = 0

for vi in v_tempos:
    inc_altura = vi * 0.01
    altura = altura + inc_altura
    alturas.append(altura)

tempos_1 = tempos [:]
v_tempos_1 = v_tempos [:]
alturas_1 = alturas [:]

print("Equação de Tsiolkovsky para um estágio: ")
print("Velocidade de exaustão: ",ve, " m/s ou ", ve/343, "Mach")
print("Velocidade máxima: ", v, " m/s ou ", v/343, "Mach")
print("Aceleração média: ", a, " m/s²")
print("Altura propulsionada: ", h, " m")
print("Altura de apogeu: ", H, " m")
print("Tempo de subida: ",t1+t2, " s")
print("\n")

## Calculo pelo impulso com massa constante

It = Isp * mp * g_medio
ve = Isp*g_medio
V = It / mt
a = V / t1
h = (a * t1**2)/2
t2 = V/g_medio
H = h + V*t2 - (g_medio*t2**2)/2

v_tempos = []
tempos =[]

for t in np.arange(0,t1,0.01):
    v_tempo = a * t
    v_tempos.append(v_tempo)
    tempos.append(t)
    t_ultimo = t
for t in np.arange(0.01,t2,0.01):
    v_tempo_queda = v_tempo - g_medio * t
    v_tempos.append(v_tempo_queda)
    tempos.append(t + t_ultimo)

alturas = []
altura = 0

for vi in v_tempos:
    inc_altura = vi * 0.01
    altura = altura + inc_altura
    alturas.append(altura)

tempos_2 = tempos [:]
v_tempos_2 = v_tempos [:]
alturas_2 = alturas [:]

print("Calculo pelo impulso com massa constante: ")
print(It)
print("Velocidade de exaustão: ",ve, " m/s ou ", ve/343, "Mach")
print("Velocidade máxima: ", v, " m/s ou ", v/343, "Mach")
print("Aceleração média: ", a, " m/s²")
print("Altura propulsionada: ", h, " m")
print("Altura de apogeu: ", H, " m")
print("Tempo de subida: ",t1+t2, " s")
print("\n")

## Cálculo por quantidade de movimento e variação de massa do foguete

R = 0
ve = Isp * g_medio


for massa in np.arange(0.1,massa_perdida,0.01):
    massa_final = massa_foguete - massa
    massa = mp
    It = Isp * mp * g_medio
    velocidade = (It - (massa_final * g_medio + R)) / massa_final
    R = Cd * 0.5 * p * area_transversal_cilindrica* (velocidade**2) #N

#velocidade = (mp * ve - massa_final * g_medio - R) / massa_final
#velocidade  = ve * np.log(massa_foguete/massa_final)

#Equação horária da velocidade
t1 = 2.697
a = velocidade / t1
h = (a * t1**2)/2
t2 = velocidade/g_medio
H = h + velocidade*t2 - (g_medio*t2**2)/2

v_tempos = []
tempos =[]

for t in np.arange(0,t1,0.01):
    v_tempo = a * t
    v_tempos.append(v_tempo)
    tempos.append(t)

for t in np.arange(0,t2,0.01):
    v_tempo_queda = v_tempo - g_medio * t
    v_tempos.append(v_tempo_queda)
    tempos.append(t + t1)

alturas = []
altura = 0

for i in np.arange(0, len(v_tempos)):
    inc_altura = v_tempos[i] * 0.01
    altura = altura + inc_altura
    alturas.append(altura)

tempos_3 = tempos [:]
v_tempos_3 = v_tempos [:]
alturas_3 = alturas [:]

print("Cálculo por quantidade de movimento e variação de massa do foguete: ")
print("Velocidade de exaustão: ",ve, " m/s ou ", ve/343, "Mach")
print("Velocidade máxima: ", velocidade, " m/s ou ", velocidade/343, "Mach")
print("Aceleração média: ", a, " m/s²")
print("Altura propulsionada: ", h, " m")
print("Altura de apogeu: ", H, " m")
print("Tempo de subida: ",t1+t2, " s")
print("A força de arrasto no corpo do foguete é: ", R, " N")
print("\n")

## Método 1

ve = 1273.338 # m/s
massa_perdida = -massa_motor_vivo + massa_motor_morto
massa_final = (massa_foguete - massa_perdida)
tempo_queima = 2.697 #s
t0 = 0 #s
R = - massa_perdida/tempo_queima
print(R)

#Tempo de queima
velocidade_queima = - g_medio * tempo_queima - ve * np.log(1-(R*tempo_queima/massa_foguete))
altura_queima = ve*tempo_queima - 0.5 * g_medio * tempo_queima**2  + ((massa_foguete - R*tempo_queima)/R) * ve * np.log(1-(R*tempo_queima/massa_foguete))

#Tempo de inércia
tempo_queima_ao_apogeu =  velocidade_queima/(g_medio)
velocidade_final = velocidade_queima - g_medio * (tempo_queima_ao_apogeu)
apogeu = altura_queima + velocidade_queima * tempo_queima_ao_apogeu - 0.5*g_medio * tempo_queima_ao_apogeu**2
tempo_subida = tempo_queima + tempo_queima_ao_apogeu

tempos = []
alturas = []
velocidades = []

for t in np.arange(0,tempo_subida, 0.01 ):
  if t <= tempo_queima:
    velocidade_queima = - g_medio * t - ve * np.log(1-(R*t/massa_foguete))
    altura_queima = ve*t - 0.5 * g_medio * t**2  + ((massa_foguete - R*t)/R) * ve * np.log(1-(R*t/massa_foguete))
    velocidades.append(velocidade_queima)
    alturas.append(altura_queima)
  else:
    velocidade_final = velocidade_queima - g_medio * (t - tempo_queima)
    apogeu = altura_queima + velocidade_queima * (t - tempo_queima) - 0.5*g_medio * (t - tempo_queima)**2
    velocidades.append(velocidade_final)
    alturas.append(apogeu)

  tempos.append(t)

print("Velocidade de exaustão: ",ve, " m/s ou ", ve/343, "Mach")
print("Velocidade de queima: ", velocidade_queima, " m/s ou ", velocidade_queima/343, "Mach")
print("Altura propulsionada: ", altura_queima, " m")
print("Tempo de quiema até o apogeu: ", tempo_queima_ao_apogeu, " s")
print("Velocidade no apogeu: ", velocidade_final, " m/s")
print("Altura de apogeu: ", apogeu, " m")
print("Tempo de subida: ",tempo_subida," s")

tempos_4 = tempos [:]
v_tempos_4 = velocidades [:]
alturas_4 = alturas [:]

### Gráficos

plt.figure(figsize = (15,5))

## Velocidades
plt.subplot(1, 2, 1)

plt.plot(tempos_1,v_tempos_1, label = 'Tsiolkovsky') ##1
plt.plot(tempos_2,v_tempos_2, label = 'Impulso') ##2
plt.plot(tempos_3,v_tempos_3, label = 'Qtd de Mov') ##3
plt.plot(tempos_4, v_tempos_4, label = 'Mtd 1') ## 4

plt.title('Variação da velocidade com o tempo')
plt.ylabel('Velocidade em m/s')
plt.xlabel('Tempo em s')
plt.grid()
plt.legend()

## Alturas
plt.subplot(1, 2, 2)

plt.plot(tempos_1,alturas_1, label = 'Tsiolkovsky') ##1
plt.plot(tempos_2,alturas_2, label = 'Impulso') ##2
plt.plot(tempos_3,alturas_3, label = 'Qtd de Mov') ##3
plt.plot(tempos_4, alturas_4, label = 'Mtd 1') ##4

plt.title('Variação da altura com o tempo')
plt.ylabel('Altura em m')
plt.xlabel('Tempo em s')
plt.grid()
plt.legend()

plt.show()