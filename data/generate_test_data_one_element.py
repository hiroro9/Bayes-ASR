import pandas as pd
import numpy as np
import random
from matplotlib import pyplot as plt
from numpy import linalg as LA

pi = np.pi
columns = ["time[h]"," CH01", "CH02", "CH03", "CH04", "CH05", "CH06", "CH07", "CH08", "CH09"]
en = pd.DataFrame(index=[],columns=columns)

"""=====  input parameters  ========"""
# file names #
file = "data-01"

# stress tensor [MPa]#
s11 = 35
s12 = 3
s13 = 5
s22 = 30
s23 = 6
s33 = 40

# pore pressure [MPa]#
p0 = 7

# Young's modulus [MPa]#
E = 75000

# Poisson's ratio #
u = 0.25

# viscosity of sample [Pa*s]#
Eta_v = 72*10**14
Eta_s = 49*10**14

"""==========================="""


"""==========  direction cosine vector  ============"""
n11 = 1
n21 = 0
n31 = 0

n12 = 0
n22 = 1
n32 = 0

n13 = 0
n23 = 0
n33 = 1

n14 = np.cos(pi/4)
n24 = np.sin(pi/4)
n34 = 0

n15 = np.cos(pi/4)
n25 = -np.sin(pi/4)
n35 = 0

n16 = 0
n26 = np.cos(pi/4)
n36 = np.sin(pi/4)

n17 = 0
n27 = np.cos(pi/4)
n37 = -np.sin(pi/4)

n18 = np.sin(pi/4)
n28 = 0
n38 = np.cos(pi/4)

n19 = -np.sin(pi/4)
n29 = 0
n39 = np.cos(pi/4)
"""==========================="""


"""======  calculated parameters (don't change this section!!!) ========="""
#bulk modulus [MPa]#
K = E/(3*(1-2*u))

#shear modulus [MPa]#
G = E/(2*(1+u))

#relaxation time [h]#
tv = Eta_v/((3*K*10**6)*3600)
ts = Eta_s/((2*G*10**6)*3600)

s = np.matrix([[s11,s12,s13],[s12,s22,s23],[s13,s23,s33]])
n = np.matrix([[n11,n12,n13,n14,n15,n16, n17, n18, n19],[n21,n22,n23,n24,n25,n26, n27, n28, n29],[n31,n32,n33,n34,n35,n36, n37, n38, n39]])
"""=================================================="""
output = file + ".csv"
figure_name = file + ".png"
txt_name = file + "_input_data.csv"

"""=====  model equations  ========"""
for min in range(0,2880,10):
  t = min/60
  Jav = (1-np.exp(-t/tv))/(3*K)
  Jas = (1-np.exp(-t/ts))/(2*G)
  
  es = (n.T@s@n - (1.0/3.0)*np.trace(s)*np.identity(9))*Jas
  ev = ((1.0/3.0)*np.trace(s) - p0)*Jav*np.identity(9)
  e = es + ev 
  e = e*10**6 + random.normalvariate(0,5)

  enn = np.array([[t,e[0,0],e[1,1],e[2,2],e[3,3],e[4,4],e[5,5],e[6,6],e[7,7],e[8,8]]])

  enp = pd.DataFrame(data=enn,columns=columns)
  en = en.append(enp)

en.to_csv(output,index=False)
"""================================"""

w, v = LA.eig(s)

"""=========  plot  ============"""
fig, ax = plt.subplots(1,1)
en.plot(x="time[h]",ax=ax)
#ax.ticklabel_format(style="sci", axis="y",scilimits=(0,0))
fig.savefig(figure_name)
"""================================"""

data_table = pd.DataFrame({
              "s11 [MPa]":[s11],
              "s12 [MPa]":[s12],
              "s13 [MPa]":[s13],
              "s22 [MPa]":[s22],
              "s23 [MPa]":[s23],
              "s33 [MPa]":[s33],
              "pore pressure [MPa]":[p0],
              "Young modulus [MPa]":[E],
              "poisson's ratio":[u],
              "bulk modulus K [MPa]":[K],
              "sheare modulus G [MPa]":[G],
              "Vol viscosity [Pa*s]":[Eta_v],
              "Shere viscosity [Pa*s]":[Eta_s],
              "Vol relax time [hour]":[tv],
              "Shere relax time [hour]":[ts],
              "sigma1":[w[0]],
              "sigma2":[w[1]],
              "sigma3":[w[2]],
              "n1":[v[0]],
              "n2":[v[1]],
              "n3":[v[2]]
              })
data_table.to_csv(txt_name,index = False)

print(w)
print(v)
