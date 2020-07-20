"""
This code is intended to solve the MC problem presented by Reinhard Hentschke in his 2019 article
Molecular Dynamics and Monte Carlo simulation of Lennard-Jones systems -a tutorial-
"""
import math as m
from random import random as ran

D = 2
dr = 0.1

class Metropolis1:

    r = []
    N = 16
    L = 8.0
    L2 = L*L

    def __init__(self):
        for i in range(self.N):
            self.r.append([self.L*ran(), self.L*ran()])

    def calculate_energy(self, n_particle, x0, y0):
        E = 0.0

        for n in range(len(self.r)):
            if n == n_particle:
                continue

            x_temp = self.r[n][0]
            y_temp = self.r[n][1]
            x_temp = x_temp - x0 - self.L*m.fabs((x_temp-x0)/self.L)
            y_temp = y_temp - y0 - self.L*m.fabs((y_temp-y0)/self.L)

            r2 = x_temp*x_temp + y_temp*y_temp

            if r2 > 0: E += 4*(m.pow(r2, -6.0) - m.pow(r2, -3.0))

        return E

    def metropolis_translation(self, beta):
        n = int(self.N*ran())
        x = self.r[n][0]
        y = self.r[n][1]

        drx = dr*(2*ran()-1)
        dry = dr*(2*ran()-1)

        x_new = x + drx
        y_new = y + dry

        if x_new < 0: x_new += self.L
        elif x_new > self.L: x_new -= self.L

        if y_new < 0: y_new += self.L
        if y_new > self.L: y_new -= self.L

        dE = self.calculate_energy(n, x_new, y_new) - self.calculate_energy(n, x, y)

        if dE <= 0:
            self.r[n][0] = x_new
            self.r[n][1] = y_new
        elif ran() < m.exp(-beta*dE):
            self.r[n][0] = x_new 
            self.r[n][1] = y_new

    def metropolis_transfer(self, beta, Box, factor):
        n = int(self.N*ran())
        x = self.r[n][0]
        y = self.r[n][1]

        drx = dr*(2*ran()-1)
        dry = dr*(2*ran()-1)

        x_new = x + drx
        y_new = y + dry

        if x_new < 0: x_new += Box.L
        elif x_new > Box.L: x_new -= Box.L

        if y_new < 0: y_new += Box.L
        if y_new > Box.L: y_new -= Box.L

        dE = Box.calculate_energy(n, x_new, y_new) - self.calculate_energy(n, x, y)

        if dE <= 0:
            Box.r.append([x_new, y_new])
            self.r.pop(n)
            self.N = len(self.r)
        elif ran() < factor*m.exp(-beta*dE):
            Box.r.append([x_new, y_new])
            self.r.pop(n)
            self.N = len(self.r)

    def volume_change_utilities(self, id):
        dV = self.L2*ran()*0.001
        V_new = self.L2 + id*dV
        change = (V_new/self.L2)**(1/3)
        temp = self.r

        for r in temp:
            r[0] *= change
            r[1] *= change

        return (self.calculate_energy(-1, 0,0), dV, change)

    def accept_volume(self, change, dV, id):
        for r in self.r:
            r[0] *= change
            r[1] *= change

        self.L2 += id*dV
        self.L = m.sqrt(self.L2)

class Metropolis2:
    
    r = []
    N = 16
    L = 8.0
    L2 = L*L

    def __init__(self):
        for i in range(self.N):
            self.r.append([self.L*ran(), self.L*ran()])

    def calculate_energy(self, n_particle, x0, y0):
        E = 0.0

        for n in range(len(self.r)):
            if n == n_particle:
                continue

            x_temp = self.r[n][0]
            y_temp = self.r[n][1]
            x_temp = x_temp - x0 - self.L*m.fabs((x_temp-x0)/self.L)
            y_temp = y_temp - y0 - self.L*m.fabs((y_temp-y0)/self.L)

            r2 = x_temp*x_temp + y_temp*y_temp

            if r2 > 0: E += 4*(m.pow(r2, -6.0) - m.pow(r2, -3.0))

        return E

    def metropolis_translation(self, beta):
        n = int(self.N*ran())
        x = self.r[n][0]
        y = self.r[n][1]

        drx = dr*(2*ran()-1)
        dry = dr*(2*ran()-1)

        x_new = x + drx
        y_new = y + dry

        if x_new < 0: x_new += self.L
        elif x_new > self.L: x_new -= self.L

        if y_new < 0: y_new += self.L
        if y_new > self.L: y_new -= self.L

        dE = self.calculate_energy(n, x_new, y_new) - self.calculate_energy(n, x, y)

        if dE <= 0:
            self.r[n][0] = x_new
            self.r[n][1] = y_new
        elif ran() < m.exp(-beta*dE):
            self.r[n][0] = x_new 
            self.r[n][1] = y_new

    def metropolis_transfer(self, beta, Box, factor):
        n = int(self.N*ran())
        x = self.r[n][0]
        y = self.r[n][1]

        drx = dr*(2*ran()-1)
        dry = dr*(2*ran()-1)

        x_new = x + drx
        y_new = y + dry

        if x_new < 0: x_new += Box.L
        elif x_new > Box.L: x_new -= Box.L

        if y_new < 0: y_new += Box.L
        if y_new > Box.L: y_new -= Box.L

        dE = Box.calculate_energy(n, x_new, y_new) - self.calculate_energy(n, x, y)

        if dE <= 0:
            Box.r.append([x_new, y_new])
            self.r.pop(n)
            self.N = len(self.r)
            print("CHANGE", end='\r', flush=True)
        elif ran() < factor*m.exp(-beta*dE):
            Box.r.append([x_new, y_new])
            self.r.pop(n)
            self.N = len(self.r)
            print("CHANGE", end='\r', flush=True)

    def volume_change_utilities(self, id):
        dV = self.L2*ran()*0.001
        V_new = self.L2 + id*dV
        change = (V_new/self.L2)**(1/3)
        temp = self.r

        for r in temp:
            r[0] *= change
            r[1] *= change

        return (self.calculate_energy(-1, 0,0), dV, change)

    def accept_volume(self, change, dV, id):
        for r in self.r:
            r[0] *= change
            r[1] *= change

        self.L2 += id*dV
        self.L = m.sqrt(self.L2)

if __name__ == "__main__":
    Box1 = Metropolis1()
    Box2 = Metropolis2()

    T = 0.1
    b = 1/T
    t_max = int(5e4)

    print(f"initial\nV1 = {Box1.L2}\nV2 = {Box2.L2}\nrho1: {len(Box1.r)/Box1.L2}\nrho2: {len(Box2.r)/Box2.L2}")

    for t in range(t_max):
        for i in range(len(Box1.r)): Box1.metropolis_translation(b)
        for i in range(len(Box2.r)): Box2.metropolis_translation(b)

        dE1, dV1, c1 = Box1.volume_change_utilities(1)
        dE2, dV2, c2 = Box2.volume_change_utilities(-1)
        dE = dE1 + dE2
        factor = (m.pow(1 + (dV1/Box1.L2), len(Box1.r)))*(m.pow(1-(dV2/Box2.L2), len(Box2.r)))

        if dE <= 0:
            Box1.accept_volume(c1, dV1, 1)
            Box2.accept_volume(c2, dV2, -1)
        elif ran() < factor*m.exp(-b*dE):
            Box1.accept_volume(c1, dV1, 1)
            Box2.accept_volume(c2, dV2, -1)

        factor1 = (Box2.L2*len(Box1.r))/((Box1.L2*len(Box2.r))+1)
        factor2 = (Box1.L2*len(Box2.r))/((Box2.L2*len(Box1.r))+1)

        Box1.metropolis_transfer(b, Box2, factor1)
        Box2.metropolis_transfer(b, Box1, factor2)

        print(t, end='\r', flush=True)

    print(f"last   \nV1 = {Box1.L2}\nV2 = {Box2.L2}\nN1 = {len(Box1.r)}\nN2 = {len(Box2.r)}")
    print(f"rho1: {len(Box1.r)/Box1.L2}\nrho2: {len(Box2.r)/Box2.L2}")


