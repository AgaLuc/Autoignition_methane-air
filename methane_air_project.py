import sys
import numpy as np
from cantera import *
import matplotlib.pyplot as plt
import csv

gas = Solution('gri30.cti')

npoints = 10
Fipoints = 10

s = 0
nt = 100000
dt = 0.000001

Pmin = 101325
Pmax = 303975
Tmin = 1300
Tmax = 2500
Fimin = 0.4
Fimax = 3.1

Pi = np.zeros(npoints, 'd')
Ti = np.zeros(npoints, 'd')
Fi = np.zeros(Fipoints, 'd')
tim = np.zeros(nt, 'd')
temp = np.zeros(nt, 'd')
dtemp = np.zeros(nt - 1, 'd')
autoignition_csv = np.zeros(npoints ** 2 * Fipoints, 'd')
finaltemp_csv = np.zeros(npoints ** 2 * Fipoints, 'd')

autoignition_temp = np.zeros(npoints, 'd')
finaltemp_temp = np.zeros(npoints, 'd')
autoignition_p = np.zeros(npoints, 'd')
finaltemp_p = np.zeros(npoints, 'd')
autoignition_Fi = np.zeros(Fipoints, 'd')
finaltemp_Fi = np.zeros(Fipoints, 'd')

for j in range(npoints):
    Ti[j] = Tmin + (Tmax - Tmin) * j / (npoints - 1)

    for k in range(npoints):
        Pi[k] = Pmin + (Pmax - Pmin) * k / (npoints - 1)

        for l in range(Fipoints):
            Fi[l] = Fimin + (Fimax - Fimin) * l / (Fipoints - 1)
            nN2 = float(0.79 * 0.5 / (Fi[l] * 0.105042))
            nO2 = float(0.21 * 0.5 / (Fi[l] * 0.105042))
            X = 'CH4:0.5 N2:{0} O2:{1}'.format(str(nN2), str(nO2))
            gas.TPX = Ti[j], Pi[k], X
            r = IdealGasReactor(gas)
            sim = ReactorNet([r])

            time = 0.0

            for n in range(nt):
                time += dt
                sim.advance(time)
                tim[n] = time
                temp[n] = r.T

            Dtmax = [0, 0.0]
            for n in range(nt - 1):
                dtemp[n] = (temp[n + 1] - temp[n]) / dt
                if dtemp[n] > Dtmax[1]:
                    Dtmax[0] = n
                    Dtmax[1] = dtemp[n]
            autoignition = (tim[Dtmax[0]] + tim[Dtmax[0] + 1]) / 2.

            # print 'For T='+str(Ti[j])+'K P='+str(Pi[k])+'Pa and Fi='+str(Fi[l])+'Autoignition time='+str(autoignition)

            autoignition_csv[s] = autoignition * 1000
            finaltemp_csv[s] = temp[nt - 1]

            # print 'finaltemp='+str(finaltemp_csv[s])

            s += 1

            if Pi[k] == 101325 and Fi[l] == 1:
                print 'For T=' + str(Ti[j]) + 'K P=' + str(Pi[k]) + 'Pa and Fi=' + str(Fi[l]) + ',autoignition time=' + str(autoignition) + '(s)'
                print 'finaltemp=' + str(finaltemp_csv[s - 1])
                finaltemp_temp[j] = temp[nt - 1]
                autoignition_temp[j] = autoignition * 1000

            if Ti[j] == 1300 and Fi[l] == 1:
                print'For T=' + str(Ti[j]) + 'K P=' + str(Pi[k]) + 'Pa and Fi=' + str(
                    Fi[l]) + ',autoignition time=' + str(autoignition) + '(s)'
                print 'finaltemp=' + str(finaltemp_csv[s - 1])
                finaltemp_p[k] = temp[nt - 1]
                autoignition_p[k] = autoignition * 1000

            if Pi[k] == 101325.0 and Ti[j] == 1300.0:
                print 'For T=' + str(Ti[j]) + 'K P=' + str(Pi[k]) + 'Pa and Fi=' + str(
                        Fi[l]) + ',autoignition time=' + str(autoignition) + '(s)'
                print 'finaltemp=' + str(finaltemp_csv[s - 1])
                finaltemp_Fi[l] = temp[nt - 1]
                autoignition_Fi[l] = autoignition * 1000

s = 0
csv_file = 'autoignition_methane_air.csv'
with open(csv_file, 'w')as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Initial temperature', 'Pressure', 'Fi', 'Autoignition time', 'Final temperature'])
    for j in range(npoints):
        writer.writerow([Ti[j]])
        for k in range(npoints):
            writer.writerow(['', Pi[k]])
            for l in range(npoints):
                writer.writerow(['', '', Fi[l], autoignition_csv[s], finaltemp_csv[s]])
                s += 1

plt.plot(Fi, autoignition_Fi, 'b^')
plt.xlabel('Fi', fontsize=14)
plt.ylabel('autoignition [ms]', fontsize=14)
plt.title('autoignition of $CH_{4}$ and air mixture, P=1atm, T=1300K', fontsize=18,
horizontalalignment='center')
plt.grid()
plt.savefig('autoignition_Fi.png', bbox_inches='tight')
plt.show()

plt.plot(Fi, finaltemp_Fi,'b^')
plt.xlabel('Fi', fontsize=14 )
plt.ylabel('final temp [K]', fontsize=14)
plt.title('autoignition of $CH_{4}$ and air mixture, P=1atm, T=1300K', fontsize=18, horizontalalignment='center')
plt.grid()
plt.savefig('finaltemp_Fi.png', bbox_inches='tight')
plt.show()

plt.plot(Ti, autoignition_temp,'b^')
plt.xlabel('temp [K]', fontsize=14)
plt.ylabel('autoignition [ms]', fontsize=14)
plt.title('autoignition of $CH_{4}$ and air mixture, P=1atm, $\Phi$=1', fontsize=18, horizontalalignment='center')
plt.grid()
plt.savefig('autoignition_initialtemperature.png',bbox_inches='tight')
plt.show()

plt.plot(Ti, finaltemp_temp,'b^')
plt.xlabel('temp [K]',fontsize=14)
plt.ylabel('final temp [K]', fontsize=14)
plt.title('autoignition of $CH_{4}$ and air mixture, P=1atm, $\Phi$=1', fontsize=18, horizontalalignment='center')
plt.grid()
plt.savefig('finaltemp_temp.png', bbox_inches='tight')
plt.show()

plt.plot(Pi/1000,autoignition_p,'b^')
plt.xlabel('pressure [kPa]', fontsize=14)
plt.ylabel('autoignition [ms]')
plt.title('autoignition of $CH_{4}$ and air mixture, T=1300K, $\Phi$=1', fontsize=18, horizontalalignment='center')
plt.grid()
plt.savefig('autoignition_pressure.png', bbox_inches='tight')
plt.show()

plt.plot(Pi/1000,finaltemp_p,'b^')
plt.xlabel('pressure [kPa]',fontsize=14)
plt.ylabel('final temp [K]',fontsize=14)
plt.title('autoignition of $CH_{4}$ and air mixture, T=1300K, $\Phi$=1', fontsize=18, horizontalalignment='center')
plt.grid()
plt.savefig('finaltemp_pressure.png',bbox_inches='tight')


