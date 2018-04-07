import pandas as pd
import numpy as np
import cmath as math
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from windrose import WindroseAxes, WindAxes, plot_windrose
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
import seaborn as sns
from scipy.optimize import fsolve

def weib(v, k, A):
    return (k/A)*(v/A)**(k-1)*np.exp(-(v/A)**k)

 #odczyt z pliku
df = pd.read_excel('zad1_DANE.xls', sheetname='Arkusz1', decimal=",")


#stworzenie przegródek do histogramu
xbins = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
ybins = [0,30,60,90,120,150, 180, 210, 240, 270, 300, 330,360]
xnames = ['<1)','<1,2)','<2,3)','<3,4)','<4,5)','<5,6)','<6,7)','<7,8)','<8,9)','<9,10)','<10,11)','<11,12)','<12,13)','<13,14)','<14,15)','<15,16)','<16,17)','<17,18)','>18']
ynames = ['<0,30)','<30,60)','<60,90)','<90,120)','<120,150)','<150,180)','<180,210)','<210,240)','<240,270)','<270,300)','<300,330)','<330,360)',]

#stworzenie 2 struktur danych dla pomiarow na 2 wysokosciach
dfDir1 = df[['Spd1', 'SD1', 'Dir1']].copy()
dfDir2 = df[['Spd2', 'SD2', 'Dir2']].copy()


#liczba pomiarów
liczbaPomiarow = dfDir1.Spd1.count()


#predkosci srednie całoroczne (niezaleznie od kierunkow)
predkoscSrednia1 = dfDir1.Spd1.mean()
predkoscSrednia2 = dfDir2.Spd2.mean()

#obliczanie szorstkosci
wys1 = 70
wys2 = 30

szorstkoscLogarytmicznie = (((wys2)**(predkoscSrednia1/predkoscSrednia2))/wys1)**(1/((predkoscSrednia1/predkoscSrednia2)-1))
szortkoscWykladniczo = wys1/math.exp(math.log(wys1/wys2,predkoscSrednia1/predkoscSrednia2))

#wsp dla 70m
#metoda odchylenia standardowego 
dfDir1['Odchylania'] = (dfDir1.Spd1 - predkoscSrednia1)**2
sigma1 = (dfDir1['Odchylania'].sum()/liczbaPomiarow)**0.5
sigmaPrzezPrSrednia1 = (sigma1/predkoscSrednia1)**2
k1 = 0.975122 * sigmaPrzezPrSrednia1**-0.555089
A1 = (-0.667514 + 3.95422*k1 -3.70524*k1**2 + 1.88678*k1**3 -0.546985*k1**4 + 0.08490210*k1**5 -0.00548212*k1**6) * predkoscSrednia1


#wsp dla 30m
#metoda odchylenia standardowego
dfDir2['Odchylania'] = (dfDir2.Spd2 - predkoscSrednia2)**2
sigma2 = (dfDir2['Odchylania'].sum()/liczbaPomiarow)**0.5
sigmaPrzezPrSrednia2 = (sigma2/predkoscSrednia2)**2
k21 = 0.975122 * sigmaPrzezPrSrednia2**-0.555089
A21 = (-0.667514 + 3.95422*k21 -3.70524*k21**2 + 1.88678*k21**3 -0.546985*k21**4 + 0.08490210*k21**5 -0.00548212*k21**6) * predkoscSrednia2

##########################################################################################################

#tabela z histogramem prędkoci i energii dla 70m
histogramPredk1, histBin1 = np.histogram(dfDir1.Spd1,bins = 19, range = [0,19] )
histogramPredk1 = pd.DataFrame(histogramPredk1)
histogramPredk1.columns = ['ni']
histogramPredk1.index = xnames
histogramPredk1.rename_axis ('Przedzialy predkosci [m/s]')
histogramPredk1['ni/N'] = histogramPredk1['ni'] / liczbaPomiarow
histogramPredk1['Predkosci obl'] = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5, 17.5, 18.5]
histogramPredk1['Pi [W/m^2]'] = 0.5*1.2* histogramPredk1['Predkosci obl']*histogramPredk1['Predkosci obl']*histogramPredk1['Predkosci obl']*histogramPredk1['ni']
histogramPredk1['Ei [kWh/m^2/rok]'] = histogramPredk1['Pi [W/m^2]'] * 8760 / liczbaPomiarow/1000
energiaCalkowitaBrutto1 = np.sum(histogramPredk1['Ei [kWh/m^2/rok]'])
histogramPredk1['Ei/Eav'] = histogramPredk1['Ei [kWh/m^2/rok]']/energiaCalkowitaBrutto1

#tabela z histogramem predkosci i energii dla 30m
histogramPredk2, histBin2 = np.histogram(dfDir2.Spd2,bins = 19, range = [0,19] )
histogramPredk2 = pd.DataFrame(histogramPredk2)
histogramPredk2.columns = ['ni']
histogramPredk2.index = xnames
x = np.arange(0., 20., 0.2)
histogramPredk2.rename_axis ('Przedzialy predkosci [m/s]')
histogramPredk2['ni/N'] = histogramPredk2['ni'] / liczbaPomiarow
histogramPredk2['Predkosci obl'] = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5, 17.5, 18.5]
histogramPredk2['Pi [W/m^2]'] = 0.5*1.2* histogramPredk2['Predkosci obl']*histogramPredk2['Predkosci obl']*histogramPredk2['Predkosci obl']*histogramPredk2['ni']
histogramPredk2['Ei [kWh/m^2/rok]'] = histogramPredk2['Pi [W/m^2]'] * 8760 / liczbaPomiarow/1000
energiaCalkowitaBrutto2 = np.sum(histogramPredk2['Ei [kWh/m^2/rok]'])
histogramPredk2['Ei/Eav'] = histogramPredk2['Ei [kWh/m^2/rok]']/energiaCalkowitaBrutto2

#zapis histogramów do exela
writer = pd.ExcelWriter('Energia Wiatru.xlsx')
histogramPredk1.to_excel(writer,'Histogram 70m')
histogramPredk2.to_excel(writer,'Histogram 30m')
writer.save()

#histogramy prędkoci całorocznych (niezależnie od kierunków)
#histogram predkosci 70m
fig, ax = plt.subplots(figsize=(10, 10), dpi=160)
ax.hist(dfDir1.Spd1, bins = 19, range = [0, 19], normed = 1, edgecolor="k", linewidth=2)
plt.ylabel('Ni/N', size = 15)
histVal1, base1 = np.histogram(dfDir1.Spd1, bins = 19, range = [0, 19], normed = 1)
#base2 = base1[:-1] + (base1[1] - base1[0])/2
#f = UnivariateSpline(base2, histVal1)
#plt.plot(base2, f(base2), c = 'red')
#sns.distplot(dfDir1.Spd1, bins = 180, hist=False, color = 'red', norm_hist = 1)
x = np.arange(0., 20., 0.2)
plt.plot(x, weib(x, k1, A1))
plt.xticks(xbins)
#sns.kdeplot(dfDir1.Spd1, bw = 1, color = 'red', legend = 0)
cum1 = np.cumsum(histVal1)
ax2 = ax.twinx()
ax2.plot(base1[:-1], cum1, c='blue')
plt.ylim(ymin = 0)
plt.ylabel('Dystrybuanta', size = 15)
plt.xlabel('Predkosci [m/s]', size = 15)
plt.title('Histogram długookresowej predkosci wiatru na wys 70m', size = 20)
plt.savefig('Histogramd ługookresowej predkosci wiatru na wys 70m.png')
plt.show()

#histogram predkosci 30m
fig, ax = plt.subplots(figsize=(10, 10), dpi=160)
ax.hist(dfDir2.Spd2, bins = 19, range = [0, 19], normed = 1, edgecolor="k", linewidth=2)
plt.ylabel('Ni/N', size = 15)
histVal2, base2 = np.histogram(dfDir2.Spd2, bins = 19, range = [0, 19], normed = 1)
#base2 = base1[:-1] + (base1[1] - base1[0])/2
#f = UnivariateSpline(base2, histVal1)
#plt.plot(base2, f(base2), c = 'red')
#sns.distplot(dfDir1.Spd1, bins = 180, hist=False, color = 'red', norm_hist = 1)
x = np.arange(0., 20., 0.2)
plt.plot(x, weib(x, k21, A21))
plt.xticks(xbins)
#sns.kdeplot(dfDir2.Spd2, bw = 1, color = 'red', legend = 0)
#sns.kdeplot()
cum2 = np.cumsum(histVal2)
ax2 = ax.twinx()
ax2.plot(base1[:-1], cum2, c='blue')
plt.ylabel('Dystrybuanta', size = 15)
plt.xlabel('Predkosci [m/s]', size = 15)
plt.title('Histogram długookresowej predkosci wiatru na wys 30m', size = 20)
#plt.xticks()
plt.savefig('Histogramd ługookresowej predkosci wiatru na wys 30m.png')
plt.show()

#histogram energii 70m
energieDoHistogramu1 = pd.DataFrame (0.5*1.2 * 8760*dfDir1.Spd1*dfDir1.Spd1*dfDir1.Spd1/1000)
histogramPredk1[['Ei/Eav','ni/N' ]].plot(kind = 'bar', edgecolor="k", linewidth=1, width = 1 , figsize = (20, 20))
plt.title ('Histogram długoterminowego rozkładu gęstoci energii wys 70m', size = 40 )
plt.xticks(size = 30, rotation=45)
plt.yticks(size = 30)
plt.legend ( fontsize = 40)
plt.xlabel('Przedziały prędkoci v [m/s]', size = 40)
plt.ylabel('Względna liczba zdarzeń/Względny przyrost energii', size = 40)
plt.savefig('Histogram długoterminowego rozkładu gęstoci energii wys 70m.png')
plt.grid()
plt.show()

#histogram energii 30m
energieDoHistogramu2 = pd.DataFrame (0.5*1.2 * 8760*dfDir2.Spd2*dfDir2.Spd2*dfDir2.Spd2/1000)
histogramPredk2[['Ei/Eav','ni/N' ]].plot(kind = 'bar', edgecolor="k", linewidth=1, width = 1, figsize = (20, 20))
plt.title('Histogram długoterminowego rozkładu gęstoci energii wys 30m', size = 40)
plt.xticks(size = 30, rotation=45)
plt.yticks(size = 30)
plt.xlabel('Przedziały prędkoci v [m/s]', size = 40)
plt.ylabel('Względna liczba zdarzeń/Względny przyrost energii', size = 40)
plt.legend ( fontsize = 40 )
plt.savefig('Histogram długoterminowego rozkładu gęstoci energii wys 30m.png')
plt.grid()
plt.show()

#histogram 30
#fig, ax = plt.subplots(figsize=(10, 10), dpi=160)
#plt.hist(dfDir2.Spd2, bins = 18, range = [0, 18], normed = 1, edgecolor="k", linewidth=2 )
#plt.ylabel('Ni/N', size = 15)
#plt.xlabel('Predkosci [m/s]', size = 15)
#plt.title('Histogramd ługookresowej predkosci wiatru na wys 30m', size = 20)
#plt.xticks(xbins)
#plt.savefig('Histogramd ługookresowej predkosci wiatru na wys 30m.png')
#plt.show()

#histogram 70 i 30
fig, ax = plt.subplots(figsize=(10, 10), dpi=160)
plt.hist(dfDir1.Spd1, bins = 18, range = [0, 19], alpha = 0.7 , label = '70m', normed = 1, edgecolor="k", linewidth=2)
plt.hist(dfDir2.Spd2, bins = 18, range = [0, 19], alpha = 0.7 ,label = '30m', normed = 1, edgecolor="k", linewidth=2)
plt.legend(loc = 'upper right')
plt.ylabel('Ni/N', size = 15)
plt.xlabel('Predkosci [m/s]', size = 15)
plt.title('Histogramy długookresowej predkosci wiatru dla 70m oraz 30m', size = 20)
plt.xticks(xbins)
plt.savefig('Histogram predkosci 70m i 30m.png')
plt.show()



#podwójny histogram dla wys1
fig, ax = plt.subplots(figsize=(10, 10), dpi=160)
hist2, xbins, ybins = np.histogram2d(dfDir1['Dir1'],dfDir1['Spd1'], bins = [ybins, xbins], normed = 1 )
plt.hist2d(dfDir1['Spd1'],dfDir1['Dir1'], bins = [ybins, xbins])
plt.xlabel('Predkosci [m/s]', size = 15)
plt.ylabel('Kierunek [stopnie]', size = 15)
plt.title('Rozkład iloci pomiarów prędkoci w danym kierunku dla wys 70m', size = 20)
plt.xticks(ybins)
plt.yticks(xbins)
plt.savefig('Histogram 2d 70m.png')
plt.show()


#podwójny histogram dla wys2
fig, ax = plt.subplots(figsize=(10, 10), dpi=160)
#hist22, xbins2, ybins2 = np.histogram2d(dfDir2['Dir2'],dfDir2['Spd2'], bins = [ybins2, xbins2], normed = 1 )
plt.hist2d(dfDir2['Spd2'],dfDir2['Dir2'], bins = [ybins, xbins])
plt.xlabel('Predkosci [m/s]', size = 15)
plt.ylabel('Kierunek [stopnie]', size = 15)
plt.title('Rozkład iloci pomiarów prędkoci w danym kierunku dla wys 30m', size = 20)
plt.xticks(ybins)
plt.yticks(xbins)
plt.savefig('Histogram 2d 30m.png')
plt.show()


#tworzenie dodatkowych kolumn koniecznych do stworzenia rózy wiatrów - wys1
dfDir1['speed_x1'] = dfDir1.Spd1*(np.sin(dfDir1.Dir1 * np.pi / 180))
dfDir1['speed_y1'] = dfDir1.Spd1*(np.cos(dfDir1.Dir1 * np.pi / 180))
#wys2
dfDir2['speed_x2'] = dfDir2.Spd2*(np.sin(dfDir2.Dir2 * np.pi / 180))
dfDir2['speed_y2'] = dfDir2.Spd2*(np.cos(dfDir2.Dir2 * np.pi / 180))

#dziwny histogram potrzebny do róży
fig, ax = plt.subplots(figsize=(8, 8), dpi=80)
#tworzenie rózy wiatru wys1
x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
ax.set_aspect('equal')
_ = dfDir1.plot(kind='scatter', x='speed_x1', y='speed_y1', alpha=0.05, ax=ax)
Vw = 20
_ = ax.set_xlim([-Vw, Vw])
_ = ax.set_ylim([-Vw, Vw])
ax = WindroseAxes.from_ax()
ax.bar(dfDir1.Dir1.values, dfDir1.Spd1.values, bins=np.arange(0.0,19), cmap=cm.hot, lw=3, normed = 1)
ax.set_legend()
plt.title('Róża wiatru na wysokoci 70m')
plt.savefig('Róża wiatru 70m.png')


#dziwny histogram potrzebny do róży 2
fig2, ax = plt.subplots(figsize=(8, 8), dpi=80)
#tworzenie rózy wiatru wys1
x0, x2 = ax.get_xlim()
y0, y2 = ax.get_ylim()
ax.set_aspect('equal')
_ = dfDir2.plot(kind='scatter', x='speed_x2', y='speed_y2', alpha=0.05, ax=ax)
Vw = 20
_ = ax.set_xlim([-Vw, Vw])
_ = ax.set_ylim([-Vw, Vw])
ax = WindroseAxes.from_ax()
ax.bar(dfDir2.Dir2.values, dfDir2.Spd2.values, bins=np.arange(0.0,19), cmap=cm.hot, lw=3, normed = 1)
ax.set_legend()
plt.title('Róża wiatru na wysokoci 30m')
plt.savefig('Róża wiatru 30m.png')

