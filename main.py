# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 17:15:09 2022

@author: Ahmad Aslan

"""

import numpy as np
from scipy.integrate import RK45, odeint
import matplotlib.pyplot as plt

import func



#==============================================================================
#Mit odeint:
    
def calculate_odeint(t0, tb, y0, multpl=1e2):
    #Löst die DGL mit dem odeint verfahren aus scipy zu einem gegebenen Anfangswert
    #Das label soll die Anfangswerte in dem Plot veranschaulichen
    
    #für den Flug wurde ein Zeitrahmen von maximal 20 sekunden geschätzt
    #die konstante Schrittweite beträgt hier 1e-2
    eta = np.linspace(t0, tb, int(tb*multpl))
    #Lösen der DGL...
    sol = odeint(func.raketen_gleichung_odeint, y0, eta)
        
    
    #Da wir nur die Werte von v und z bis zur Landung bewerten wollen sind dies 
    #alle Werte wo die Höhe z >= 0 ist. 
    count_zv = 0
    for  i in sol[:, 0]:
        if i >= 0 :
            count_zv += 1
        else:
            break
            
    #Suche den Hochpunkt von z, um diesen später zu markieren
    high_point_z = 0
    for i in range(len(sol[:, 0])):
        if sol[i, 0] < sol[i+1, 0]:
            high_point_z += 1
        else:
            break
    
    #Suche den Hochpunkt von v, um diesen später zu markieren
    high_point_v = 0
    for i in range(len(sol[:, 1])):
        if sol[i, 1] < sol[i+1, 1]:
            high_point_v += 1
        else:
            break
    
    #Der Wasserspiegel soll nur da betrachtet werden, wo sich auch wasser in Ihm
    #befindet
    count_x = 0
    for i in sol[:, 2]:
        if i >= 0:
            count_x += 1
            
        else:
            break
            
    z = [eta[0:count_zv], sol[0:count_zv, 0]]
    v = [eta[0:count_zv], sol[0:count_zv, 1]]
    x = [eta[0:count_x], sol[0:count_x, 2]]
    
    z_hp = [eta[high_point_z], sol[high_point_z, 0]]
    v_hp = [eta[high_point_v], sol[high_point_v, 1]]

    
    return z, v, x, z_hp, v_hp



#==============================================================================
#mit RK45
def calculate_RK45(t0, tb, y0, rtol=1e-10, atol=1e-10):
    #Löst die DGL mit dem RK45 Verfahren zu einem bestimmten Anfangswert und geg.
    #Anfangs- und Endzeitpunkt.verwendet wird hier der Algorithmus von scipy.
    
    #Die nachfolgenden ergebnisse werden in Listen gespeichert
    values_z = []
    values_v = []
    values_x = []
    
    #die Liste der Zeitangaben
    times = []
    
        
    sol = RK45(func.raketen_gleichung_RK45, t0, y0, tb, rtol=rtol, atol=atol)
    
    values_z.append(sol.y[0]) 
    values_v.append(sol.y[1])
    values_x.append(sol.y[2])
    times.append(sol.t)
        
    while True:
        #RK45 sucht sich jetzt dynamisch, passende Zeitschritte und löst die DGL
        #mit jedem 'step' wird der nächste Zeitschirtt berechnet
        sol.step()
        #nach jedem Step werden die etwaigen Ergebnisse der Liste hinzugefügt
        values_z.append(sol.y[0]) 
        values_v.append(sol.y[1])
        values_x.append(sol.y[2])
        times.append(sol.t) 
        
        #wenn der Status der Gleichung "finished" erreicht hat, ist die DGL
        #fertig gelöst und der Loop kann beendet werden
        if sol.status == "finished":
            print(sol.step_size)
            break

    #Da wir nur die Werte von v und z bis zur Landung bewerten wollen sind dies 
    #alle Werte wo die Höhe z >= 0 ist. 
    count_zv = 0
    for i in values_z:
        if i >= 0:
            count_zv += 1
            
    #Suche den Hochpunkt von z, um diesen später zu markieren
    high_point_z = 0
    for i in range(len(values_z)):
        if values_z[i] < values_z[i+1]:
            high_point_z += 1
        else: 
            break
    
    #Suche den Hochpunkt von v, um diesen später zu markieren
    high_point_v = 0
    for i in range(len(values_v)):
        if values_v[i] < values_v[i+1]:
            high_point_v += 1
        else: 
            break
    
    #Der Wasserspiegel soll nur da betrachtet werden, wo sich auch wasser in Ihm
    #befindet
    count_x = 0
    for i in values_x:
        if i >= 0:
            count_x += 1
        else:
            break
        
    #die Wasserspiegel werte von m in mm umrechnen
    for i in range(len(values_x)):
        values_x[i] *= 1000
    
    
    z = [times[0:count_zv], values_z[0:count_zv]]
    v = [times[0:count_zv], values_v[0:count_zv]]
    x = [times[0:count_x], values_x[0:count_x]]
    z_hp = [times[high_point_z], values_z[high_point_z]]
    v_hp = [times[high_point_v], values_v[high_point_v]]
    
    
    return z, v, x, z_hp, v_hp

#==============================================================================
#plotter methode
def plot(x, y, param, y0, proce, color="blue", high_point=False, xHP=None, yHP=None):
    
    label = "$V_{W0}=%s * V_{ges}[m^3]$ \n$z_0=%s[m]$ \n$v_0=%s[m/s]$ \n$x_0=%s[m]$"%(
        round(func.anteil_VW0,2),round(y0[0],3),round(y0[1],3),round(y0[2],3)
        )
    
    if param == "z" or param == "Z":
        plt.figure()
        plt.plot(x, y, label=label, c=color)
        plt.plot(x, np.zeros(len(x)), "--", c=color)
        plt.plot(x[-1], y[-1], "o", c="black", label=str(round(x[-1],2))+"[s]")
        plt.xlabel("Zeit in $s$")
        plt.ylabel("Höhe der Wasserrakete in $m$")
        plt.title("Höhenverlauf der Wasserrakete mit ${}$".format(proce))
        if high_point:
            plt.plot(xHP, yHP, "x", c="black", label="high point")
        plt.legend()
        
    elif param == "v" or param == "V":
        plt.figure()
        plt.plot(x, y, label=label, c=color)
        plt.plot(x, np.zeros(len(x))+y[-1], "--", c=color)
        plt.plot(x[-1], y[-1], "o", c="black", label=str(round(x[-1],2))+"[s]")
        plt.xlabel("Zeit in $s$")
        plt.ylabel("Geschw. der Wasserrakete in $m/s$")
        plt.title("Geschwindigkeitsverlauf der Wasserrakete mit ${}$".format(proce))
        if high_point:
            plt.plot(xHP, yHP, "x", c="black", label="high point")
        plt.legend()
        
    elif param == "x" or param == "X":
        plt.figure()
        plt.plot(x, y, label=label, c=color)
        plt.plot(x, np.zeros(len(x))+y[-1], "--", c=color)
        plt.plot(x[-1], y[-1], "o", c="black", label=str(round(x[-1],2))+"[s]")
        plt.ylabel("Höhe des Wasserspiegels in $m$")
        plt.title("Höhenverlauf des Wasserspiegels mit ${}$".format(proce))
        if high_point:
            plt.plot(xHP, yHP, "x", c="black", label="high point")
        plt.legend()
        
    else:
        print("Ungültiger Parameter!")
        
 
#==============================================================================
#Volumen verändern
def setWaterVol(faktor):
    func.anteil_VW0 = faktor
    func.V_W0 = faktor * func.V_ges
    func.V_L0 = func.V_ges - func.V_W0
    

#==============================================================================
#main code...
if __name__ == "__main__":
    #Validierung...
    #Was passiert mit einer Anfangseswchw. von v0=10m/s und einer Wasserfüllmenge von 0L?
    #mit einer Anfangshöhe von 35 m
    
    #Wasservolumen auf 0 setzen:
    setWaterVol(0)
    
    #Anfangswerte
    z0 = 35
    v0 = 10
    x0 = 0
    
    #Lösen der DGL mit y0=[0,10,0] ==>
    y0 = [z0,v0,x0]
    #anfangszeit
    t0 = 0 #[s]
    #endzeit (ungef abschätzen)
    tb = 4 #[s]
    
    
    #Lösen der DGL mit odeint
    z1, v1, x1, z_hp1, v_hp1 = calculate_odeint(t0, tb, y0)
    #plotten z, v, x mit odeint
    plot(z1[0], z1[1], "z", y0, "odeint", high_point=True, xHP=z_hp1[0], yHP=z_hp1[1])
    plot(v1[0], v1[1], "v", y0, "odeint", high_point=True, xHP=v_hp1[0], yHP=v_hp1[1], color="red")
    plot(x1[0], x1[1], "x", y0, "odeint", color="green")
    
    #RK45
    z2, v2, x2, z_hp2, v_hp2 = calculate_RK45(t0, tb, y0)
    #plotten z, v, x mit RK45
    plot(z2[0], z2[1], "z", y0, "RK45", high_point=True, xHP=z_hp2[0], yHP=z_hp2[1])
    plot(v2[0], v2[1], "v", y0, "RK45", high_point=True, xHP=v_hp2[0], yHP=v_hp2[1], color="red")
    plot(x2[0], x2[1], "x", y0, "RK45", color="green")
    
    
    
    
    #Was passiert mit einer anfangsgeschwindigkeit von 0 und einem Wasservolumen
    #von einem drittel der gesamtmenge
    
    #Wasservolumen auf 1/3 setzen
    setWaterVol(1/3)
    
    #Anfangswerte
    z0 = 0
    v0 = 0
    x0 = 1/3 * func.h_R
    
    y0 = [z0,v0,x0]
    #Anfangszeit:
    t0 = 0 #[s]
    #Endzeit (ungef abschätzen)
    tb = 7 #[s]
    
    #Lösen der DGL mit odeint
    z1, v1, x1, z_hp1, v_hp1 = calculate_odeint(t0, tb, y0)
    #plotten z, v, x mit odeint
    plot(z1[0], z1[1], "z", y0, "odeint", high_point=True, xHP=z_hp1[0], yHP=z_hp1[1])
    plot(v1[0], v1[1], "v", y0, "odeint", high_point=True, xHP=v_hp1[0], yHP=v_hp1[1], color="red")
    plot(x1[0], x1[1], "x", y0, "odeint", color="green")
    
    #RK45
    z2, v2, x2, z_hp2, v_hp2 = calculate_RK45(t0, tb, y0)
    #plotten z, v, x mit RK45
    plot(z2[0], z2[1], "z", y0, "RK45", high_point=True, xHP=z_hp2[0], yHP=z_hp2[1])
    plot(v2[0], v2[1], "v", y0, "RK45", high_point=True, xHP=v_hp2[0], yHP=v_hp2[1], color="red")
    plot(x2[0], x2[1], "x", y0, "RK45", color="green")
    
