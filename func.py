# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:30:32 2022

@authors: Ahmad Aslan, Ali Akbey, Kamil Hammoud
"""

import numpy as np

#Konstanten:
pi = np.pi

rA1 = 1.35e-2 #m ==> Radius der Düsenöffnung
A1 = pi * rA1**2 #m^2 ==> Fläche der Düsenöffnung
rA0 = 4.5e-2 #m ==> Radius des Wassertanks
A0 = pi * rA0**2 #m^2 ==> Fläche des Wassertanks

V_ges = 0.002 #m^3 ==> Gesamtvoloumen der Wasserrakete
anteil_VW0 = 1/3 #Dimensionslos ==> Anteil an Wasser im gesamtvolumen
V_W0 = anteil_VW0 * V_ges #m^3 ==> Volumen des Wassers zum Zeitpunkt 0
V_L0 = V_ges - V_W0 #m^3 ==> Volumen der Luft zum Zeitpunkt 0
h_R = V_ges/A0 #m ==> Höhe der Wasserrakete


rho_W = 997 #kg/m^3 ==> Dichte des Wassers
rho_L = 1.2041 #kg/m^3 ==> Dichte des Wassers

P_L0 = 8e5 #Pa ==> Anfangsluftdruck der Wasserrakete
P_out = 101325 #Pa ==> Atmosphärischer Außendruck

m0 = 0.15 #kg ==> trockene Masse der Wasserrakete
g = 9.81 #m/s^2 ==> Erdbeschleunigung

C_w = 1/2 #Dimensionslos ==> Proportionalitätskonstante
zeta = 1 #Dimensionslos ==> Verlustziffer
kappa = 1.4 #Dimensionslos ==> Isentropen Exponent



#Differntialgleichung
def raketen_gleichung_odeint(y, xx):
    #Lösung der DGL mit odeint
    
    #Anfangswert entpacken
    z, v, x = y
    
    if x > 0:
        #Wenn Treibstoff vorhanden, dann...
        #V(x), w(x), m(x) und mu(x) evaluieren
        V_Wx = A0 * x
        V_Lx = V_ges - V_Wx
        w_x = np.sqrt( 2 * ( (P_L0 * ((V_L0/V_Lx)**kappa)) - P_out) / rho_W ) * zeta
        m_x = (m0 + (V_Wx * rho_W))
        mu_x = rho_W * A1 * w_x
    else:
        #wenn kein Treibstoff vorhandeb, dann w=0, mu=0, m=const, V_W=0
        #und daraus folgt: dz=0, dv=0, dx=0
        w_x = 0
        m_x = m0
        mu_x = 0
        
    
    #Lösung der dreidimensionalen Differentialgleichung #Siehe Mathematischer Anteil
    dz = v
    #mit luftwiderstand
    dv = ( ( - (m_x * g) - ( 0.5 * A0 * rho_L * v * abs(v) * C_w ) )/m_x ) + ( (mu_x/m_x) * w_x ) 
    dx =  -A1/A0 * w_x
    
    #in Liste (Vektor) verpacken
    dydt = [dz, dv, dx]
    
    #return der endgültigen Lösung der DGL
    return dydt

def raketen_gleichung_RK45(xx, y):
    #Lösung der DGL mit RK45
    
    #Anfangswert entpacken 
    z, v, x = y
    
    if x > 0:
        #wenn Treibstoff vorhanden, dann...
        #V(x), w(x), m(x) und mu(x) evaluieren
        V_Wx = A0 * x
        V_Lx = V_ges - V_Wx
        w_x = np.sqrt( 2 * ( (P_L0 * ((V_L0/V_Lx)**kappa)) - P_out) / rho_W ) * zeta
        m_x = (m0 + (V_Wx * rho_W))
        mu_x = rho_W * A1 * w_x
    else:
        
        #wenn kein Treibstoff vorhandeb, dann w=0, mu=0, m=const, V_W=0
        #und daraus folgt: dz=0, dv=0, dx=0
        w_x = 0
        m_x = m0
        mu_x = 0
        
    
    #Lösung der dreidimensionalen Differentialgleichung #Siehe Mathematischer Anteil
    dz = v
    #mit Luftwiderstand
    dv = ( ( - (m_x * g) - ( 0.5 * A0 * rho_L * v * abs(v) * C_w ) )/m_x ) + ( (mu_x/m_x) * w_x )
    dx =  -A1/A0 * w_x
    
    #in Liste (Vektor) verpacken
    dydt = [dz, dv, dx]
    
    #return der endgültigen Lösung der DGL
    return dydt