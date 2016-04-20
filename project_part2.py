'''It is desired to absorb 90% of the acetone in a gas containing 1.0 mol %
acetone in air in a countercurrent stage tower. The total inlet gas flow to the
tower is 30.0 kg mol/h, and the total inlet pure water flow to be used to
absorb the acetone is 90 kg mol H2O/h. The process is to operate isothermally
at 300 K and a total pressure of 101.3 kPa. The equilibrium
relation for the acetone in the gas-liquid is y = 2.53x. Determine the
number of theoretical stages required for this separation.'''

#Absorption of Acetone in a Countercurrent Stage Tower
import numpy as np
import matplotlib.pyplot as plt


#Variable Declaration
VN1 = 30.0       #Inlet Gas Flow Rate, kg mol/hr
L0 = 90.0        #Inlet Liquid Flow Rate, kg mol/hr
xA0 = 0.00       #Mole fraction of Acetone in entering Water
yAN1 = 0.01     #Mole fraction of Acetone in entering Air
AAAds = 90.00    #Percent of acetone absorbed
H = 2.53


#Calculations
        

AVi = yAN1*VN1 # moles of acetone entering in gas phase
Ai = (1-yAN1)*VN1 #moles of air entering in gas phase
AVo = (1-AAAds/100)*AVi #moles of acetone exiting in gas phase
ALo = AAAds*AVi/100 #moles of acetone exiting liquid phase
V1 = AVo + Ai  #total moles of exiting gas phase
yA1 = AVo/Ai #mole fraction of acetone in exiting gas phase
LN = L0 + ALo # total moles of liquid exiting the column
xAN = ALo/LN #mole fraction of acetone exiting the column

fig= plt.subplots()
y = np.arange(0,0.012,0.0005)
x = y/H

plt.plot(x,y)
plt.text(.0035, .007, 'Equilibrium Curve')
plt.text(.001, .009, 'Operating Line')
plt.plot(xAN,yAN1,'ro')
plt.annotate('$(x_{AN},y_{AN+1})$', xy=(xAN,yAN1), xytext=(xAN,yAN1))
plt.plot(xA0,yA1,'ro')
plt.annotate('$(x_{A0},y_{A1})$', xy=(xA0,yA1), xytext=(xA0,yA1+0.001))

m = (yAN1-yA1)/(xAN-xA0)

def ypos(xop):
    return m*xop + yA1

def xpos(yop):
    return yop/H

yy = m*x + yA1
plt.plot(x,yy,'r-')
plt.xlabel('Mole fraction of acetone in water, $x_A$')
plt.ylabel('Mole fraction of acetone in air, $y_A$')
plt.plot([xAN,xAN],[0.0,yAN1],'b--')
plt.plot([0.0,xAN],[yAN1,yAN1],'b--')
x1 = xA0
y1 = yA1
n = 0
while x1 <= xAN:
    x2 = xpos(y1)
    y2 = y1
    if x2 > xAN:
        plt.plot([x1,x2], [y1,y2], 'k-', lw=1)  #Draw Horizontal line to equilibrium curve
        dxt = x2-x1
        dx = xAN-x1
        dxbydxt = dx/dxt
        n = n + dxbydxt
        break
    plt.text(x2, y2-0.0007, str(n+1))
    plt.plot([x1,x2], [y1,y2], 'k-', lw=1)      #Draw Horizontal line to equilibrium curve
    plt.pause(0.01)     
    x1 = x2
    y1 = y2
    y2 = ypos(x1) 
    plt.plot([x1, x2], [y1, y2], 'k-', lw=1)    #Draw a vertical line to operating line
    plt.pause(0.01)
    n = n+1
    x1 = x2
    y1 = y2 
    
  
print "Acetone laden air rate from the absorber:",V1, "kgmol air/hr"
print "Acetone concentration in leaving stream", round(yA1,6)
print "Acetone + Water rate from the absorber:", LN, "kmol/hr"
print "Acetone concentrqtion in liquid leaving absorber:",round(xAN,5)
print "Amount of acetone etering with air", round(AVi,4),"kmol acetone/hr"
print "Amount of air entering ", round(Ai,4), " kgmol air/hr"
print "Amount of acetone leaving with treated air", round(AVo,4), "kgmol/hr"
print "Amount of acetone leaving with Water", round(ALo,4), "kgmol/hr"
print "Number of Stages:", round(n,1)