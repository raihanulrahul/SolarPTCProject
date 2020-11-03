from CoolProp.CoolProp import PropsSI
from tabulate import tabulate
import math

mf=1.25
Nr=96
Cpf=1495
Te=774.65
Ti=563
Mf=mf*Nr
c=10**-3

p=[0]*11
h=[0]*11
s=[0]*11
t=[0]*11
d=[0]*11
v=[0]*11
q=[0]*11

t[7]=Te

p[6],p[7]=9*10**6,9*10**6
p[4],p[5],p[8]=0.5044*p[6],0.5044*p[6],0.5044*p[6]
p[2],p[3],p[9]=0.2289*p[6],0.2289*p[6],0.2289*p[6]
p[10],p[1]=8*1000,8*1000

Ta=303


#State 1 [defined: p1,q1]
h[1]=PropsSI('H', 'P', p[1], 'Q', q[1], 'WATER')
s[1]=PropsSI('S', 'P', p[1], 'Q', q[1], 'WATER')
t[1]=PropsSI('T', 'P', p[1], 'Q', q[1], 'WATER')
d[1]=PropsSI('D', 'P', p[1], 'Q', q[1], 'WATER')
v[1]=1/d[1]


w_pump_1=v[1]*(p[2]-p[1]) #Pump 1 work output
h[2]=h[1]+w_pump_1
s[2]=s[1]

#Sate 2 [defined: h2, s2, p2]
t[2]=PropsSI('T', 'P', p[2], 'H', h[2], 'WATER')
d[2]=PropsSI('D', 'P', p[1], 'H', h[2], 'WATER')
v[2]=1/d[2]
q[2]=PropsSI('Q', 'P', p[2], 'H', h[2], 'WATER')


#State 3 [defined: p3, q-0]

h[3]=PropsSI('H', 'P', p[3], 'Q', q[3], 'WATER')
s[3]=PropsSI('S', 'P', p[3], 'Q', q[3], 'WATER')
t[3]=PropsSI('T', 'P', p[3], 'Q', q[3], 'WATER')
d[3]=PropsSI('D', 'P', p[3], 'Q', q[3], 'WATER')
v[3]=1/d[3]


w_pump_2=v[3]*(p[4]-p[3]) #Pump 2 work output
h[4]=h[3]+w_pump_2
s[4]=s[3]

#State 4 [defined: p4, s4, h4]
t[4]=PropsSI('T', 'P', p[4], 'H', h[4], 'WATER')
d[4]=PropsSI('D', 'P', p[4], 'H', h[4], 'WATER')
v[4]=1/d[4]
q[4]=PropsSI('Q', 'P', p[4], 'H', h[4], 'WATER')

#State 5 [defined: p5, q-0]
t[5]=PropsSI('T', 'P', p[5], 'Q', q[5], 'WATER')
h[5]=PropsSI('H', 'P', p[5], 'Q', q[5], 'WATER')
s[5]=PropsSI('S', 'P', p[5], 'Q', q[5], 'WATER')
d[5]=PropsSI('D', 'P', p[5], 'Q', q[5], 'WATER')
v[5]=1/d[5]

w_pump_3=v[5]*(p[6]-p[5]) #Pump 2 work output
h[6]=h[5]+w_pump_3
s[6]=s[5]

#State 6 [defined: p6,s6]
t[6]=PropsSI('T', 'P', p[6], 'S', s[6], 'WATER')
d[6]=PropsSI('D', 'P', p[6], 'S', s[6], 'WATER')
v[6]=1/d[6]
q[6]=PropsSI('Q', 'P', p[6], 'S', s[6], 'WATER')

#State 7 [defined: p7, T7]
h[7]=PropsSI('H', 'P', p[7], 'T', t[7], 'WATER')
s[7]=PropsSI('S', 'P', p[7], 'T', t[7], 'WATER')
d[7]=PropsSI('D', 'P', p[7], 'T', t[7], 'WATER')
v[7]=1/d[7]
q[7]=PropsSI('Q', 'P', p[7], 'T', t[7], 'WATER')

s[8]=s[7]
#State 8 [defined: p8, s8]
h[8]=PropsSI('H', 'P', p[8], 'S', s[8], 'WATER')
t[8]=PropsSI('T', 'P', p[8], 'S', s[8], 'WATER')
d[8]=PropsSI('D', 'P', p[8], 'S', s[8], 'WATER')
v[8]=1/d[8]
q[8]=PropsSI('Q', 'P', p[8], 'S', s[8], 'WATER')


s[9]=s[8]
#State 9 [defined: p9, s9]
h[9]=PropsSI('H', 'P', p[9], 'S', s[9], 'WATER')
t[9]=PropsSI('T', 'P', p[9], 'S', s[9], 'WATER')
d[9]=PropsSI('D', 'P', p[9], 'S', s[9], 'WATER')
v[9]=1/d[9]
q[9]=PropsSI('Q', 'P', p[9], 'S', s[9], 'WATER')

s[10]=s[9]
#State 10 [defined: p10, s10]
h[10]=PropsSI('H', 'P', p[10], 'S', s[10], 'WATER')
t[10]=PropsSI('T', 'P', p[10], 'S', s[10], 'WATER')
d[10]=PropsSI('D', 'P', p[10], 'S', s[10], 'WATER')
v[10]=1/d[10]
q[10]=PropsSI('Q', 'P', p[10], 'S', s[10], 'WATER')


y=(h[5]-h[4])/(h[8]-h[4])  #fraction of steam for regeneration
z=((h[3]-h[2])/(h[9]-h[2]))*(1-y)

#Solar Salt properties required for water mass flowrate calculation

m=(Mf*Cpf*(Te-Ti))/(h[7]-h[6])
print("\nm [kg/s]={}\ny= {}\nz= {} ".format(m,y,z))

#Work and Heat calculations
W_pump_1=(1-y-z)*m*w_pump_1*c
W_pump_2=(1-y)*m*w_pump_2*c
W_pump_3=m*w_pump_3*c
Q_in=m*(h[7]-h[6])*c
Q_out=m*(1-y-z)*(h[10]-h[1])*c
W_turb=m*((h[7]-h[8])+(1-y)*(h[8]-h[9])+(1-y-z)*(h[9]-h[10]))*c
W_net=(Q_in-Q_out)

FWH_1_in=y*m*((h[8]-h[5]))*c
FWH_1_out=(1-y-z)*m*((h[3]-h[2]))*c
FWH_2_in=y*m*((h[8]-h[5]))*c
FWH_2_out=(1-y)*m*((h[5]-h[4]))*c

print(FWH_1_in,FWH_1_out,FWH_2_in,FWH_2_out,FWH_1_in/FWH_1_out,FWH_2_in/FWH_2_out)


 #Thermal efficiency
n_thI=W_net/Q_in

#Exergy calculation

X_HX_in=Mf*Cpf*(Te-Ti-Ta*math.log(Te/Ti))*c
X_HX_out=m*((h[7]-h[6])-Ta*(s[7]-s[6]))*c
X_PUMPI_in,X_PUMPI_out=W_pump_1,W_pump_1
X_PUMPII_in,X_PUMPII_out=W_pump_2,W_pump_2
X_PUMPIII_in,X_PUMPIII_out=W_pump_3,W_pump_3
X_TURB_in,X_TURB_out=W_turb,W_turb
X_FWH1_in=z*m*((h[9]-h[3])-Ta*(s[9]-s[3]))*c
X_FWH1_out=(1-y-z)*m*((h[3]-h[2])-Ta*(s[3]-s[2]))*c
X_FWH2_in=y*m*((h[8]-h[5])-Ta*(s[8]-s[5]))*c
X_FWH2_out=(1-y)*m*((h[5]-h[4])-Ta*(s[5]-s[4]))*c
X_COND_LOSS=Ta*((1-y-z)*m*(s[1]-s[10])*c+Q_out/Ta)

n_thII_HX=X_HX_out/X_HX_in
n_thII_PUMP1=X_PUMPI_out/X_PUMPI_in
n_thII_PUMP2=X_PUMPII_out/X_PUMPII_in
n_thII_PUMP3=X_PUMPIII_out/X_PUMPIII_in
n_thII_TURB=X_TURB_out/X_TURB_in
n_thII_FWH1=X_FWH1_out/X_FWH1_in
n_thII_FWH2=X_FWH2_out/X_FWH2_in


X_dest_total=(X_HX_in-X_HX_out)+(X_FWH1_in-X_FWH1_out)+(X_FWH2_in-X_FWH2_out)+X_COND_LOSS #Total exergy destroyed

mass_flow=[0,(1-y-z)*m,(1-y-z)*m,(1-y)*m,(1-y)*m,m,m,m,y*m,z*m,(1-y-z)*m]

state_num=list(range(0,11))

## Properties at different states
print("Properties at different states:\n")
State_val= [(state_num[x],p[x]/1000,t[x],h[x]/1000,s[x]/1000,v[x],mass_flow[x],q[x])
for x in range(1, 11)]

print(tabulate(State_val, headers=["State","Pressure [kPa]",
"Temperature [K]","Enthalpy [kJ/kg]",
"Entropy[kJ/K]","Specific Volume [m3/kg]","Mass Flowrate [kg/s]",
"Quality"],tablefmt="plain"),"\n\n\n")


#First Law Analysis
print("First Law Analysis:\n")


print(tabulate([["Pump I Work Input [kW]",W_pump_1],["Pump II Work Input [kW]",W_pump_2],
["Pump III Work Input [kW]",W_pump_3],
["Turbine Work Output [kW]",W_turb],
["Net Work Output [kW]",W_net],["HX Heat Input [kW]",Q_in],
["Condenser Heat Output [kW]",Q_out],
["Thermal Efficiency [%]",n_thI]],floatfmt=".3f"),"\n\n\n")

#Second Law Analysis
print("Second Law Analysis:\n")

print("Exergetic Power Input and Output:\n")

print(tabulate([["Heat Exchanger Input",X_HX_in],["Heat Exchanger Output",X_HX_out],["Heat Exchanger Loss",X_HX_in-X_HX_out],
["Pump 1 Input",X_PUMPI_in],["Pump 1 Output",X_PUMPI_out],["Pump 1 Loss",X_PUMPI_in-X_PUMPI_out],
["Pump 2 Input",X_PUMPII_in],["Pump 2 Output",X_PUMPII_out],["Pump 2 Loss",X_PUMPII_in-X_PUMPII_out],
["Pump 3 Input",X_PUMPIII_in],["Pump 3 Output",X_PUMPIII_out],["Pump 3 Loss",X_PUMPIII_in-X_PUMPIII_out],
["Turbine Input",X_TURB_in],["Turbine Output",X_TURB_out],["Turbine Loss",X_TURB_in-X_TURB_out],
["FWH 1 Input",X_FWH1_in],["FWH 1 Output",X_FWH1_out],["FWH 1 Loss",X_FWH1_in-X_FWH1_out],
["FWH 2 Input",X_FWH2_in],["FWH 2 Output",X_FWH2_out],["FWH 2 Loss",X_FWH2_in-X_FWH2_out],
["Condenser Exergetic Power Loss",X_COND_LOSS]],floatfmt=".3f"))

print("\n Exergetic Efficiencies\n")
print(tabulate([["HX",n_thII_HX],["PUMP 1",n_thII_PUMP1,],["PUMP 2",n_thII_PUMP2],
["PUMP 3",n_thII_PUMP3],["Turbine",n_thII_TURB],["FWH 1",n_thII_FWH1],["FWH 2",n_thII_FWH2]],floatfmt=".3f"))