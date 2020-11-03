import math
import matplotlib.pyplot as plt
import pandas as pd
from CoolProp.CoolProp import PropsSI
from tabulate import tabulate
     


df=pd.read_csv("GFG.csv")
values=df.value.to_list()
names=df.var_name.to_list()
values_converted=[]
for c in values:
    values_converted.append(eval(c))
data=dict(zip(names,values_converted))

pi=math.pi
Ib=data['Ib']
theta_deg=data['theta_deg']
theta=math.radians(theta_deg)
Ta=data['Ta']
Tsky=data['Tsky']
Ts=data['Ts']
Vw=data['Vw']
rho_w=data['rho_w']
k_air=data['k_air']
mu_w=data['mu_w']
sigma=data['sigma']
#Collector
Nm,Nr,W,module_length=data['Nm'],data['Nr'],data['W'],data['module_length']
L,f,Do,Di=data["L"],data['f'],data['Do'],data['Di']
Dco,Dci,epsilon_c,kc=data['Dco'],data['Dci'],data['epsilon_c'],data['kc']
kr,IF,Yr,Tg=data['kr'],data['IF'],data['Yr'],data['Tg']
alpha_a,nd=data['alpha_a'],data['nd']
#Solar Salt
kf,rho_f,Cpf,muf=data['kf'],data['rho_f'],data['Cpf'],data['muf']

#Assumptions
Ti=data['Ti']
mf_initial=data['mf_initial']
def ReHTFTest(mf,Ti):
    T=Ti-273
    rho_f=2090-0.636*T
    muf=(22.714-0.120*T+2.281*10**-4*T**2-1.474*10**-7*T**3)*10**-3
    Vf=(4*mf)/(rho_f*pi*Di**2)
    Re_f=(rho_f*Vf*Di)/muf
    return Re_f

def solarblock(mf,Ti):
    T=Ti-273
    rho_f=2090-0.636*T
    Cpf=1443+0.172*T
    muf=(22.714-0.120*T+2.281*10**-4*T**2-1.474*10**-7*T**3)*10**-3
    kf=0.443+1.9*10**-4*T
    
    Aa=(W-Dco)*L
    Ar=pi*Di*L
    QI=Ib*Aa*Nm*Nr*math.cos(theta)

    K_theta=math.cos(theta)-(2.859621*10**-5)*theta_deg**2-(5.25097*10**-4)*theta_deg
    Lsc=module_length*Nm
    EndLoss=1-(f*math.tan(theta))/Lsc
    Qa=Ib*K_theta*Aa*Nr*Nm*Tg*Yr*alpha_a*IF*nd*EndLoss

    Re_w=(rho_w*Vw*Dco)/mu_w
    if Re_w<1000:
        Nu_w=0.40+0.54*Re_w**0.52
    elif Re_w>=1000 and Re_w<50000:
        Nu_w=0.3*Re_w**0.6
    else:
        print("Reynold's Number for wind invalid")
    hw=(Nu_w*k_air)/Dco
    Tr=Ti+60
    epsilon_r=0.04795+0.0002331*(Tr-273)
    def QlossError(Tco):
        Qloss=pi*Dco*L*hw*(Tco-Ta)+epsilon_c*pi*Dco*L*sigma*(Tco**4-Tsky**4)
        Tci=((Qloss*math.log(Dco/Dci))/(2*pi*kc*L))+Tco
        test_Qloss=(pi*Do*L*sigma*(Tr**4-Tci**4))/((1/epsilon_r)+((1-epsilon_c)/epsilon_c)*(Do/Dci))
        error=Qloss-test_Qloss
        return error,Qloss,Tci

    Tco1=450
    err1,_,_=QlossError(Tco1)
    Tco2=310
    err2,_,_=QlossError(Tco2)
    Tco=(Tco1+((0-err1)*(Tco2-Tco1))/(err2-err1))

    _,Qloss,Tci=QlossError(Tco)
    UL=Qloss/(pi*Do*L*(Tr-Ta))

    Vf=(4*mf)/(rho_f*pi*Di**2)
    Re_f=(rho_f*Vf*Di)/muf
    if Re_f>4000:   #Turbulent
        Nu_f=(0.023*Re_f**0.8)*((muf*Cpf)/kf)**0.4
    elif Re_f<2300:  #Laminar
        Nu_f=3.66
    else:  #Transient
        print('Error: Transient flow not supported')

    hf=(Nu_f*kf)/Di
    F1=1/(UL*(1/UL+Do/(Di*hf)+(Do/(2*kr))*math.log(Do/Di)))
    FR=((mf*Cpf)/(UL*Ar))*(1-math.exp(-(UL*Ar*F1)/(mf*Cpf)))
    Qu_module=FR*(Qa/(Nr*Nm)-Ar*UL*(Ti-Ta))
    Qu_row=Qu_module*Nm
    del_T=Qu_row/(mf*Cpf)
    Te=Ti+del_T
    Qu=Qu_row*Nr
    Mf=mf*Nr
    XI=QI*(1-(4/3)*(Ta/Ts)+(1/3)*(Ta/Ts)**4)
    Xa=Qa*(1-Ta/Tr)
    Xu=Mf*Cpf*(Te-Ti-Ta*math.log(Te/Ti))
    package={'Aa':Aa,'Ar':Ar,'QI':QI,'K_theta':K_theta,'Lsc':Lsc,'EndLoss':EndLoss,
    'Qa':Qa,'hw':hw,'Tr':Tr,'epsilon_r':epsilon_r,'Tco':Tco,'Qloss':Qloss,'Tci':Tci,
    'UL':UL,'F1':F1,'FR':FR,'Qu_module':Qu_module,'Qu_row':Qu_row,'hf':hf,
    'Vf':Vf,'Re_f':Re_f,'Te':Te,'Qu':Qu,'Mf':Mf,'XI':XI,'Xa':Xa,'Xu':Xu}

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

    #Work and Heat calculations
    W_pump_1=(1-y-z)*m*w_pump_1*c
    W_pump_2=(1-y)*m*w_pump_2*c
    W_pump_3=m*w_pump_3*c
    Q_in=m*(h[7]-h[6])*c
    Q_out=m*(1-y-z)*(h[10]-h[1])*c
    W_turb=m*((h[7]-h[8])+(1-y)*(h[8]-h[9])+(1-y-z)*(h[9]-h[10]))*c
    W_net=(Q_in-Q_out)

    #Thermal efficiency
    n_thI=W_net/Q_in
    n_thI_overall=W_net/(QI*c)

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

    n_thII_overall=(X_TURB_in-X_PUMPI_in-X_PUMPII_in-X_PUMPIII_in)/(XI*c)


    X_dest_total=(X_HX_in-X_HX_out)+(X_FWH1_in-X_FWH1_out)+(X_FWH2_in-X_FWH2_out)+X_COND_LOSS #Total exergy destroyed

    mass_flow=[0,(1-y-z)*m,(1-y-z)*m,(1-y)*m,(1-y)*m,m,m,m,y*m,z*m,(1-y-z)*m]

    powerBlockOutputs={'w_pump_1':w_pump_1,'w_pump_2':w_pump_2,'w_pump_3':w_pump_3,
    'Q_in':Q_in,'Q_out':Q_out,'W_turb':W_turb,'W_net':W_net,'n_thI':n_thI,'X_HX_in':X_HX_in,
    'X_HX_out':X_HX_out,'X_PUMPI_in':X_PUMPI_in,'X_PUMPII_in':X_PUMPII_in,'X_PUMPIII_in':X_PUMPIII_in,
    'X_PUMPI_out':X_PUMPI_out,'X_PUMPII_out':X_PUMPII_out,'X_PUMPIII_out':X_PUMPIII_out,'X_FWH1_in':X_FWH1_in,
    'X_FWH1_out':X_FWH1_out,'X_FWH2_in':X_FWH2_in,'X_FWH2_out':X_FWH2_out,'X_TURB_in':X_TURB_in,
    'X_TURB_out':X_TURB_out,'X_COND_LOSS':X_COND_LOSS,'n_thII_HX':n_thII_HX,'n_thII_PUMP1':n_thII_PUMP1,
    'n_thII_PUMP2':n_thII_PUMP2,'n_thII_PUMP3':n_thII_PUMP3,'n_thII_TURB':n_thII_TURB,'n_thII_FWH1':n_thII_FWH1,
    'n_thII_FWH2':n_thII_FWH2,'X_dest_total':X_dest_total,'mass_flow':mass_flow,'n_thI_overall':n_thI_overall,
    'n_thII_overall':n_thII_overall}
    powerBlockOutputs.update(package)
    return powerBlockOutputs