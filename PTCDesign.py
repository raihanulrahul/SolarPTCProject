import math
import matplotlib.pyplot as plt
pi=math.pi
Ib=800
theta_deg=23.314
theta=math.radians(theta_deg)
Ta=303
Tsky=289.12
Ts=5600
Vw=5
rho_w=1.164
k_air=0.02662
mu_w=18.6*10**-6
sigma=5.67*10**-8
#Collector
Nm,Nr,W,module_length=12,96,5.77,12.27
L,f,Do,Di=11.9,1.71,0.07,0.065
Dco,Dci,epsilon_c,kc=0.115,0.109,0.88,1.4
kr,IF,Yr,Tg=16,0.92,0.92,0.945
alpha_a,nd=0.94,0.98
#Solar Salt
kf,rho_f,Cpf,muf=0.55,1899,1495,0.00326

#Assumptions
Ti=563
mf_initial=0.1

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
def coverOuterTemp(mf):
    Vf=(4*mf)/(rho_f*pi*Di**2)
    Re_f=(rho_f*Vf*Di)/muf
    if Re_f>4000:
        Nu_f=(0.023*Re_f**0.8)*((muf*Cpf)/kf)**0.4
       # print("Turbulent")
    elif Re_f<2300:
        Nu_f=3.66
        #print("Laminar")
    else:
       # print("############Transient#########")
        return None,Re_f,None,None,None,None,None,None

    hf=(Nu_f*kf)/Di
    F1=1/(UL*(1/UL+Do/(Di*hf)+(Do/(2*kr))*math.log(Do/Di)))
    FR=((mf*Cpf)/(UL*Ar))*(1-math.exp(-(UL*Ar*F1)/(mf*Cpf)))
    Qu_module=FR*(Qa/(Nr*Nm)-Ar*UL*(Ti-Ta))
    Qu_row=Qu_module*Nm
    del_T=Qu_row/(mf*Cpf)
    return del_T,Re_f,F1,FR,Qu_module,Qu_row,hf,Vf
#Selecting mass flowrate
def showGraph():
    mf=mf_initial
    M_turb=[]
    M_lam=[]
    M_tran=[]   
    T_turb=[]
    T_lam=[]
    plt.figure(dpi=210)
    for _ in range(0,200):
        del_T,Re_f,F1,FR,Qu_module,Qu_row,hf,Vf=coverOuterTemp(mf)
        if Re_f>4000:
            M_turb.append(mf)
            T_turb.append(del_T)
        elif Re_f<2300:
            M_lam.append(mf)
            T_lam.append(del_T)
        else:
            M_tran.append(mf)
        mf=mf+0.01
    plt.plot(M_turb,T_turb,color='k',marker='.')
    plt.plot(M_lam,T_lam,color='k',marker='.')
    plt.axvspan(min(M_lam),max(M_lam),facecolor='silver')
    plt.axvspan(min(M_tran),max(M_tran),facecolor='darkgrey')
    plt.axvspan(min(M_turb),max(M_turb),facecolor='dimgrey',alpha=0.9)

    plt.text(0.18, 2000, "Laminar",fontsize=8)
    plt.text(0.45, 2000, "Transient",fontsize=8)
    plt.text(1.25, 2000, "Turbulent",fontsize=8)
    plt.text(0.45, 150, "Stability Range",fontsize=8)
    plt.axhspan(1,310,facecolor='mediumseagreen',alpha=0.8)
    #plt.plot(M_lam+M_tran+M_turb,[303]*len(M_lam+M_tran+M_turb),'k-')
    plt.xlabel('Solar Salt mass flowrate [kg/s]')
    plt.ylabel("Temperature gain [K]")
    plt.show()
    plt.savefig('fig1.png')

prompt1=input("Show del T vs mf graph? (y/n): ")
if prompt1=='y':
    showGraph()


mf=float(input("Enter corrected mass flowrate,mf [kg/s] (obtained using the graph)= "))

del_T,Re_f,F1,FR,Qu_module,Qu_row,hf,Vf=coverOuterTemp(mf)
Te=Ti+del_T
Qu=Qu_row*Nr
Mf=mf*Nr

print("""
Calculated design parameters:
Re_windflow= {}
Nu_windflow= {}
hw= {}

Calculated collector parameters:
Aa={}
K_theta= {}
Lsc= {}
EndLoss= {}
Ar= {}
Tr= {}
Tco= {}
Tci= {}
qloss= {}
UL= {}
F1= {}
FR= {}
epsilon_r= {}

HTF calculated properties:
mf= {}
vf= {}
hf= {}

Major derived quantities:
Mf= {}
Te= {}
Qu_module= {}
Qu_row= {}
Qu= {}
Qa= {}
QI= {}
""".format(Re_w,Nu_w,hw,Aa,K_theta,Lsc,EndLoss,Ar,Tr,Tco,Tci,Qloss,UL,F1,FR,epsilon_r,mf,
Vf,hf,Mf,Te,Qu_module/1000,Qu_row/1000,Qu/1000,Qa/1000,QI/1000)
)


XI=QI*(1-(4/3)*(Ta/Ts)+(1/3)*(Ta/Ts)**4)
Xa=Qa*(1-Ta/Tr)
Xu=Mf*Cpf*(Te-Ti-Ta*math.log(Te/Ti))

print("""
XI= {}
Xa= {}
Xu= {}

""".format(XI/1000,Xa/1000,Xu/1000))
