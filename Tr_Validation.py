import PTCRankine as PR
import matplotlib.pyplot as plt
color=['black']
mf=1.25
Ti=563
label_list=['Ti=563']
i=0

x=10
M_turb=[]
M_lam=[]
M_tran=[]   
T_turb=[]
Tx_turb=[]
Tx_lam=[]
T_lam=[]

fig,ax = plt.subplots()

for _ in range(0,20):
    Re_f=PR.ReHTFTest(mf,Ti)
    if Re_f>4000:
        M_turb.append(x)
        solar=PR.solarblock(mf,Ti,x)
        T_turb.append((solar['n_thI']*100))
        Tx_turb.append((solar['Xa']-solar['Xu'])/10**3)
        #print(solar['Te'])
    elif Re_f<2300:
        M_lam.append(x)
        solar=PR.solarblock(mf,Ti,x)
        T_lam.append((solar['n_thI']*100))
        Tx_lam.append((solar['Xa']-solar['Xu'])/10**3)
        #print(solar['Te'])
    else:
        M_tran.append(x)
    x=x+10
ax.plot(M_turb,T_turb,color='navy',marker='s')
ax.plot(M_lam,T_lam,color='navy',marker='s')
ax2=ax.twinx()
ax2.plot(M_turb,Tx_turb,color='black',marker='D')

plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00))

ax.set_xlabel('Tr-Ti [K]',color='red',fontsize=14)
ax.set_ylabel("First Law Overall Efficiency [%]",color='navy',fontsize=14)
ax2.set_ylabel("Exergy Loss in Receiver [kW]",color='black',fontsize=14)
plt.show()
