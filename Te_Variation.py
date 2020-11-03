import PTCRankine as PR
import matplotlib.pyplot as plt
color=['gold','mediumvioletred','indigo','crimson','darkred']
Temp_list=[535,563,600,650,700]
label_list=['Ti=535K','Ti=563K','Ti=600K','Ti=650K','Ti=700K']
i=0
plt.figure(dpi=210)
for Ti in Temp_list:
    mf=0.1
    M_turb=[]
    M_lam=[]
    M_tran=[]   
    T_turb=[]
    T_lam=[]
    maxTe={}
    for _ in range(0,200):
        Re_f=PR.ReHTFTest(mf,Ti)
        if Re_f>4000:
            M_turb.append(mf)
            solar=PR.solarblock(mf,Ti)
            if solar['Te']<873:
                maxTe.update({solar['Te']:mf})
            T_turb.append(solar['Te'])
            #print(solar['Te'])
        elif Re_f<2300:
            M_lam.append(mf)
            solar=PR.solarblock(mf,Ti)
            if solar['Te']<873:
                maxTe.update({solar['Te']:mf})
            T_lam.append(solar['Te'])
            #print(solar['Te'])
        else:
            M_tran.append(mf)
        mf=mf+0.01
    print('For Ti={} :Highest mf={}'.format(Ti,maxTe[max(maxTe)]))
    plt.plot(M_turb,T_turb,color=color[i],marker='.',label=label_list[i])
    plt.plot(M_lam,T_lam,color=color[i],marker='.')
    i=i+1
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00))
#plt.axvspan(min(M_lam),max(M_lam),facecolor='silver')
#plt.axvspan(min(M_tran),max(M_tran),facecolor='darkgrey')
#plt.axvspan(min(M_turb),max(M_turb),facecolor='dimgrey',alpha=0.9)

#plt.text(0.18, 2000, "Laminar", fontsize=12)
#plt.text(0.45, 2000, "Transient", fontsize=12)
#plt.text(1.25, 2000, "Turbulent", fontsize=12)
plt.text(0.45, 400, "Stability Range", fontsize=12)
plt.axhspan(0,873,facecolor='royalblue',alpha=0.5)
#plt.plot(M_lam+M_tran+M_turb,[303]*len(M_lam+M_tran+M_turb),'k-')
plt.xlabel('Solar Salt mass flowrate [kg/s]')
plt.ylabel("Outlet temperature [K]")
plt.show()
