import PTCRankine as PR
import matplotlib.pyplot as plt
color=['gold','mediumvioletred','indigo','crimson','darkred']
Temp_list=[535,563,600,650,700]
mf_list=[0.88,0.86,0.95,1.11,1.34]
label_list=['Ti=535K','Ti=563K','Ti=600K','Ti=650K','Ti=700K']
i=0
plt.figure(dpi=210)
for Ti in Temp_list:
    mf=mf_list[i]
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
            T_turb.append(solar['n_thII_overall']*100)
            #print(solar['Te'])
        elif Re_f<2300:
            M_lam.append(mf)
            solar=PR.solarblock(mf,Ti)
            T_lam.append(solar['n_thII_overall']*100)
            #print(solar['Te'])
        else:
            M_tran.append(mf)
        mf=mf+0.01
    plt.plot(M_turb,T_turb,color=color[i],marker='.',label=label_list[i])
    plt.plot(M_lam,T_lam,color=color[i],marker='.')
    i=i+1
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00))

plt.xlabel('Solar Salt mass flowrate per row [kg/s]')
plt.ylabel("Second Law Overall Efficiency (%)")
plt.show()
