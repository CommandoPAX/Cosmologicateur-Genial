import numpy as np
import pandas as pd

snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
labels = [r"$\Lambda$CDM", r"${\rm mixed~DM}$",  r"$f_{\rm NL} = -500$", r"$f_{\rm NL} = 500$ ", r"${\rm EDE}$"]

tab_minko = r"""
\begin{table}[]
    \centering
    \begin{tabular}{|c|ccc|ccc|ccc|ccc|}
    \hline"""

for i in range(1,len(snapshots)):
    
    tab_minko += r"& \multicolumn{3}{c|}{"+labels[i]+"} "

tab_minko += r"""
 \\
       
         &$z = 3$ & $z= 1$ & $z = 0$&$z = 3$ & $z= 1$ & $z = 0$&$z = 3$ & $z= 1$ & $z = 0$&$z = 3$ & $z= 1$ & $z = 0$\\
         \hline 
"""

tab_minko +="\n"

tab_corr = ""

tab_rarete = ""

##### Connectivite

tab_minko += r"$\langle \kappa \rangle $"

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]

for i in range(1,len(snapshots)) :
    for j in [1,2,4]:
            z_k = indices_z[j-1]
                    
            connect_lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
            try : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
            except : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{z_k}_densite_smooth2_c0.1_connect_fil.txt.npy")
            moyenne = np.mean(connect)
            moyennes_lcdm = np.mean(connect_lcdm)
            tab_minko +=" &"+str(round((moyenne-moyennes_lcdm)/moyennes_lcdm *100, 1))+ r" \%"

tab_minko += r"\\" +"\n"

##### Longueur
tab_minko += r"$\langle l \rangle $"


for i in range(1,len(snapshots)) :
    for j in [1,2,4]:
            z_k = indices_z[j-1]
            
            long_lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{j}_densite_smooth2_c0.1_len_fil.txt.npy")
            try : long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_len_fil.txt.npy")
            except :     long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{z_k}_densite_smooth2_c0.1_len_fil.txt.npy")
            moyenne = np.mean(long)
            moyenne_lcdm = np.mean(long_lcdm)
            tab_minko +=" &"+str(round((moyenne-moyenne_lcdm)/moyenne_lcdm *100, 1))+ r" \%"

tab_minko += r"\\"

##### Minkowski

for p in range(4):
    tab_minko += rf"$v_{p}$"

    lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{i}.txt.npy")
    if i in [2,4]:lcdmzoom = np.load(f"/data100/fcastillo/RESULT/benchM/{i}_minkowski_zoom.txt.npy")

    X = np.linspace(-3,3,61)
    if i in [2,4]: 
        X1 = X[X<-1]
        X2 = np.linspace(-1,1,101)
        X3 = X[X>1]


        X = np.concatenate((X1,X2,X3))

    if p == 0 and i in [2,4] : X = X2


    for i in range(1,len(snapshots)) :
        for j in [1,2,4]:
        
                z_k = indices_z[j-1]
                lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{i}.txt.npy")

                try:
                    data = np.load(f"/home/fcastillo/minkowski/minkowski_{j}_{i}.txt.npy")
                except :
                    try : data= np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_minkowski.txt.npy")
                    except: data= np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{z_k}_minkowski.txt.npy")
                if i in [2,4]:
                    try :
                        datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{i}_minkowski_zoom.txt.npy")
                    except :
                        datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[j]}/{z_k}_minkowski_zoom.txt.npy")
                    data1 = data[:,:20]
                    data2 = data[:,41:]

                    data = np.concatenate([data1,datazoom, data2],axis=1)


                    lcdm = np.concatenate([lcdm[:,:20],lcdmzoom, lcdm[:,41:]],axis=1)

                    if p == 0 :
                        data = datazoom
                        lcdm = lcdmzoom

                max_delta = np.argmax(np.abs(data[p]-lcdm[p]))

                tab_minko +=" &"+str(round((data[p][max_delta]-lcdm[p][max_delta])/np.max(np.abs(lcdm[p])) *100, 1))+ r" \%"

    tab_minko += r"\\"

        



tab_minko+=r"""
\hline
\end{tabular}
\end{table}

"""

print(tab_minko)