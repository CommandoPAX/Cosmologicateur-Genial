import numpy as np
import pandas as pd

snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
labels = [r"$\Lambda$CDM", r"${\rm mixed~DM}$",  r"$f_{\rm NL} = -500$", r"$f_{\rm NL} = 500$ ", r"${\rm EDE}$"]


#snapshots = ["benchM", "NsPNG_EDE_F500","NsPNG_EDE_F1833", "NG_ViVi","NG_Fminus500_ViVi"]
#labels = [r"$\Lambda {\rm CDM}$", r"$f_{\rm NL} = -300~\&~ \rm{EDE}$",  r"$f_{\rm NL} = -1100~\&~{\rm EDE}$", r"$f_{\rm NL} = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL} = 500~\&~{\rm mixed~DM}$"]


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

tab_corr = r"""
\begin{table}[]
    \centering
    \begin{tabular}{|c|cc|cc|cc|cc|}
    \hline"""

for i in range(1,len(snapshots)):
    
    tab_corr += r"& \multicolumn{2}{c|}{"+labels[i]+"} "

tab_corr += r"""
 \\
       
         & $z= 1$ & $z = 0$  & $z= 1$ & $z = 0$  & $z= 1$ & $z = 0$  & $z= 1$ & $z = 0$\\
         \hline 
"""

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
tab_minko+="\n"

##### Minkowski

for p in range(4):
    tab_minko += rf"$v_{p}$"



    for i in range(1,len(snapshots)) :
        for j in [1,2,4]:
        
                z_k = indices_z[j-1]
                lcdm = np.load(f"/home/fcastillo/minkowski/minkowski_{0}_{j}.txt.npy")
                if j in [2,4]:lcdmzoom = np.load(f"/data100/fcastillo/RESULT/benchM/{j}_minkowski_zoom.txt.npy")

                try:
                    data = np.load(f"/home/fcastillo/minkowski/minkowski_{i}_{j}.txt.npy")
                except :
                    try : data= np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_minkowski.txt.npy")
                    except: data= np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{z_k}_minkowski.txt.npy")
                if j in [2,4]:
                    try :
                        datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_minkowski_zoom.txt.npy")
                    except :
                        datazoom = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{z_k}_minkowski_zoom.txt.npy")
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
    tab_minko+="\n"
        



tab_minko+=r"""
\hline
\end{tabular}
\end{table}

"""

print(tab_minko)

##### Fonctions de correlations


snapshots = ["benchM","NG_F500","G_m500","NG_F500_m500","NG_Fminus500","NG_Fminus500_m500", "G_ViVi","NG_ViVi" , "NG_Fminus500_ViVi","NEDE", "NsPNG_EDE_F500"]
labels = [r"$\Lambda$CDM", r"$f_{\rm NL} = -500$", "m = 500 eV", "WDM & fnl = -500", r"$f_{\rm NL} = 500$", "WDM & fnl = 500", r"$m_{\rm WDM} = 10~{\rm eV}, f_{\rm WDM} = 2~\%$", r"$f_{\rm NL} = -500~\&~{\rm mixed~DM}$", r"$f_{\rm NL} = 500~\&~{\rm mixed~DM}$", r"${\rm EDE}$", r"$f_{\rm NL} = -300~\&~{\rm EDE}$", r"$f_{\rm NL} = -1100~\&~{\rm EDE}$"]

indices_hdm = [7,8,10,11]
#indices_hdm = [0,2,6]
indices_hdm = [1,4,6,9]

##### Autocorrelations

R = 2
P = 5

for p in range(4):
    tab_corr += r"$\mathcal {" + (["P", "F", "W", "V"][p])*2 + r"}$"
    for i in indices_hdm :
        for j in [2,4]:
            z_k = indices_z[j-1]
           
            lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{j}_zeta_{p}_{p}_s{R}_P{P}.txt.npy")

            try : zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{i}_{j}_zeta_{p}_{p}_s{R}_P{P}.txt.npy")
            except : zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{i}_{z_k}_zeta_{p}_{p}_s{R}_P{P}.txt.npy")

            max_delta = np.argmax(np.abs(zeta-lcdm))

            tab_corr +=" &"+str(round(100*((zeta[max_delta]-lcdm[max_delta])/(1+np.max(np.abs(lcdm)))),1))+ r" \%"

    tab_corr+=r"\\"
    tab_corr+="\n"

for a in range(4):
    for b in range(a+1,4):
        tab_corr += r"$\mathcal {" + ["P", "F", "W", "V"][b]+["P", "F", "W", "V"][a] + r"}$"
        for i in indices_hdm :
            for j in [2,4]:
                z_k = indices_z[j-1]
            
                lcdm = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{0}_{j}_zeta_{b}_{a}_s{R}_P{P}.txt.npy")

                try : zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{i}_{j}_zeta_{b}_{a}_s{R}_P{P}.txt.npy")
                except : zeta = np.load(f"/data100/fcastillo/RESULT/extrema/snapshot_{i}_{z_k}_zeta_{b}_{a}_s{R}_P{P}.txt.npy")

                max_delta = np.argmax(np.abs(zeta-lcdm))

                tab_corr +=" &"+str(round(100*((zeta[max_delta]-lcdm[max_delta])/(np.max(1+np.abs(lcdm)))),1))+ r" \%"

        tab_corr+=r"\\"
        tab_corr+="\n"

tab_corr+=r"""
\hline
\end{tabular}
\end{table}

"""
print(tab_corr)