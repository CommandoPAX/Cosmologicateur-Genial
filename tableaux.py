import numpy as np
import pandas as pd

snapshots = ["benchM", "G_ViVi","NG_F500", "NG_Fminus500","NEDE"]
labels = [r"$\Lambda$CDM", r"${\rm mixed~DM}$",  r"$f_{\rm NL} = -500$", r"$f_{\rm NL} = 500$ ", r"${\rm EDE}$"]

tab_minko = """
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



tab_corr = ""

tab_rarete = ""

##### Connectivite

tab_minko += r"$\langle \kappa \rangle $"

z= [15,12, 10, 8, 5,3,1,0.5,0.25,0]
indices_z = [5,6,8,9]

for i in range(len(snapshots)) :
    for j in [1,2,4]:
            z_k = indices_z[j-1]
                    
            connect_lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
            try : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_connect_fil.txt.npy")
            except : connect = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{z_k}_densite_smooth2_c0.1_connect_fil.txt.npy")
            moyenne = np.mean(connect)
            moyennes_lcdm = np.mean(connect_lcdm)
            tab_minko +=" &"+str(round((moyenne-moyennes_lcdm)/moyennes_lcdm *100, 1))+ r" \%"

tab_minko += r"\\"

##### Longueur


for i in range(len(snapshots)) :
    for j in [1,2,4]:
            z_k = indices_z[j-1]
            
            long_lcdm = np.load(f"/data100/fcastillo/RESULT/{snapshots[0]}/{k}_densite_smooth2_c0.1_len_fil.txt.npy")
            try : long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{k}_densite_smooth2_c0.1_len_fil.txt.npy")
            except :     long = np.load(f"/data100/fcastillo/RESULT/{snapshots[i]}/{j}_densite_smooth2_c0.1_len_fil.txt.npy")
            moyennes = np.mean(long)
            moyennes_lcdm = np.mean(long_lcdm)
            tab_minko +=" &"+str(round((moyenne-moyennes_lcdm)/moyennes_lcdm *100, 1))+ r" \%"

tab_minko += r"\\"