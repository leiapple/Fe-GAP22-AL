#------------------------------------------------------------------
# Plot the max. predicted error as a function of K.
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------
def plot_error(output_error_all, output_error_k):
    
    fig, (ax, bx) = plt.subplots(nrows=1, ncols=2, figsize=(15,7))
    plt.rcParams['font.size'] = '16'
    for inpfile in [output_error_all, output_error_k]:
        inp = open(inpfile,'r')
        lines = inp.readlines()
        alldata = []
        max_error = []
        for line in lines:
            line = line.split()
            singledata = []
            for everyele in line:
                ele = float(everyele)
                singledata.append(ele)
            merror = max(singledata)
            max_error.append(merror)
            alldata.append(singledata)
        # x = np.arange(len(max_error))/240 * 0.7 + 1.1
        if inpfile == output_error_all:
            x = np.arange(len(max_error))
            ax.plot(x,max_error)
            ax.scatter(x,max_error)
        else:
            x = 1.1 + np.arange(len(max_error))*0.01
            bx.plot(x,max_error)
            bx.scatter(x,max_error)

    ax.set_title("GAP uncertainty VS. Simulation Step", fontsize=16)
    ax.set_xlabel("Loading step", fontsize=14)
    ax.set_ylabel("GAP Predicted Error, [$eV$]", fontsize=14)
    ax.grid(color='lightgrey',ls='--',lw=0.5)
    ax.tick_params(axis='both', labelsize=14)
    plt.axhline(y=0.01,ls='--',c='grey')
    
    bx.set_title("GAP uncertainty VS. SIF", fontsize=16)
    bx.set_xlabel("$K_I$, [$Mp\sqrt{m}$]", fontsize=14)
    bx.set_ylabel("GAP Predicted Error, [$eV$]", fontsize=14)
    bx.grid(color='lightgrey',ls='--',lw=0.5)
    bx.tick_params(axis='both', labelsize=14) 
    plt.axhline(y=0.01,ls='--',c='grey')
    plt.savefig('./gap_error/GAP_predicted_error.png',dpi=300)
    
# Plot data
plot_error('./gap_error/output_error_all','./gap_error/output_error_k')