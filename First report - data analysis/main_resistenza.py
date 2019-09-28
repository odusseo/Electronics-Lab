#!/usr/bin/env python3
import numpy as np
import array as arr
import csv
import matplotlib.pyplot as plt
import pdb
from uncertainties import ufloat
import csv
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from functions import *

###################################
# Import data function
###################################

def import_data(file):
    with open(file) as csvfile:
        y=np.empty([15,2])
        csvfile.seek(0)
        reader=csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC)
        for i,row in enumerate(reader):
            if row:
                y[i]=row
        x=y[:,1]
        y=y[:,0]
        return x,y

def import_data2(file):
    with open(file) as csvfile:
        y=np.empty([4,2])
        csvfile.seek(0)
        reader=csv.reader(csvfile,quoting=csv.QUOTE_NONNUMERIC)
        for i,row in enumerate(reader):
            if row:
                y[i]=row
        x=y[:,1]
        y=y[:,0]
        return x,y
###################################
# Random function
###################################
def trow():
    return (np.random.rand()*0.46)

# Main
def main():
    y = np.array([16814.6158214537,	427.130169405182,	161.604697431015,	104.487160095613,	89.2648800814352,	72.5740172397344])
    x = np.array([1002, 42457, 118850, 197600, 239300, 297000])
    dy = np.array([124.110262042707,	121.000578508906,	121.002984633184,	121.000324986300,	121.000446871527,	121.000078245646])
    B, dB=linear_regression_B(x,y,dy)
    print('Linear regression only B: volt =',repr(ufloat(B,dB)))

    ###################################
    # Random voltage
    ###################################
    # out=np.empty(15,dtype=np.float)
    # for i in range(15):
    #     out[i] = trow()
    # print('out:', out)

    ###################################
    # Import data and adjustments
    ###################################
    curr, volt = import_data('bobina.csv')
    curr2, volt2 = import_data('valle_grande.csv')
    pesi, res = import_data2('pesi.csv')

    curr =curr*1e-3
    curr2 =curr2*1e-3


    # Uncertainties

    # # # FS grande
    # dx = (0.0001/np.sqrt(12))
    # dy = (1/np.sqrt(12))
    # # FS piccolo
    # dx = (0.00001/np.sqrt(12))
    # dy = (0.04/np.sqrt(12))

    # # # bobina
    dx = (0.001/np.sqrt(12))
    dy = (0.002/np.sqrt(12))

    dy_trasf = 0
    dy_tot = np.sqrt(dy**2+dy_trasf**2)
    dy_weights = 1/(dy_tot**2)

    ###################################
    # Data anlysis
    ###################################
    # Linear regression
    A1, B1, dA1, dB1 = linear_regression_AB(curr, volt, dy_weights)
    A2, B2, dA2, dB2 = linear_regression_AB(curr2, volt2, dy_weights)

    dy_trasf = B1*dx
    dy_tot = np.sqrt(dy**2+dy_trasf**2)
    dy_weights = 1/(dy_tot**2)

    A1, B1, dA1, dB1 = linear_regression_AB(curr, volt, dy_weights)

    print('Linear regression 1: volt =',repr(ufloat(A1,dA1)),'+', repr(ufloat(B1,dB1)),'* curr')
    print('Linear regression 2: volt =',repr(ufloat(A2,dA2)),'+', repr(ufloat(B2,dB2)),'* curr')

    # Linear regression
    chi2_1 = chi2(volt, dy_tot, A1+B1*curr)
    print('Chi2 1:', chi2_1)


    B1 = (np.sum( curr*volt )/np.sum( curr*curr ))
    inc_m = (dy_tot / np.sqrt(np.sum(curr*curr)))

    ##### Print debugger
    print('out:', B1, inc_m)
    chi2_1 = chi2(volt, dy_tot, B1*curr)
    print('Chi2 1:', chi2_1)


    ###################################
    # Plots
    ###################################
    MONO_FIGSIZE=(6,5)
    DOUBLE_FIGSIZE=(11.5,7.3)
    curr = curr*1000
    dx = dx*1000
    curr2 = curr2*1000


    ######### Plot 1
    fig1=plt.figure(figsize=DOUBLE_FIGSIZE)
    ax11=fig1.add_axes((.1,.38,.85,.55))
    ax11.errorbar(curr,volt,yerr=dy,xerr=dx,fmt='r.',label='Dati sperimentali')
    ax11.set_xlabel('$I_{ICE}\ [mA]$', labelpad = -8)
    ax11.set_ylabel('$\Delta V_{ICE}\ [V]$')
    ax11.grid()

    # grande
    # x1 = np.linspace(0.0, 5.0)
    # piccolo
    # x1 = np.linspace(0.0, 0.5)
    # bobina
    x1 = np.linspace(0.0, 50.0)

    y1 = B1*(x1/1000)
    plt.plot(x1, y1, '-b', label = 'Retta di regressione lineare')
    ax11_zoom=fig1.add_axes((.75,.47,.19,.19))
    ax11_zoom.errorbar(curr[8], volt[8], yerr=dy,xerr=dx, fmt = 'r.')
    ax11_zoom.set_xlabel('$[mA]$')
    ax11_zoom.set_ylabel('$[V]$')
    ax11_zoom.grid()
    mark_inset(ax11, ax11_zoom, loc1=1, loc2=3, fc="none", ec="0.5")
    # ax11.set_xticklabels([])

    ax11_res = fig1.add_axes((.1,.08,.85,.2))
    ax11_res.errorbar(curr, volt-(B1*(curr/1000)), yerr=dy, xerr=dx, fmt = 'r.')
    ax11_res.set_xlabel('$I_{ICE}\ [mA]$')
    ax11_res.set_ylabel('$\Delta V_{ICE}\ [V]$')
    ax11_res.plot(x1, 0*x1, '-b', label = 'Retta di regressione lineare 2')
    ax11_res.grid()
    plt.title('Grafico dei residui rispetto alla regressione', fontsize=15)

    fig1.suptitle('Tensione misurata in funzione della corrente misurata - monte 50 V',fontsize=18)
    legend1 = ax11.legend(loc='upper left', shadow=True, prop={'size': 15})
    legend1.get_frame().set_facecolor('#FFFFFF')
    fig1.savefig('fig8.png', transparent=False, dpi=180, )


    ######## Plot 2
    fig1=plt.figure(figsize=(10,5))
    ax11=fig1.add_axes((.09,.1,.86,.81))
    ax11.errorbar(curr,volt,yerr=dy,xerr=dx,fmt='r.',label='Amperometro a monte')
    ax11.errorbar(curr2,volt2,yerr=dy,xerr=dx,fmt='b.',label='Amperometro a valle')
    ax11.set_xlabel('$I_{ICE}\ [mA]$', labelpad = -8)
    ax11.set_ylabel('$\Delta V_{ICE}\ [V]$')
    ax11.grid()

    # grande
    x1 = np.linspace(0.0, 5.0)
    # piccolo
    # x1 = np.linspace(0.0, 0.5)

    y1 = A1 + B1*(x1/1000)
    y2 = A2 + B2*(x1/1000)
    plt.plot(x1, y1, '-r', label = 'Regressione lineare amp a monte')
    plt.plot(x1, y2, '-b', label = 'Retta di regressione amp a valle')
    ax11_zoom=fig1.add_axes((.74,.2,.19,.19))
    ax11_zoom.errorbar(curr[10], volt[10], yerr=dy,xerr=dx, fmt = 'r.')
    ax11_zoom.set_xlabel('$[mA]$')
    ax11_zoom.set_ylabel('$[V]$')
    ax11_zoom.grid()
    mark_inset(ax11, ax11_zoom, loc1=1, loc2=3, fc="none", ec="0.5")
    # ax11.set_xticklabels([])

    fig1.suptitle('Tensione misurata in funzione della corrente misurata - bobina',fontsize=18)
    legend1 = ax11.legend(loc='upper left', shadow=True, prop={'size': 13})
    legend1.get_frame().set_facecolor('#FFFFFF')
    fig1.savefig('fig5.png', transparent=False, dpi=180, )


    # ##PLOT media pesata
    # fig7=plt.figure(figsize=DOUBLE_FIGSIZE)
    # ax71=fig7.add_subplot(1,1,1)
    # ax71.errorbar(x,g,xerr=dx,yerr=dg,fmt='b.',label='Pendolo semplice')
    # ax71.errorbar(x,g_1,xerr=dx,yerr=dg_1,fmt='r.',label='Pendolo fisico')
    # ax71.axhline(y=g0_w,color='b',label='Media pesata pendolo semplice')
    # ax71.axhline(y=g0_1w,color='r',label='Media pesata pendolo fisico')
    # ax71.axhline(y=9.806,color='#000000',label='Tabulato')
    # ax71.axhspan(ymin=g0_w-dg0_w,ymax=g0_w+dg0_w,color='b',alpha=.2)
    # ax71.axhspan(ymin=g0_1w-dg0_1w,ymax=g0_1w+dg0_1w,color='r',alpha=.2)
    # ax71.set_title('g=g(l)')
    # ax71.set_ylabel('g ['+l_units+'*'+t_units+'^-2]')
    # ax71.set_xlabel('L ['+l_units+']')
    # legend71 = ax71.legend(loc='lower right', shadow=True)
    # #legend71.get_frame().set_facecolor('#00FF69')
    # #finishing
    # fig7.suptitle('Dipendenza di g da l',fontsize=16)
    # fig6.savefig('fig7.png', transparent=False, dpi=180, )

    plt.show()























if __name__ == '__main__':
    main()
