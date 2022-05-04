import numpy as np
import matplotlib.pyplot as plt
import os 
from funs import *
import scipy
from scipy import optimize
import time
from scipy import special

m = 50
# number of timesteps
t = 400
# Set magnetic field to 0
H = 0

# Create temperature array with different spacings in different regions
T = np.arange(0.1,2,0.2)
T = np.concatenate((T, np.arange(2,2.8,0.02)))
T = np.concatenate((T,np.arange(2.8,3.6,0.2)))

N = [100,150,200]
tc = []
tce = []
tx = []
txe = []

T_model = np.linspace(T[0],T[-1],1000)
plt.plot(T_model,modelC(T_model),ls = '--', c = 'k')
for n in N:
    
    [E,Ee,M,Me] = getData(n, t, m, H,T,0, False)

    mark = 200
    ET = np.array([ meanW(E[i][-mark:],Ee[i][-mark:]) for i in range(len(E))])
    ETe = np.array([ np.mean(x[-mark:]**2)**0.5 for x in Ee])
    
    MT = np.array([ meanW(M[i][-mark:],Me[i][-mark:]) for i in range(len(M))])
    MTe = np.array([ np.mean(x[-mark:]**2)**0.5 for x in Me])
    
    C = np.array([ n**2 *ETe[i]**2/T[i]**2 for i in range(len(ET))])
    Ce = np.array([ n**2/m * 2/T[i]**2 * np.sqrt( np.sum( (E[i][-mark:]-ET[i])**2 *Ee[i][-mark:]**2)) for i in range(len(ET))])
    
    X = np.array([ n**2 *MTe[i]**2/T[i] for i in range(len(MT))])
    Xe = np.array([ n**2/m * 1/T[i] * np.sqrt( np.sum( (M[i][-mark:]-MT[i])**2 *Me[i][-mark:]**2)) for i in range(len(MT))])
    

    ''' Plot magnetisation ''' #========================
    # plt.figure()
    #plt.plot(T_model,modelM(T_model),ls = '--', c = 'k')
    # plt.errorbar(T,MT, yerr = MTe, ls='',capsize = 4)
    # plt.ylabel("<m>")
    # plt.xlabel("T")
    # plt.legend(['model','n = 100', 'n = 150', 'n = 200'])
    # plt.savefig('mbad.png', dpi = 300)
    #===================================================
    
    ''' Plot energy ''' #===============================
    # plt.figure()
    # plt.plot(T_model,modelE(T_model),ls = '--', c = 'k')
    # plt.errorbar(T,ET, yerr = ETe,ls='',capsize = 4)
    # plt.ylabel("<E>")
    # plt.xlabel('T')
    # plt.legend(['model','data'])
    # plt.savefig('E(T).png', dpi = 300)
    #===================================================
    
    ''' plot Specific heat ''' #========================
    # plt.figure()
    # plt.plot(T_model,modelC(T_model),ls = '--', c = 'k')
    plt.errorbar(T,C,yerr = Ce,ls='',capsize = 4)
    plt.ylabel("<C>")
    plt.xlabel('T')
    plt.legend(['model','n = 100', 'n = 150', 'n = 200'])
    plt.savefig('Cbad.png', dpi = 300)
    #===================================================
    
    ''' plot Susceptibility ''' #=======================
    # plt.figure()
    # plt.errorbar(T,X,yerr = Xe,ls='',capsize = 4)
    # plt.ylim([min(X), max(X)])
    # plt.ylabel("X")
    # plt.xlabel('T')
    # plt.savefig('X100.png', dpi = 300)    
    #===================================================
    
    ''' plot autocorrelation ''' #======================
    # y = []
    # ye = []
    # Tt = []
    # b = 50
    # for i in range(5,len(T)):
    #     [x,xe] = autoCov(M[i],Me[i],200,t-b-1,b)
    #     if( np.any (xe == 0) ): 
    #         print('xe = 0 for ' + str(T[i]))
    #         continue
    #     popt,pcov = scipy.optimize.curve_fit(lambda x,a: np.exp(-x/a), np.arange(b), x,sigma = xe)
    #     y.append(popt[0])
    #     ye.append(np.sqrt(np.diag(pcov))[0] )
    #     Tt.append(T[i])
    
    # plt.errorbar(Tt,y,yerr=ye,ls='', marker = 'x')
    # plt.ylabel('decorelation time')
    # plt.xlabel('Temperature')
    # plt.legend(['10 x 10','100 x 100'])
    # plt.ylim([0,500])
    #===================================================
    
    ''' plot m(T) at some temperatures ''' #============
    # plt.figure()
    # leg = []
    # for i in [5,26,30,50]:
    #     plt.errorbar(np.arange(t), M[i], Me[i])
    #     leg.append('T = ' + str(round(T[i],2)))
    # plt.legend(leg)
    #===================================================
    
    ''' Fit approx C to data C ''' #====================
    # i = np.where( C == max(C) )[0][0]
    # s = i-10
    # e = i+10
    # C_fit = C[s:e]
    # Ce_fit = Ce[s:e]
    # T_fit = T[s:e]
    
    # f = lambda x,a,b,c: a/( (x-b)**2 +c)
    # pc,pce = scipy.optimize.curve_fit(f, T_fit,C_fit, sigma=Ce_fit)
    # pce = np.sqrt(np.diag(pce))
    # tc.append(pc[1])
    # tce.append(pce[1])
    
    # plt.figure()
    # plt.errorbar(T_fit, C_fit, yerr = Ce_fit, ls = '', marker = 'x',capsize=4) 
    # plt.plot(T_fit,f(T_fit,p[0],p[1],p[2]), 'k--' )
    
    #=================================================
    
    ''' Fit approx X to data X ''' #====================
    # i = np.where( X == max(X) )[0][0]
    # s = i-10
    # e = i+10
    # X_fit = X[s:e]
    # Xe_fit = Xe[s:e]
    # T_fit = T[s:e]
    # tx.append( T[ np.where( X == max(X))[0][0]] ) 
    # txe.append(0.02)
    # f = lambda x,a,b,c: a/( (x-b)**2 +c)
    # px,pxe = scipy.optimize.curve_fit(f, T_fit,X_fit, sigma=Xe_fit)
    # pxe = np.sqrt(np.diag(pxe))
    #tx.append(px[1])
    #txe.append(pxe[1])
    # plt.figure()
    # plt.errorbar(T_fit, X_fit, yerr = Xe_fit, ls = '', marker = 'x',capsize=4) 
    # plt.plot(T_fit,f(T_fit,px[0],px[1],px[2]), 'k--' )
    # plt.errorbar(T_fit, X_fit, yerr = Xe_fit, ls = '', marker = 'x',capsize=4) 
    # plt.plot(T_fit,f(T_fit,px[0],2.37,0.01), 'k--' )
    #=================================================
    

# g = lambda x,a,b,c :a/(1+b*x**(-1/c))
# p,pe = scipy.optimize.curve_fit(g, N,tc, sigma=tce)
# pe = np.sqrt(np.diag(pe))
# plt.figure()
# plt.errorbar(N, tc, yerr =tce, ls = '', capsize = 3, marker = 'x')
# plt.plot(N,g(N,p[0],p[1],p[2]), 'k--')
# plt.ylabel('Tc')
# plt.xlabel('N')
# plt.legend(['fit a + b N^(-1/c)','data points'])
# #plt.savefig('Tc_C.png', dpi = 300)
# print(p)
# print(pe)


# g = lambda x,a,b,c :a/(1+b*x**(-1/c))
# p,pe = scipy.optimize.curve_fit(g, N,tx, sigma=txe)
# pe = np.sqrt(np.diag(pe))
# plt.figure()
# plt.errorbar(N, tx, yerr =txe, ls = '', capsize = 3, marker = 'x')
# plt.plot(N,g(N,p[0],p[1],p[2]), 'k--')
# plt.ylabel('Tc')
# plt.xlabel('N')
# plt.legend(['fit a + b N^(-1/c)','data points'])
# #plt.savefig('Tc_X.png', dpi = 300)
# print(p)
# print(pe)

