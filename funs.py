import os
import numpy as np
import scipy 
from scipy import integrate

def getData(n,t,m,H,T,random,h,write_new = False):
    # if h == true, T = const and H varies
    #Check if n has already been simulated, if so load
    name = str(n) + '/' # name of directory to save data
    x = os.popen('[ -d "' + str(n) + '" ] && echo "1"').read()
   
    if( len(x) == 0 ):  os.system('mkdir ' + str(n))
   
    if (not write_new and len(x) != 0):
        E = np.loadtxt( name + 'E.csv', delimiter = ',')
        Ee =  np.loadtxt( name + 'Es.csv', delimiter = ',')
        M = np.loadtxt( name + 'M.csv', delimiter = ',')
        Me = np.loadtxt( name + 'Me.csv', delimiter = ',')
    # if not, generate data using base executable
    else:
        
        temps =''
        for x in T:
            temps += ' ' + str(x)
        
        if( h == False):
            os.system('./base ' + str(n) + ' ' + str(t) + ' ' + str(m) + ' ' + str(random) +' ' +str(H) + temps )
        else:
            os.system('./base_h ' + str(n) + ' ' + str(t) + ' ' + str(m) + ' ' + str(random) +' ' +str(H) + temps )
     
        Ms = np.loadtxt('Mags.txt',delimiter =',',dtype='str')
        Mse = np.loadtxt('MagsSTD.txt',delimiter =',',dtype='str')
        
        En = np.loadtxt('Energies.txt' ,delimiter =',',dtype='str')
        Ene = np.loadtxt('EnergiesSTD.txt' ,delimiter =',',dtype='str')
            
        if( h == False):
            E = np.array([ [ float(x) for x in En[i]] for i in range(len(T)) ])
            Ee = np.array([ [ float(x) for x in Ene[i]] for i in range(len(T)) ])
            
            M = np.array([ [ float(x) for x in Ms[i]] for i in range(len(T)) ])
            Me = np.array([ [ float(x) for x in Mse[i]] for i in range(len(T)) ])
        else:
            E = np.array([ float(x) for x in En])
            Ee = np.array([  float(x) for x in Ene ])
            
            M = np.array([  float(x) for x in Ms ])
            Me = np.array([  float(x) for x in Mse ])

        # Delete created text files
        os.system("rm *.txt")      
        # Save generated data
        np.savetxt( name + 'E.csv',E, delimiter = ',')
        np.savetxt( name + 'Es.csv',Ee, delimiter = ',')
        np.savetxt( name + 'M.csv',M, delimiter = ',')
        np.savetxt( name + 'Me.csv',Me, delimiter = ',')

    return [E,Ee,M,Me]

def autoCov(arr,err, start,end,i_max):
    # function to return autocorrelation and its errors
    # wit maximum translation of i_max, [start:end] used as reference
    x = np.zeros(i_max)
    xe = np.zeros(i_max)
    for i in range(i_max):
        a1 = arr[start : end]
        a2 = arr[start + i: end + i]
        e1 = err[start : end]
        e2 = err[start + i: end + i]
        x[i] = np.mean(a1*a2) - np.mean(a1)*np.mean( a2 ) 
        xe[i] = np.sum( (a2*e1/i_max - np.mean(a2)*e1/i_max)**2 ) + np.sum( (a1*e2/i_max - np.mean(a1)*e2/i_max)**2 )
        xe[i] = np.sqrt(xe[i])
    if (x[0] != 0): return [x/x[0],xe/x[0]]
    return [x,xe]

def meanW(data, error):
    # returns weighted average
    if( np.any( error == 0) ): return np.mean(data)
    return np.sum(data*(error**(-1)))/np.sum(error**(-1))

def modelM(T):
    # returns model magnetisation per site for array T
    x = np.zeros(len(T))
    for i in range(len(x)):
        if (1 - np.sinh(2/T[i])**(-4) > 0):
            x[i] = (1 - np.sinh(2/T[i])**(-4))**(1/8)
        else:
            x[i] = 0
    return x

def modelE(T):
    # returns model energy per site for array T
    x = np.zeros(len(T))
    for i in range(len(T)):
        k1 = 2* np.tanh(2/T[i]) /np.cosh(2/T[i])
        k2 = 2*(np.tanh(2/T[i]))**2 -1
        K,_ = scipy.integrate.quad(lambda x: ( 1 - (k1*np.sin(x))**2 )**(-1/2), 0, np.pi/2)
        x[i] = - np.tanh( 2/T[i])**-1 *( 1 + 2/np.pi*k2*K ) 
    return x

def modelC(T):
    # returns model heat capacity for array T
    x = np.zeros(len(T))
    for i in range(len(T)):
        k1 = 2* np.tanh(2/T[i]) /np.cosh(2/T[i])
        k2 = 2*(np.tanh(2/T[i]))**2 -1
        K,_ = scipy.integrate.quad(lambda x: ( 1 - (k1*np.sin(x))**2 )**(-1/2), 0, np.pi/2)
        E,_ = scipy.integrate.quad(lambda x: ( 1 - (k1*np.sin(x))**2 )**(1/2), 0, np.pi/2)
        x[i] = 2/np.pi *(T[i]*np.tanh(2/T[i]))**(-2) * (2*K -2*E -(1-k2)*(0.5*np.pi + k2*K))
    return x

def approxC(T,Tc = 2/np.log(1+2**0.5)):
    # returns approximate heat capacity for array T ( good around Tc )
    K1 = lambda x: np.log( 2**0.5/abs(1/x-1/Tc) )
    return 2/np.pi * np.log(np.tan(np.pi/8))**2 * (K1(T) - 1 - np.pi/4) 
    
def approxE(T,Tc = 2/np.log(1+2**0.5)):
    # returns approximate Energy for array T ( good around Tc )
    k2 = lambda x: 2*(np.tanh(2/x))**2 -1
    K1 = lambda x: np.log( 2**0.5/abs(1/x-1/Tc + 1e-12) )
    return -np.tanh(2*T**-1)**-1 *( 1 + 2/np.pi*k2(T)*K1(T) )