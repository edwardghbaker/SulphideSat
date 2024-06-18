#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ChemCal_v2 import CC
#from ThermoEngine.thermoengine.thermoengine import phases

#%%
EC_wt = pd.DataFrame(index = ['O','Si','Mg','Al','Ca','Fe','Ti','Mn','Cr','Na','K'],
                            data = [48.0,29.7,19.5,1.36,0.17,0.33,0.01,0.03,0.13,1.04,0.09],
                            columns = ['Wt. %'])



#%%

class SCSS:

    '''
    input example for enstaite chondrite

    EC_wt = pd.DataFrame(index = ['O','Si','Mg','Al','Ca','Fe','Ti','Mn','Cr','Na','K'],
                            data = [48.0,29.7,19.5,1.36,0.17,0.33,0.01,0.03,0.13,1.04,0.09],
                            columns = ['Wt. %'])
    
    '''

    def __init__(self,comp,Tc,Pgpa,log_fo2):
            
            self.comp = comp
            self.Tc = Tc
            self.Pgpa = Pgpa
            self.log_fo2 = log_fo2
    
            self.Xm = self.computeXm()
            self.Am = self.computeActivities()
            self.a_FeS = self.computeFeS()
            self.a_FeO = self.computeFeO()
            self.B_prime = self.computeBprime()
            self.C = self.computeC()
            self.S_ppm_sat = self.computeSCSS()


    def computeSCSS(self):
        '''
        This equation is taken directly from the Smythe paper.
        '''

        Tk = self.Tc + 273.15
        B_prime = self.B_prime
        C = self.C
        P = self.Pgpa*1e9
        Xm = self.Xm
        Am = self.Am
        a_FeS = self.a_FeS
        a_FeO = self.a_FeO
        
        S_ppm  = 14695/Tk - 9.656 + 1.02*np.log(Tk) + B_prime + ((C*P)/Tk)*np.sum(np.dot(Xm*Am)/Tk) + np.log(a_FeS) - np.log(a_FeO)
        return S_ppm 


# %%

def SCSS_oneill(Xca=0,Xmg=0,Xna=0,Xk=0,Xti=0,Xfe=0,Xal=0,fo2=10**(-11.73),fs2=10**(-2)):
     '''
     This function computes the sulphide capacity (ppm) for a given composition. As reported by O'neill, Hugh S. T. C., and John A. Mavrogenes. 2002. “The Sulfide Capacity and the Sulfur Content at Sulfide Saturation of Silicate Melts at 1400°C and 1 Bar.” Journal of Petrology 43 (6): 1049-87.
     
     this relationship is derived from 1atm experiments done at 1400C. The composition is in wt. %.
     '''
     
     A0 = -5.018
     Aca = 7.56
     Amg = 4.48
     Anak = 4.24
     Ati = 5.20
     Afe = 26.31
     Aal = 1.06
     Bfeti = 48.48


     return np.exp(A0+Xca*Aca+Xmg*Amg+(Xna+Xk)*Anak+Xti*Ati+Xfe*Afe+Xal*Aal+Xfe*Xti*Bfeti)*((fo2/fs2)**(-1/2))

def CS_approx_Oneill(wtFeO,fo2=1,fs2=1):
     Cs = 0.0003*(100-wtFeO)*np.exp(0.21*wtFeO)
     return Cs*((fo2/fs2)**(-1/2))


def SCSS_Liu(Tc,Pbar,xFeO,xTiO2,xCaO,xSiO2):
     
     lnxS = -1.76-(0.474*(10e4)/Tc)-0.021*Pbar+5.559*xFeO+2.565*xTiO2+2.709*xCaO-3.192*xSiO2
     xS = np.exp(lnxS)

     return xS



# %%
EC = CC(EC_wt)
EC_at = EC.get_at()
EC_ox_wt = EC.get_ox_wt()

SCSS_oneill(EC_at.loc['Ca'][0],EC_at.loc['Mg'][0],EC_at.loc['Na'][0],EC_at.loc['K'][0],EC_at.loc['Ti'][0],EC_at.loc['Fe'][0],EC_at.loc['Al'][0])

CS_approx_Oneill(EC_ox_wt.loc['FeO1.0'][0],10**(-11.73),10**(-2))

#%%
