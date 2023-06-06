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

def SCSS_oneill(Xca=0,Xmg=0,Xna=0,Xk=0,Xti=0,Xfe=0,Xal=0):
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

     return np.exp(A0+Xca*Aca+Xmg*Amg+(Xna+Xk)*Anak+Xti*Ati+Xfe*Afe+Xal*Aal+Xfe*Xti*Bfeti)
     

def SCSS_Liu(Tc,Pbar,xFeO,xTiO2,xCaO,xSiO2):
     
     lnxS = -1.76-(0.474*(10e4)/Tc)-0.021*Pbar+5.559*xFeO+2.565*xTiO2+2.709*xCaO-3.192*xSiO2
     xS = np.exp(lnxS)

     return xS



# %%
