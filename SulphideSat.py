#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from ChemCal_v2 import CC
from ThermoEngine.thermoengine.thermoengine import phases

#%%

class SCSS:

    '''
    input example for enstaite chondrite

    EC_wt = pd.DataFrame(index = ['O','Si','Mg','Al','Ca','Fe','Ti','Mn','Cr','Na','K'],
                            data = [48.0,29.7,19.5,1.36,0.17,0.33,0.01,0.03,0.13,1.04,0.09],
                            columns = ['Wt. %'])
    
    '''

    def __init__(self,CC,inputDF,Tc,Pgpa,log_fo2):

    

    def computeActivities(self):

        

        return Am


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
        
        S_ppm  = 14695/Tk - 9.656 + 1.02*np.log(Tk) + B_prime + ((C*P)/Tk)*np.sum((Xm*Am)/Tk) + np.log(a_FeS) - np.log(a_FeO)
        return S_ppm 


# %%
