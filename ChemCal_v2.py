# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 14:04:55 2022

@author: Ed
"""
#%% import packages
import mendeleev as m
import pandas as pd 
import numpy as np
from molmass import Formula
from string import digits

#%% define functions

class CC:
    
    '''
    class for converting to and from any of at%, wt% or ox_wt%.
    
    input is a dataframe with the chemicals as the index, eg 'Fe' or FeO (Currently only accepts FeO not Fe2O3). Must have only one input column and the header of that column must be 'Wt. %', 'Ox. wt. %' or 'At. %'. 
    
    once read the output may be selected with the appropriate function.   
    
    ****
    
    input example:
    EC_wt = pd.DataFrame(index = ['O','Si','Mg','Al','Ca','Fe','Ti','Mn','Cr','Na','K'],
                             data = [48.0,29.7,19.5,1.36,0.17,0.33,0.01,0.03,0.13,1.04,0.09],
                             columns = ['Wt. %'])
    
    OxSates must also be a dataframe with the element and the oxidation state 
    '''
    suggestedOxiStates = pd.DataFrame(data = {'Charges':[4,2]},index=['Ti','Mn'])
    
    def __init__(self,
                 input_DF,FeCharge=2,
                 imposeOxidation=True,OxStates=suggestedOxiStates):
        
        if type(input_DF) != pd.core.frame.DataFrame:
            print(f'Incorect input format, you input as {type(input_DF)}. see documentation.')
        
        self.input_DF = input_DF 
        self.FeCharge = FeCharge
        self.imposeOxidation = imposeOxidation
        self.OxStates = OxStates
        
        if input_DF.columns[0] == 'Wt. %':
            self.at_DF = CC.wt2at(self)
            self.wt_DF = input_DF
        
        elif input_DF.columns[0] == 'Ox. Wt. %':
            self.at_DF = CC.ox_wt2at(self)
            self.ox_wt_DF = input_DF
        
        elif input_DF.columns[0] == 'At. %':
            self.ox_wt_DF = CC.at2ox_wt(self)
            self.wt_DF = CC.at2wt(self)
            self.at_DF = input_DF
          
    def mid(abc):
        return abc[(len(abc)//2)]
    
    def normalise(array,value=100):
        array = np.array(array)
        total = np.sum(array)
        norm_arr = value*array/total
        return norm_arr
    
    def wt2at(self,index=None,weight=None):
        if index is None or weight is None:
            index = self.input_DF.index
            weight = self.input_DF['Wt. %']
        Atomic = []
        
        for i,j in zip(index,weight):
            Atomic.append(j/m.element(i).atomic_weight)
        normAtomic = CC.normalise(Atomic)
        at_DF = pd.DataFrame(data = {'At. %' : normAtomic},index=index)
        #print(at_DF)
        return at_DF
    
    def ox_wt2at(self,ox_wtIndex=None,ox_weight=None):
        if ox_wtIndex is None or ox_weight == None:
            ox_wtIndex = self.input_DF.index
            ox_weight = self.input_DF['Ox. Wt. %']

        wtIndex_coeffs = [i.rsplit('O')[0] for i in ox_wtIndex]

        wtIndex = np.array([i.rstrip(digits) for i in wtIndex_coeffs])
        conversion_coeffs = []
        for i,j in zip(ox_wtIndex,wtIndex):
            conversion_coeffs.append((Formula(j).mass/Formula(i).mass))
        
        conversion_coeffs = np.array(conversion_coeffs)
        weight = np.array(ox_weight)*conversion_coeffs
        atomic = CC.wt2at(self,index=wtIndex,weight=weight)
        charges = np.array([CC.mid(m.element(i).oxistates) for i in wtIndex])
        Charges = pd.DataFrame(data={'Oxidation State':charges},index=wtIndex)
        Charges.loc['Fe'] = self.FeCharge
        if self.imposeOxidation == True:
            for i in self.OxStates.index:
                try: Charges.loc[i] = self.OxStates.loc[i][0]
                except: next
        atomic_O = np.sum(Charges['Oxidation State']*atomic['At. %'])/2
        atomic = np.append(atomic,atomic_O)
        norm_atomic = CC.normalise(atomic)
        wtIndex = np.append(wtIndex,'O')
        at_DF = pd.DataFrame(data = {'At. %' : norm_atomic}, index = wtIndex)
        return at_DF
        
    def at2wt(self,atIndex=None,Atomic=None):
        if atIndex == None or Atomic == None:
            atIndex = self.input_DF.index
            Atomic = self.input_DF['At. %']
        weight = []
        for i,j in enumerate(atIndex):
            weight.append(Atomic[i]*Formula(j).mass)
        weight = np.array(weight)  
        norm_weight = CC.normalise(weight)
        wt_DF = pd.DataFrame(data = {'Wt. %' : norm_weight}, index = atIndex)
        return wt_DF
            
        
    def at2ox_wt(self,atIndex=None,Atomic=None):
        if any(atIndex) == None or any(Atomic) == None:
            atIndex = self.input_DF.index
            Atomic = self.input_DF['At. %']
        atIndex = atIndex.drop('O')
        
        charges = np.array([CC.mid(m.element(i).oxistates) for i in atIndex])
        Charges = pd.DataFrame(data={'Oxidation State':charges},index=atIndex)
        Charges.loc['Fe'] = self.FeCharge
        if self.imposeOxidation == True:
            for i in self.OxStates.index:
                try: Charges.loc[i] = self.OxStates.loc[i][0]
                except: next
        print(Charges)
        Charges['O ions'] = Charges['Oxidation State']/2
        
        ox_wtIndex = [str(i+'O'+str(j)) for i,j in zip(atIndex,Charges['O ions'])] 
        
        ox_wt = []
        
        for i,j in zip(atIndex,ox_wtIndex):
            ox_wt.append(float(Atomic[i])*(Formula(j.split('O')[0]).mass + 
                                           (m.O.mass*float(j.split('O')[1]) )))
            
        normOx_Wt = CC.normalise(ox_wt)
        ox_wt_DF = pd.DataFrame(data = {'Ox. Wt. %' : normOx_Wt}, index = ox_wtIndex)
        
        return ox_wt_DF 
   
    def get_ox_wt(self):
        if self.input_DF.columns[0] == 'Wt. %':
            self.ox_wt_DF = CC.at2ox_wt(self,atIndex=self.at_DF.index, 
                                        Atomic=self.at_DF['At. %'])
        return self.ox_wt_DF
    
    def get_wt(self):
        if self.input_DF.columns[0] == 'Ox. Wt. %':
            self.wt_DF = CC.at2wt(atIndex=self.at_DF.index, 
                                        Atomic=self.at_DF['At. %'])
        return self.wt_DF
    
    def get_at(self):
        return self.at_DF


# %%

"""

EC_wt = pd.DataFrame(index = ['O','Si','Mg','Al','Ca','Fe','Ti','Mn','Cr','Na','K'],
                             data = [48.0,29.7,19.5,1.36,0.17,0.33,0.01,0.03,0.13,1.04,0.09],
                             columns = ['Wt. %'])

EC = CC(EC_wt)
EC.get_ox_wt()
EC.get_at()

"""