# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 12:17:39 2018

@author: rayg

Stores common functions:
    1. Create a periodic table dictionary
    2. standardize units and normalize 
    3. Interpolate spectrum
    
"""
import pandas as pd
import re
#import numpy as np
#import matplotlib.pyplot as plt
#from scipy import interpolate
#import my_jcamp as jcamp

#1. Create a periodic table dictionary
#create an ordered dictionary. In this case, order of elements is important.

import collections
# periodic_table['element_name']=['atomic_number','atomic_mass']
periodic_table = collections.OrderedDict()  # initialize dictionary.

def create_periodic_table(file_name):
    """ Read a text file line-by-line and extract information for: atomic number, atomic symbol, and relative atomic mass 

        The author is aware of at least one python periodic elements package (e.g., https://pypi.python.org/pypi/periodictable),
        but chose to implement a different approach.
    """

    file = open(file_name, 'r')  # open a file to read

    for line in file:  # read file line-by-line, extract values, and store into a dictionary

        if line.startswith('Atomic Number'):
            line = line.replace(" ", "")
            value1 = line.split('=', 1)

        if line.startswith('Atomic Symbol'):
            line = line.replace(" ", "")
            key = line.split('=', 1)

        if line.startswith('Relative Atomic Mass'):
            line = line.replace(" ", "")
            value2 = line.split('=', 1)

            # Store values into dictionary only if the element is unique.
            # Otherwise, the atomic mass will be override by isotope's atomic
            # mass.
            if key[1].strip() not in list(periodic_table.keys()):
                periodic_table[key[1].strip()] = [
                    value1[1].strip(), value2[1].strip()]
            else:
                pass

    return periodic_table

def extract_unique_elements(df):
    """Given a pandas DataFrame containing a column of chemical formula, return the DataFrame 
        with an additional column containing unique elements in the formula called df['Elements']
    
        Args:
        df_Formula=a column of chemical formula in a pandas DataFrame
        Return:
        a list of unique elements
    """

    #1a. check for formulas that contain ')n' and remove it
    ##df[df['Formula'].str.contains('\)n')==True]['Formula'].str.replace('\)n','')   
    df.loc[df['Formula'].str.contains('\)n')==True,'Formula']=df[df['Formula'].str.contains('\)n')==True]['Formula'].str.replace('\)n','') 
    #1b. check for formulas that contain '.xW' and remove it
    df.loc[df['Formula'].str.contains('xW')==True,'Formula']=df[df['Formula'].str.contains('xW')==True]['Formula'].str.replace('x','')
    #1c. check for formulas that contain '.XH2O' and remove it
    df.loc[df['Formula'].str.contains('XH2O')==True,'Formula']=df[df['Formula'].str.contains('XH2O')==True]['Formula'].str.replace('X','')
    #1d. check for formulas that contain 'xC6' and remove it
    df.loc[df['Formula'].str.contains('xC6')==True,'Formula']=df[df['Formula'].str.contains('xC6')==True]['Formula'].str.replace('x','')

    
    #2. apply regex to extract unique elements in chemical formula
    df['Elements'] = df['Formula'].str.replace('[^a-z,A-Z]+', '') #Extract only alphabets and drop numbers
    #3. df['Elements'] = df['Elements'].str.findall('[A-Z][^A-Z]*') 
    df['Elements'] = df['Elements'].str.findall('[A-Z][^A-Z]*') 
    return df

def shorten_df_by_elements_list(df_el, elements_list, all_or_any):
    """
    Reduce number of entries in DataFrame by keeping only those that contain elements indicated on elements_list
    
    Args:
    df_el=the DataFrame containing the distinct elements in each entry
    all_or_any=specify if the output DataFrame should contain all or any of the elements indicated on elements_list
    Returns:
    df_el_filt=pandas DataFrame
    """
    df_el_filt=pd.DataFrame()
    
    if all_or_any=='all':
        #create an empty pandas DataFrame 

        #1. Precise search: match all elements in the list
        for item in df_el.Elements.iteritems():
            #print(item[1])
            if set(item[1])==set(elements_list):
                #Add to dataframe
                df_el_filt=pd.concat([df_el_filt,df_el[df_el.index==item[0]]]) #when appending to a DataFrame, remember to reassign!
                print(df_el[df_el.index==item[0]])
    else:
        #2. Search for any elements in the list 
        ##df_el_filt=filter_by_elements(df_el,elements_list)
        
        for item in df_el.Elements.iteritems(): #if a subset of elements_list, append to the new DataFrame
            if set(item[1]).issubset(elements_list):
                df_el_filt=pd.concat([df_el_filt,df_el[df_el.index==item[0]]]) #when appending to a DataFrame, remember to reassign!
                print(df_el[df_el.index==item[0]])
        
    return df_el_filt

def calc_molec_weight(nist_chem_list, periodic_table):
    """ Calculate molecular weight, Mw, given a chemical formula, df.Formula 
    #ref.: https://stackoverflow.com/questions/41818916/calculate-molecular-weight-based-on-chemical-formula-using-python

    Args:
    nist_chem_list=a pandas DataFrame from which Mw will be calculated.
    periodic_table=a periodic table dictionary from which the atomic weight of each element can be retrieved.

    Returns:
    the same dataframe, but with an extra column added. This extra column is labelled as 'Mw'.

    """

    #1a. check for formulas that contain '.xW' and remove it
    nist_chem_list.loc[nist_chem_list['Formula'].str.contains('xW'),'Formula']=nist_chem_list[nist_chem_list['Formula'].str.contains('xW')==True]['Formula'].str.replace('x','')
    #1b. check for formulas that contain '.XH2O' and remove it
    nist_chem_list.loc[nist_chem_list['Formula'].str.contains('XH2O')==True,'Formula']=nist_chem_list[nist_chem_list['Formula'].str.contains('XH2O')==True]['Formula'].str.replace('X','')
    #1c. check for formulas that contain 'xC6' and remove it
    nist_chem_list.loc[nist_chem_list['Formula'].str.contains('xC6')==True,'Formula']=nist_chem_list[nist_chem_list['Formula'].str.contains('xC6')==True]['Formula'].str.replace('x','')

    
    for i, each in enumerate(nist_chem_list.Formula):

        # separate atomic symbols from atomic ratios
        total_weight = 0  # initialize total weight at the beginning of each loop
        s = re.findall('([A-Z][a-z]?)([0-9]*)', each)

        for elem, count in s: #unpack s into elemental symbol and an integer
            if count == '':  # For singular elements (i.e., contains no integers)
                count = 1
            else:
                pass
            try: #calculate molecular weight by multiplying no. of elements to element's atomic mass
                total_weight += int(count) * \
                    float(re.sub('[(#)]', '', periodic_table[elem][1]))
            except: 
                print("error:", nist_chem_list.iloc[i,1])
                continue

        # Append total weight to a new column in pandas DataFrame
        nist_chem_list.loc[nist_chem_list.index[i], 'Mw'] = total_weight
        
    return nist_chem_list
    




    


