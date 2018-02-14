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
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


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
    
    print('Work in progress. Please wait...')
    if all_or_any=='all':
        #create an empty pandas DataFrame 

        #1. Precise search: match all elements in the list
        for item in df_el.Elements.iteritems():
            #print(item[1])
            if set(item[1])==set(elements_list):
                #Add to dataframe
                df_el_filt=pd.concat([df_el_filt,df_el[df_el.index==item[0]]]) #when appending to a DataFrame, remember to reassign!
                #print(df_el[df_el.index==item[0]])
    else:
        #2. Search for any elements in the list 
        ##df_el_filt=filter_by_elements(df_el,elements_list)
        
        for item in df_el.Elements.iteritems(): #if a subset of elements_list, append to the new DataFrame
            if set(item[1]).issubset(elements_list):
                df_el_filt=pd.concat([df_el_filt,df_el[df_el.index==item[0]]]) #when appending to a DataFrame, remember to reassign!
                #print(df_el[df_el.index==item[0]])
    
    print('Work is completed. ')    
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
    
import requests
import urllib

def download_jcamp_from_nist(df, key, minimum_file_size=1000, no_of_files=20):
    """ Download jcamp files from NIST website by calling out CAS values (i.e., the df index) from chemicals DataFrame, df 

    Args:
    df=a pandas DataFrame consisting a CAS no, column as the index.
    size=the minimum file size at which the file will be downloaded.
    no_of_files=the number of files to download.

    Returns:
    automatically saves files into the same directory as this python program.
    print the name of the downloaded files.

    """
    print('Please wait..')
    # NOTE THAT SOME OF THESE FILES TURN out to be empty.
    # Need to check file size before proceeding. Don't download if empty.
    # Files are saved in the same folder as this python script.

    # Download data based on CAS number partially
    
    ##count=0 #This variable keeps track of number of files that have been successfully downloaded
    
    count=0 #initialize counts
    for cas_no in df.index: #the index column of the DataFrame is the CAS number
        # specify CAS no. separately
        url = "http://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C%s&Index=0&Type=IR" % (
            cas_no)

        try:

            # Download only if size is significant. 1000 seems optimum.
            if len(requests.get(url).content) >= minimum_file_size:
                urllib.request.urlretrieve(url, "reference/%s" %(cas_no+'_'+key+'.jcamp'))  # save file according to its cas_no
                                
                count+=1 #increment count if download is successful
                
                #print(count)
                if count==no_of_files:
                    break #exit the for loop
                else: 
                    continue
        except:
            #print('Download unsuccessful. Please check.')
            continue
    print(str(count) + ' ' + key + ' files downloaded' )

"""
#spectra treatent functions
"""

#1. standardize units and normalize 

def normalize_peaks(jcamp_dict):
    # uniformize data, #2:
    # normalize absorbance peaks (y-values) to between 0 and 1.
    jcamp_dict['y'] = (jcamp_dict['y']-min(jcamp_dict['y'])) / max(jcamp_dict['y'])

    return jcamp_dict
    
def standardize_units(jcamp_dict):
        """Given a jcamp dictionary, standardize x and y units to 1/cm and absorbance
           , respectively and return a modified dictionary. 
           Args:
           jcamp_dict=a jcamp dictionary    
           Returns:
           jcamp_dict=a modified jcamp dictionary   
        """

        if jcamp_dict['yunits'] == "ABSORBANCE":
            pass 

        elif jcamp_dict['yunits'] == "TRANSMISSION" or jcamp_dict['yunits'] == "TRANSMITTANCE":
            jcamp_dict['y'] = 2 - np.log10(100*(jcamp_dict['y']+0.001)) #avoid division by zero
            #normalize
            jcamp_dict['y'] = (jcamp_dict['y']-min(jcamp_dict['y'])) / max(jcamp_dict['y'])


        # check xunits: if in microns, change to 1/cm
        if jcamp_dict['xunits'] == "MICROMETERS":
            jcamp_dict['x'] = 10000 / jcamp_dict['x']
            #also change the first and lastx
            jcamp_dict['firstx']=10000/jcamp_dict['firstx']
            jcamp_dict['lastx']=10000/jcamp_dict['lastx']

        # uniformize data, #2:
        # normalize absorbance peaks (y-values) to between 0 and 1.
        #jcamp_dict['y'] = (jcamp_dict['y']-min(jcamp_dict['y'])) / max(jcamp_dict['y'])

        return jcamp_dict

#2. Interpolate spectrum
def interpolate_spectrum(jcamp_dict,xmin=800,xmax=3000,res=0.5):
    """
    Given a jcamp dictionary, interpolate the data points to the specified 
       xmin, xmax, and resolution.
       Args:
       jcamp_dict=a jcamp dictionary    
       xmin=minimum x-value
       xmax=maximum x-value
       res=x-inteval
       Returns:
       returns back a dictionary    

    """    
    #First, check the unit for x-axis before applying interpolation.
    #if it is in micron, change it to wavenumber.
    # check xunits: if in microns, change the limits and resolution to 1/cm
    #if jcamp_dict['xunits'] == "MICROMETERS":
    #    res = 10000 / float(res)
    #    xmin= 10000 / float(xmin)
    #    xmax = 10000 / float(xmax)
        #print('microns')
    
    #dx=(jcamp_dict['lastx']-jcamp_dict['firstx'])/jcamp_dict['npoints']
    dx=(jcamp_dict['x'][-1]-jcamp_dict['x'][0])/jcamp_dict['npoints']
    
    my_list_x=[] #generate x values
    for i in range(jcamp_dict['npoints']):
        #my_list_x.append(jcamp_dict['firstx']+(dx*i))
        my_list_x.append(jcamp_dict['x'][0]+(dx*i))
        
    tck = interpolate.splrep(my_list_x, jcamp_dict['y'], s=0) # s is smoothing
    

    jcamp_dict['x'] = np.arange(xmin,xmax,res) #the second limit is excluded.
    jcamp_dict['y'] = interpolate.splev(jcamp_dict['x'], tck, der=0) #der refers to the derivative
    
    return jcamp_dict



"""
#TEST    
#read jcamp of two spectra data sets
spectrum1=jcamp.JCAMP_reader("1-Ethoxy-4-nitrobenzene_100-29-8") #1st data set
spectrum2=jcamp.JCAMP_reader("Acetanilide 2-chloro-4-tert-butyl-_100141-30-8") #2nd data set

calculate_euclidean_dist(spectrum1,spectrum2)

"""
#baseline subtraction function
#ref.: https://raw.githubusercontent.com/zmzhang/airPLS/master/airPLS.py

#!/usr/bin/python
'''
airPLS.py Copyright 2014 Renato Lombardo - renato.lombardo@unipa.it
Baseline correction using adaptive iteratively reweighted penalized least squares

This program is a translation in python of the R source code of airPLS version 2.0
by Yizeng Liang and Zhang Zhimin - https://code.google.com/p/airpls
Reference:
Z.-M. Zhang, S. Chen, and Y.-Z. Liang, Baseline correction using adaptive iteratively reweighted penalized least squares. Analyst 135 (5), 1138-1146 (2010).

Description from the original documentation:

Baseline drift always blurs or even swamps signals and deteriorates analytical results, particularly in multivariate analysis.  It is necessary to correct baseline drift to perform further data analysis. Simple or modified polynomial fitting has been found to be effective in some extent. However, this method requires user intervention and prone to variability especially in low signal-to-noise ratio environments. The proposed adaptive iteratively reweighted Penalized Least Squares (airPLS) algorithm doesn't require any user intervention and prior information, such as detected peaks. It iteratively changes weights of sum squares errors (SSE) between the fitted baseline and original signals, and the weights of SSE are obtained adaptively using between previously fitted baseline and original signals. This baseline estimator is general, fast and flexible in fitting baseline.


LICENCE
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
'''

from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve

def WhittakerSmooth(x,w,lambda_,differences=1):
    '''
    Penalized least squares algorithm for background fitting
    
    input
        x: input data (i.e. chromatogram of spectrum)
        w: binary masks (value of the mask is zero if a point belongs to peaks and one otherwise)
        lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background
        differences: integer indicating the order of the difference of penalties
    
    output
        the fitted background vector
    '''
    X=np.matrix(x)
    m=X.size
    i=np.arange(0,m)
    E=eye(m,format='csc')
    D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
    W=diags(w,0,shape=(m,m))
    A=csc_matrix(W+(lambda_*D.T*D))
    B=csc_matrix(W*X.T)
    background=spsolve(A,B)
    return np.array(background)

def airPLS(x, lambda_=100, porder=1, itermax=15):
    '''
    Adaptive iteratively reweighted penalized least squares for baseline fitting
    
    input
        x: input data (i.e. chromatogram of spectrum)
        lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background, z
        porder: adaptive iteratively reweighted penalized least squares for baseline fitting
    
    output
        the fitted background vector
    '''
    m=x.shape[0]
    w=np.ones(m)
    for i in range(1,itermax+1):
        z=WhittakerSmooth(x,w,lambda_, porder)
        d=x-z
        dssn=np.abs(d[d<0].sum())
        if(dssn<0.001*(abs(x)).sum() or i==itermax):
            if(i==itermax): print('WARING max iteration reached!') 
            break
        w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
        w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
        w[0]=np.exp(i*(d[d<0]).max()/dssn) 
        w[-1]=w[0]
    return z

def baseline_subtract(data_dict):
    """
    Given a jcamp dictionary, this function performs a baseline subtraction
    Args:
    data_dict=a dictionary that is generated using my_jcamp.py from a jcamp file.
    Return:
    data_dict=a dictionary of treated spectrum.    
    """
    
    data_dict['y'] = data_dict['y'] - airPLS(data_dict['y'])
    
    
    return data_dict #assign new values to the dictionary



#create a function to treat spectra
def treat_spectra(data_dict):
    """Uniformize spectra units (x- and y-axes) to make them all comparable with one another.
    
        Args:
        data_dict=a dictionary that is generated using my_jcamp.py from a jcamp file.
        Return:
        data_dict=a dictionary of treated spectrum.
    """
    
    if data_dict['yunits'] == "ABSORBANCE":
            pass 

    elif data_dict['yunits'] == "TRANSMISSION" or "TRANSMITTANCE":
        data_dict['y'] = 2 - np.log10(data_dict['y'])
        #normalize
        data_dict['y'] = (data_dict['y']-min(data_dict['y'])) / max(data_dict['y'])

    # check xunits: if in microns, change to 1/cm
    if data_dict['xunits'] == "MICROMETERS":
        data_dict['x'] = 10000 / data_dict['x'] #ref_spectra[i]['x']

    # uniformize data, #2:
    # normalize absorbance peaks (y-values) to between 0 and 1.
    data_dict['y'] = (data_dict['y']-min(data_dict['y'])) / max(data_dict['y'])

    # !Still need to uniformize x-axis range and uniformize x-axis intervals.

    return data_dict

    
def pick_peaks(compound_name, x, y):
    from peakutils.peak import indexes as index_utils

    """ Identify peaks maxima. This algorithm is a 1-D search. It doesn't account for the x-values.
        ref.: #https://stackoverflow.com/questions/31016267/peak-detection-in-python-how-does-the-scipy-signal-find-peaks-cwt-function-work
        
        Args:
        compound_name=to be used asa title in the generated plot (string).
        x=x-values (list)
        y=y-values (list)
        Returns:
        uses peakutils package to pick peaks in 1-D and plot out the x,y values, along with the identified peak maxima.
    """

    index_utils = index_utils(np.array(y), thres=0.025 / max(y), min_dist=5)
    # print('Peaks are: %s' % (indexes))

    x_max = []
    y_max = []

    for each in index_utils:
        x_max.append(x[each])
        y_max.append(y[each])

    # plot the x,y pair of the identified maxima data points, i.e., indexes,
    # max_peaks:
    plt.plot(x, y, 'b-')
    plt.plot(x_max, y_max, 'rD', alpha=0.3)

    plt.title(compound_name)
    plt.ylabel("Absorbance")  # All units are transformed to absorbance
    plt.xlabel("Wavenumber, cm-1")  # all units are transformed into 1/cm
    plt.show()
    
# End of functions



    


