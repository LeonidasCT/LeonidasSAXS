""" 
A module to create a SAXS scattering object consisting of a Pandas DataFrame set
of Q and I(Q) data for the full time length of a scattering experiment. The 
DataFrame will contain multiple curves, each consisting of a set of Q and I(Q)
data. The module then applies built-in methods to find the characteristic length,
invariant, intensity, I(Q) maximum, and perform curve fitting to a gaussian 
curve. 

Note: For now, data input must be as a Pandas Data Frame, with Q as the column
labels, aging time as the index labels, and each set of I(Q) data as a row 
corresponding to an aging time index. The format looks as follows:
    for two scattering curves at aging time 8 minutes and 12 minutes, with Q 
    values 0.003 ad 0.005, the DataFrame is:
        
        0.003  0.005
    8   7000   9000
    12  7483   9126

The data frame must be saved as a .csv file using df.to_csv('filename.csv') to 
be exported from Jupyter or another server, but can also be input as a data frame.

It is also recommended that you import the data that you wish to analyze using 
this module as a Pandas DataFrame. This will allows you to examine the data, 
indices, and row/columns lengths simultaneously, which is helpful when defining
limits in some of the limit-bound functions present in this module.

The SAXS module is intialized as follows:

import SAXS

The SAXS object is called as follows:

data = SAXS.ScatteringData(data = variable name or None, csv = False/True, csv_path = 'C:\csv\path_example.csv' if csv=True)

Future functionality will allow for multiple data type inputs with conversion 
to Pandas DF, and eventually a conversion to server-based database storage 
using JSON and Mondo DB.
"""

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt

class ScatteringData: 
    """Defines a class with methods to analyze a database of scattering data."""
    def __init__(self, data = None, csv = False, csv_path = ''):
        """ Creates a new set of data from a Data Frame to use for analysis. 
        If csv is false, the function looks for a python DataFrame. If csv is 
        True, the function takes a file name as an argument and uses that to create
        the data frame for analysis.
        """
        if csv == False:
            self._data = data
        elif csv == True:
            self._data = pd.read_csv(csv_path,header=0,index_col=0,dtype=np.float64)
        else:
            raise ValueError("csv boolean must be set to either True or False (default False).")
    
    def intensity(self):
        """ Computes the intensity curves in a given data frame using the 
        intensity equation, I* = I(Q)*Q**2 and returns them in a data frame.
        """
        intensity = np.empty((self._data.shape[0],self._data.shape[1]))
        q = np.empty(self._data.shape[1])
        for i in range(self._data.shape[1]):
            q[i] = float(self._data.columns[i])
        for i in range(self._data.shape[0]): #loop over rows (tage)
            for j in range(self._data.shape[1]): #loop over columns (Q)
                intensity[i,j] = self._data.iloc[i,j]*q[j]**2 #fills intensity array with intensities computed from original data.
        dfi = pd.DataFrame(intensity, index = self._data.index, columns = q)
        return dfi
    
    def IQmax(self, lowlim = 0, upplim = 1000, peak_positions = False,char_length=False):
        """
        Finds the I(Q) maximum of a set of curves given input data. 
        
        Uses the np.amax method to find the maximum of the intensity curve at a 
        given aging time, then returns a 2-columns numpy array with aging time
        and the I(Q) values from the input data corresponding to the intensity 
        maximum. This function also calculates the peak positions (Q values)
        corresponding to each found I(Q) maximum. This function can also calculate
        the characteristic length (CL) from the equation CL = (2*pi)/Q. 
        
        The order of output for the results is I(Q)max,peak_positions,char_length.
        
        Note : Computationally heavy code. This code currently hard codes
        the aging time limits to include the entire range. Users will need to
        slice or index the resulting numpy arrays. 

        Parameters
        ----------
        lowlim : int, optional
            The lower limit Q index to search for the maximum. The default is 0.
        upplim : int, optional
            The upper limit Q index to search for the maximum. The default is 1000.
        peak_positions : bool, optional
            If True, outputs a 2 x len(aging time) numpy array that contains the 
            Q values corresponding to the calculated I(Q) maximum values. The default is False.
        char_length : bool, optional
            If True, divides 2*pi by the second column (index 1) of the peak 
            positions result array and outputs a 2 x len(aging time) numpy
            array containing the results.
            
        Returns
        -------
        IQmax : numpy array
            A 2 x len(aging time) numpy array that has the aging time in column 0 
            and the corresponding I(Q)  maximum in column 1.
        pp : numpy array, optional
            A 2 x len(aging time) numpy array that has the aging time in column 0
            and the corresponding peak position (Q value) at which there is an 
            I(Q) maximum in column 1.
        cl : numpy array, optional
            A 2 x len(aging time) numpy array that has the aging time in column 0
            and the corresponding characteristic length value at which there is an 
            I(Q) maximum in column 1.
        """
        intensity = self.intensity()
        IQmax = np.empty([self._data.shape[0],2])
        if peak_positions == True:
            pp = np.empty([self._data.shape[0],2]) #initialize peak position array
        if char_length == True:
            cl = np.empty([self._data.shape[0],2]) # initialize characteristic length array
        for i in range(self._data.shape[0]): #loop over rows corresponding to aging times
            IQmax[i,0] = float(self._data.index[i]) #save aging time to first column of results array
            if peak_positions == True:
                pp[i,0] = float(self._data.index[i]) #save aging time to first column of peak position array
            if char_length == True:
                cl[i,0] = float(self._data.index[i]) # save aging time to first column of characteristic length array
            for j in range(self._data.shape[1]): #loop over columns corresponding to Q values
                if intensity.iloc[i,j] == np.amax(intensity.iloc[i,lowlim:upplim]): #if row value is a maximum of the whole row in a given range
                    IQmax[i,1] = self._data.iloc[i,j] #save the maximum I(Q) value corresponding to the maximum intensity value's indices. 
                    if peak_positions == True:
                        pp[i,1] = float(self._data.columns[j]) #save the Q value corresponding to the j (column) index where there exists an I(Q) maximum in the defined limits.
                    if char_length == True:    
                        cl[i,1] = (2*np.pi)/float(self._data.columns[j]) # save 2*pi/Q value corresponding to the j (column) index where there exists an I(Q) maximum in the defined limits.
        if char_length == True and peak_positions == True:
            return IQmax,pp,cl
        elif peak_positions == True:
            return IQmax,pp
        elif char_length == True:
            return IQmax,cl
        else:
            return IQmax
                
    def invariant(self, lowlim = 0, upplim = 1000):
        """
        Uses the intensity to calculate the invariant from the equation 
        invariant = integral(intensity) from zero to infinity. 
        Optionally takes arguments for the lower and upper limits of integration.

        Parameters
        ----------
        lowlim : int, optional
            The Q data index to use as the lower integration limit when calculating the invariant. The default is 0.
        upplim : int, optional
            The Q data index to use as the upper integration limit when calculating the invariant. The default is 1000.
            
        Returns
        -------
        invariant : a 2D numpy array with column 0 as the aging time and column 1 as the corresponding invariant.

        """
        intensity = self.intensity()
        try:
            from scipy import integrate
        except:
            raise ImportError("You need to have the scipy module installed for this code to function properly.")
        invariant = np.empty([self._data.shape[0],2])
        for i in range(self._data.shape[0]):
            invariant[i,0] = self._data.index[i]
            invariant[i,1] = integrate.simps(intensity.iloc[i,lowlim:upplim],intensity.columns[lowlim:upplim]) # calculating the integral using simpson's 1/3 method confined to the limits defined using lowlim, upplim.
        return invariant
    
    def curve_fit(self,lowlim = 0, upplim = 1000, fitfunc = 'Guinier', parameters = False):
        """
        Uses the scipy.optimize curve_fit function to fit the input data to some
        commonly-used models for fitting scattering data to known structures. 

        Parameters
        ----------
        lowlim : int, optional
            The Q data index to use as the lower integration limit when fitting the curve. The default is 0.
        upplim : TYPE, optional
            The Q data index to use as the upper integration limit when fitting the curve. The default is 1000.
        fitfunc : str, optional
            The function to use for fitting. This method also supports the ability to use a custom function. 
            The default is 'Guinier'.
        parameters : bool, optional
            If True, outputs a numpy array that contains the fitting parameters 
            found using scipy.optimize.curve_fit.The default is False.

        Returns
        -------
        fit : numpy array
            A variable size numpy array that contains the fitted data points for plotting
            or analysis.
        params : numpy array, optional
            A variable size numpy array that contains the fitted parameters associated
            with the selected (or input) function.
        """
        #check for scipy install
        try:
            from scipy.optimize import curve_fit
        except:
            raise ImportError("You need to have the scipy module installed for this code to function properly.")
            
        def guinier(q, I0, Rg):
            """
            defining the guinier equation for use in curve fitting code.
            """
            I = I0*np.exp(-(Rg**2*q**2)/3)
            return I
        #data_bound = len(self._data.columns[lowlim:upplim])
        ## Initialize numpy arrays 
        q = np.empty(self._data.shape[1])
        fit = np.empty([self._data.shape[0],self._data.shape[1]])
        for i in range(self._data.shape[1]):
            q[i] = float(self._data.columns[i]) #create q vector (x data for fitting)
        for j in range(self._data.shape[0]):
            if fitfunc == 'Guinier':
                if j == 0:
                    popt,pcov = curve_fit(guinier,q[lowlim:upplim],self._data.iloc[j,lowlim:upplim])
                    if parameters == True:
                        params = np.empty([self._data.shape[0],len(popt)+1]) #define an empty array that has the exact number of rows as fitting params from the function, plus 1 additional column for the aging time 
                        for k,kk in enumerate(popt):
                            params[j,0] = self._data.index[j]
                            params[j,k+1] = kk
                else:
                    popt,pcov = curve_fit(guinier,q[lowlim:upplim],self._data.iloc[j,lowlim:upplim])
                    if parameters == True:
                        for k,kk in enumerate(popt):
                            params[j,0] = self._data.index[j]
                            params[j,k+1] = kk
            for k in range(self._data.shape[1]): #loop over columns
                fit[j,k] = guinier(q[k],*params[j,1:]) #use the calculated parameters to fit the data and save to array
        if parameters == True:
            #return fit,params
            return fit[:,lowlim:upplim],params
        elif parameters != True:
            #return fit
            return fit[lowlim:upplim]
        
    def get_data(self, aging_time=False):
        """
        Simple method to extract the Q values, IQ values, and aging times as
        numpy arrays for analysis. Numpy arrays are extracted in the lowest
        order possible (for example, tage is a 1D numpy vector).

        Parameters
        ----------
        aging_time : bool, optional
        If True, returns the aging time vector.
        
        Returns
        -------
        q : numpy array
            An vector of Q values.
        iq : numpy array
            An array of I(Q) values as a function of aging time.
        tage : numpy array, optional
            A vector of aging time values.
        """
        q = np.empty(self._data.shape[1])
        iq = np.empty([self._data.shape[0],self._data.shape[1]])
        tage = np.empty(self._data.shape[0])
        for j in range(self._data.shape[0]):
            tage[j] = float(self._data.index[j])
            for k in range(self._data.shape[1]):
                iq[j,k] = self._data.iloc[j,k]
        for i in range(self._data.shape[1]):
            q[i] = float(self._data.columns[i])
        if aging_time == True:
            return q,iq,tage
        elif aging_time == False:
            return q,iq