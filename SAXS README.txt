From the SAXS object documentation:

""" 
SAXS.py is a module to create a SAXS scattering object consisting of a Pandas DataFrame set
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

The github repository also contains "testing SAXS class.py", a module that demonstrates the capabilities of the SAXS object. This module uses the included files 
"80C scattering data.csv" and "RT scattering data.csv", so be sure to download those as well and keep them in the same directory as the "testing SAXS class.py" module. 