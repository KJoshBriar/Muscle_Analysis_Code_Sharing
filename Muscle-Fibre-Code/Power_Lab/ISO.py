
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import openpyxl
from tkinter import filedialog

plt.rcParams.update({
    'font.family': 'times new roman',
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'figure.dpi': 200,
    'legend.edgecolor': 'white',
    'figure.figsize': [6, 4],
    'figure.autolayout': True})

#edit this for how long each thing occurs within the trial
sample_rate = 10000
stiffness_time:int = 65 * sample_rate


class Error(Exception):
    pass

def ReadFile(File: str = None, model: str = None, test: str = None):
    Subject: str = None
    Fibre: str = None
    Muscle: str = None
    def GetTextData(skiprows: int, columns: dict) -> pd.DataFrame():
        return pd.read_table(
                    File,
                    encoding = 'latin1',
                    header = None,
                    delim_whitespace = True,
                    low_memory = False,
                    on_bad_lines='skip',
                    memory_map = True,
                    skiprows = skiprows,
                    usecols = columns.keys(),
                    names = columns.values())
        

    Filename: str = os.path.basename(File)
    #-----Get characteristic info from filename-----+

    if Filename.__contains__('.dat'):
        Filename: list = Filename.split("_")
        Filename = Filename[0:4]
    else:
        Filename = Filename.split("_")

    #-----Grab info from filename-----+
    # Ex: if filename = vcd01_1_SOL
    ## Subject = vcd01
    ## Fibre = 1
    ## Muscle = SOL
    Subject: str = Filename[0]
    Muscle: str = Filename[1]
    Fibre: str = Filename[2].split('.')[0]

    Data = GetTextData(skiprows = 0, columns = {0: 'Time', 1: 'Length', 2: 'Stim', 3: 'Force'})
    # Grab desired fibre characteristic info
    # If any of the following fail, an error message will be printed to screen to warn user
    # Example of failure:
    ## if df['Stim'][10] doesn't equal a number (i.e., 1.50, 0.90, etc) 
    ## then Fibre_Length will fail and a warning message will appear
    #Might be the rows below this that are drawing the error if it isnt pulling the correct row
    try:
        FibreLength = float(Data['Stim'][10])
        SarcomereLength = float(Data['Force'][11])
        Diameter = float(Data['Length'][12])
        CSA = np.pi*(Diameter/2)**2
    except:
        print(Error(f"Physical characteristics for {File} not calculated correctly. Check text file for weird lines of data at beginning of file."))

#change this to the appropriate start time
    data_start = Data.index[Data['Time'] == '0.00'][0]
    Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
    Baseline_force = Data['Force'][200000:205000].mean()
    Data['Force'] = Data['Force'] - Baseline_force
    Data['Time'] = Data['Time'].div(1000)
    #Data['Normalized Length'] = Data['Stim'].div(FibreLength)
            
    return Data, Subject, Muscle, Fibre, FibreLength, SarcomereLength, CSA

def ktrAnalysis(Data: pd.DataFrame = None, Filename: str = None, CSA: float = None, Graph: bool = False) -> tuple[float, float, float, pd.Series, pd.Series, float]:

    def StiffnessAnalysis(Data: pd.DataFrame, stiffness_time: int = stiffness_time):
        StiffnessWindow = range(int(stiffness_time) - 100, int(stiffness_time) + 200)
        ForceWindow = range(int(stiffness_time) - 5001, int(stiffness_time) - 1)
        dF = (Data['Force'][StiffnessWindow]).max() - (Data['Force'][ForceWindow]).mean()
        #dLo = (Data['Normalized Length'][StiffnessWindow]).max() - (Data['Normalized Length'][ForceWindow]).mean()
        dLo = (Data['Stim'][StiffnessWindow]).max() - (Data['Stim'][ForceWindow]).mean()
        Strain = ((Data['Stim'][StiffnessWindow]).max() - (Data['Stim'][ForceWindow]).mean())/(Data['Stim'][ForceWindow].mean())
        Stiffness = dF/dLo
        #could make strain here a set value, regardless it should pull whatever is in the length column
        Modulus = (dF/CSA)/Strain


        return Stiffness, Strain, Modulus

    Stiffness, Strain, Modulus = StiffnessAnalysis(Data = Data)
    Peak_force = Data['Force'][600000:610000].mean()

    if Graph == True:
        fig = plt.figure()
        
        plt.plot(Data['Time'], Data['Force'], color = 'black', label = 'Raw')
        plt.plot(Data['Time'][600000:610000], Data['Force'][600000:610000], color = 'blue', label = 'Peak')
        plt.plot(Data['Time'][200000:205000], Data['Force'][200000:205000], color = 'Yellow', label = 'Baseline')
        plt.plot(Data['Time'][650000:655000], Data['Force'][650000:655000], color = 'Orange', label = 'Stiff')

        plt.ylabel('Force (mN)')
        plt.xlabel('Time (s)')
        plt.text(
            x = 0.5, y = 0.1,
            s = f'Active Specific Force = {(Peak_force/CSA):.2f}',
            transform = plt.gca().transAxes,
            horizontalalignment = 'center',
            verticalalignment = 'center')
        plt.title(Filename)
        plt.legend()
        plt.show()

    return Stiffness, Modulus, Peak_force/CSA
