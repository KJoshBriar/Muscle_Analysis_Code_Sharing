
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import openpyxl
from tkinter import filedialog
import tkinter as tk
import tkinter.simpledialog as simpledialog

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

class Error(Exception):
    pass

def ReadFile(File: str = None, model: str = None, test: str = None):
    Subject: str = None
    Fibre: str = None
    Muscle: str = None
    def GetTextData(skiprows: int, columns: dict) -> pd.DataFrame:
        return pd.read_table(
                    File,
                    encoding = 'latin1',
                    header = None,
                    sep=r'\s+',
                    engine='python',
                    on_bad_lines='skip',
                    memory_map = True,
                    skiprows = skiprows,
                    usecols = columns.keys(),
                    names = columns.values())
        

    Filename: str = os.path.basename(File)
    #-----Get characteristic info from filename-----+
    Fibre = simpledialog.askstring(title="Fibre Input", prompt=f"Enter fibre number for: {Filename}")

    Data = GetTextData(skiprows = 0, columns = {0: 'Time', 1: 'Length', 2: 'Stim', 3: 'Force'})

    try:
        FibreLength = float(Data['Stim'][10])
        SarcomereLength = float(Data['Force'][11])
        Diameter = float(Data['Length'][12])
        CSA = np.pi*(Diameter/2)**2
    except:
        print(Error(f"Physical characteristics for {File} not calculated correctly. Check text file for weird lines of data at beginning of file."))

#change this to the appropriate start time
    data_start = Data.index[Data['Time'] == '5000.00'][0]
    Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
    Baseline_force = Data['Force'][1:500].mean()
    Data['Force'] = Data['Force'] - Baseline_force
    Data['Time'] = Data['Time']/(1000)
    #Data['Normalized Length'] = Data['Stim'].div(FibreLength)
            
    return Data, Subject, Muscle, Fibre, FibreLength, SarcomereLength, CSA

def SFAnalysis(Data: pd.DataFrame = None, Filename: str = None, CSA: float = None, Graph: bool = False) -> tuple[float, float, float, pd.Series, pd.Series, float]:

    Peak_force = Data['Force'][600000:610000].mean()

    if Graph == True:
        fig = plt.figure()
        
        plt.plot(Data['Time'], Data['Force'], color = 'black', label = 'Raw')
        plt.plot(Data['Time'][600000:610000], Data['Force'][600000:610000], color = 'Orange', label = 'Peak')
        plt.plot(Data['Time'][1:500], Data['Force'][1:500], color = 'Yellow', label = 'Baseline')

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
        plt.close('all') 

    return Peak_force/CSA
