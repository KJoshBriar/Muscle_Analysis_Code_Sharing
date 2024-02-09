import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import openpyxl
import KTR_STIFF
import Unloaded_Short
from scipy.optimize import differential_evolution, curve_fit
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

class Error(Exception):
    pass

def ReadFile(File: str = None, model: str = None, Filename = None, test: str = None):
    
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
        

    # Filename: str = os.path.basename(File)
    # #-----Get characteristic info from filename-----+

    # if Filename.__contains__('.dat'):
    #     Filename: list = Filename.split("_")
    #     Filename = Filename[0:4]
    # else:
    #     Filename = Filename.split("_")

    #-----Grab info from filename-----+
    print(Filename)
    Subject: str = Filename[0]
    Fibre: int = Filename[1]
    Muscle: str = Filename[2].split('.')[0]

    Data = GetTextData(skiprows = 0, columns = {0: 'Time', 1: 'Length', 2: 'Stim', 3: 'Force'})
  
    try:
        FibreLength = float(Data['Stim'][10])
        SarcomereLength = float(Data['Force'][11])
        Diameter = float(Data['Length'][12])
        CSA = np.pi*(Diameter/2)**2
    except:
        print(Error(f"Physical characteristics for {File} not calculated correctly. Check text file for weird lines of data at beginning of file."))
    if test == 'ktr':
        data_start = Data.index[Data['Time'] == '25000.00'][0]
        Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
        Baseline_force = Data['Force'][1:500].mean()
        Data['Force'] = Data['Force'] - Baseline_force
        Data['Time'] = Data['Time'].div(1000)
        results = {
                'data': Data,
                'Subject': Subject,
                'Muscle': Muscle,
                'Fibre Length': FibreLength,
                'Fibre': Fibre,
                'Sarcomere Length': SarcomereLength,
                'CSA': CSA}

    elif test == 'unloaded':
        Velocity: str = Filename[3]
        data_start = Data.index[Data['Time'] == '0.00'][0]
        Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
        Baseline_force = Data['Force'][1:500].mean()
        Data['Force'] = Data['Force'] - Baseline_force
        Data['Time'] = Data['Time'].div(1000)
        results = {
                'data': Data,
                'Subject': Subject,
                'Muscle': Muscle,
                'Velocity': Velocity,
                'Fibre Length': FibreLength,
                'Fibre': Fibre,
                'Sarcomere Length': SarcomereLength,
                'CSA': CSA}
    #Data['Normalized Length'] = Data['Stim'].div(FibreLength)
            
    return results

def main():
    Directory = filedialog.askopenfilenames()
    #Directory = 'FILE DIRECTORY HERE'

    # AllFiles = []
    # for root, subdirs, files in os.walk(Directory):
    #    for file in files:
    #        if 'Store' in file:
    #            continue
    #        AllFiles.append(os.path.join(root, file))
    # AllFiles = sorted(AllFiles)


    final_results = []
    excel_results=pd.DataFrame()
    for file in Directory:
        individual_test_info = {}
        Filename = os.path.basename(file)
        if Filename.__contains__('.dat'):
            Filename: list = Filename.split("_")
            Filename = Filename[0:4]
        else:
            Filename = Filename.split("_")

        if len(Filename) == 4:
            test='unloaded'
        else:
            test='ktr'
        
        results = ReadFile(File = file, Filename = Filename, test=test)
        if test == 'unloaded':
            results['Time of force redevelopment'], results['Isometric Force']=Unloaded_Short.VelocityAnalysis(Data = results['data'], Filename = os.path.basename(file), CSA=results['CSA'], Graph = False)
        if test == 'ktr':
            results['Stiffness (pCa 4.5)'], results['Relaxing Stiffness'], results['Modulus (4.5)'], results['Modulus (Relaxing)'], results['Normalized Modulus'], results['ktr (pCa 4.5)'], results ['Goodness of fit'], results ['Active Specific Force'], results['KTR Delay'] =KTR_STIFF.ktrAnalysis(Data = results['data'], Filename = os.path.basename(file), CSA=results['CSA'], Graph = False)
        individual_test_info[file] = results.copy()
        del individual_test_info[file]['data']

        excel_results.loc[file, individual_test_info[file].keys()] = individual_test_info[file].values()

        # return Final_results
    return excel_results 

#call the main () function to get final_results
excel_results = main()

# excel_results.to_excel('/Users/keaton.briar/Desktop/resultsbalh.xlsx')
#excel_file_path = filedialog.asksaveasfilename(defaultextension= '.xlsx')
excel_file_path = '/Volumes/KJBRIAR/AR02/Output.xlsx'


excel_results.to_excel(excel_file_path, index = False)

os.system(f'open "{excel_file_path}"')

if __name__ == "__main__":
    main()
