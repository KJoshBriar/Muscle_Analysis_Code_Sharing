import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import openpyxl
import ISO
import RFE
import combine
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


    #-----Grab info from filename-----+
    print(Filename)
    Subject: str = Filename[0]
    Matching: str = Filename[2].split('.')[0]

    Data = GetTextData(skiprows = 0, columns = {0: 'Time', 1: 'Length', 2: 'Stim', 3: 'Force'})
  
    try:
        FibreLength = float(Data['Stim'][10])
        SarcomereLength = float(Data['Force'][11])
        Diameter = float(Data['Length'][12])
        CSA = np.pi*(Diameter/2)**2
    except:
        print(Error(f"Physical characteristics for {File} not calculated correctly. Check text file for weird lines of data at beginning of file."))

    if test == 'ISO':
        Fibre: str = Filename[1]
        data_start = Data.index[Data['Time'] == '0.00'][0]
        Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
        Baseline_force = Data['Force'][200000:205000].mean()
        Data['Force'] = Data['Force'] - Baseline_force
        Data['Time'] = Data['Time'].div(1000)
        results = {
                'data': Data,
                'Subject': Subject,
                'Matching': Matching,
                'Fibre Length': FibreLength,
                'Fibre': Fibre,
                'Sarcomere Length': SarcomereLength,
                'CSA': CSA}

    elif test == 'RFe':
        Fibre: str = Filename[1]
        data_start = Data.index[Data['Time'] == '0.00'][0]
        Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
        Baseline_force = Data['Force'][90000:95000].mean()
        Data['Force'] = Data['Force'] - Baseline_force
        Data['Time'] = Data['Time'].div(1000)
        results = {
                'data': Data,
                'Subject': Subject,
                'Matching': Matching,
                'Fibre Length': FibreLength,
                'Fibre': Fibre,
                'Sarcomere Length': SarcomereLength,
                'CSA': CSA}
        
    elif test == 'combined':
        Fibre: str = Filename[1]
        data_start = Data.index[Data['Time'] == '0.00'][0]
        Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
        Baseline_force = Data['Force'][200000:205000].mean()
        Data['Force'] = Data['Force'] - Baseline_force
        Data['Time'] = Data['Time'].div(1000)
        results = {
                'data': Data,
                'Subject': Subject,
                'Matching': Matching,
                'Fibre Length': FibreLength,
                'Fibre': Fibre,
                'Sarcomere Length': SarcomereLength,
                'CSA': CSA}              

    return results

def main():
    Directory = filedialog.askopenfilenames()

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

        if 'ISO' in Filename:
            test = 'ISO'
        elif 'rFE' in Filename:
            test = 'RFe'
        else: 
            test = 'combined'
    
        
        results = ReadFile(File = file, Filename = Filename, test=test)
        if test == 'RFe':
            results['Stiffness (pCa 4.5)'], results['Modulus (pCa 4.5)'], results ['Active Specific Force']=RFE.ktrAnalysis(Data = results['data'], Filename = os.path.basename(file), CSA=results['CSA'], Graph = True)
        if test == 'ISO':
            results['Stiffness (pCa 4.5)'], results['Modulus (pCa 4.5)'], results ['Active Specific Force']=ISO.ktrAnalysis(Data = results['data'], Filename = os.path.basename(file), CSA=results['CSA'], Graph = True)
        if test == 'combined':
            results['Stiffness (pCa 4.5)'], results['Modulus (pCa 4.5)'], results ['Active Specific Force'], results['RFEStiffness (pCa 4.5)'], results['RFEModulus (pCa 4.5)'], results ['RFEActive Specific Force']=combine.ktrAnalysis(Data = results['data'], Filename = os.path.basename(file), CSA=results['CSA'], Graph = True)
        individual_test_info[file] = results.copy()
        del individual_test_info[file]['data']

        excel_results.loc[file, individual_test_info[file].keys()] = individual_test_info[file].values()

        # return Final_results
    return excel_results 

#call the main () function to get final_results
excel_results = main()

excel_file_path = '/Users/keaton.briar/Desktop/Outputd4.xlsx'

excel_results.to_excel(excel_file_path, index = False)

os.system(f'open "{excel_file_path}"')

if __name__ == "__main__":
    main()
