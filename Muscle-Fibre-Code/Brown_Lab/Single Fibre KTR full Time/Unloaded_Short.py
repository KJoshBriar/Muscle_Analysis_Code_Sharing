
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings
import openpyxl
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

#edit this for how long each thing occurs within the trial
sample_rate = 10000
UnloadedShortening_Time:int = 8 * sample_rate

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
    print(Filename)
    Subject: str = Filename[0]
    Fibre: int = Filename[1]
    Muscle: str = Filename[2]
    Velocity: str = Filename[3].split('.')[0]

    Data = GetTextData(skiprows = 0, columns = {0: 'Time', 1: 'Length', 2: 'Stim', 3: 'Force'})
  
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
    Baseline_force = Data['Force'][1:500].mean()
    Data['Force'] = Data['Force'] - Baseline_force
    Data['Time'] = Data['Time'].div(1000)
    #Data['Normalized Length'] = Data['Stim'].div(FibreLength)
            
    return Data, Subject, Muscle, Fibre, Velocity, FibreLength, SarcomereLength, CSA

def VelocityAnalysis(Data: pd.DataFrame = None, Filename: str = None, CSA: float = None, Graph: bool = False) -> tuple[float, float, float, pd.Series, pd.Series, float]:
    def Force_Positive(Data: pd.DataFrame, UnloadedShortening_Time: int = UnloadedShortening_Time):
        positive_window = range(UnloadedShortening_Time - 1, UnloadedShortening_Time + 10000)
        def handle_dem_mouse_clickies(): #Code to define interactable graph
            click_x=round(coords[0][0]*10000)
         #print (click_x-250000)
            click_y=round(coords[0][1])

            ax.plot(Data['Time'][click_x:UnloadedShortening_Time+10000], Data['Force'][click_x:UnloadedShortening_Time+10000], color='red') #plot graph data from the pre-selected model data of the KTR time
            ax.axvline(x=Data['Time'][click_x], ymin = 0, ymax = 1, color = 'Green', linestyle = '-', linewidth = 2) #show cursor line
            fig.canvas.draw()
            plt.show()

            return click_x
        for index in positive_window:
            fig=plt.figure()
            ax=fig.add_subplot()

            ax.plot(Data['Time'][positive_window], Data['Force'][positive_window], color='black')
            coords=plt.ginput(n=1, show_clicks=True, timeout=9999) #plots the data and then allows you to click a point
            if coords:
                Positive_Force=handle_dem_mouse_clickies() #given a mouse click, we pull the function from above for the interactable graph to get a defined start time to the KTR data
                #Pos_Time=VelocityAnalysis(Positive_Force-80000) #uses the manually selected ktr time to run the KTR fitting
            #add a print
                return Positive_Force  # Return the index when force becomes positive
        return -1  # Return -1 if no positive force is found

    Positive_force = Force_Positive(Data = Data)
    Peak_force = Data['Force'][70000:75000].mean()
   
    if Graph == True:
        fig = plt.figure()
        
        plt.plot(Data['Time'], Data['Force'], color = 'black', label = 'Raw')
        plt.plot(Data['Time'][70000:75000], Data['Force'][70000:75000], color = 'blue', label = 'Peak')
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


    return (Positive_force-80000)/10, Peak_force/CSA


# def main():
#     Directory = filedialog.askopenfilenames()
#     #Directory = 'FILE DIRECTORY HERE'

#     # AllFiles = []
#     # for root, subdirs, files in os.walk(Directory):
#     #    for file in files:
#     #        if 'Store' in file:
#     #            continue
#     #        AllFiles.append(os.path.join(root, file))
#     # AllFiles = sorted(AllFiles)


#     final_results = []
#     for file in Directory:

#         Data, Subject, Muscle, Fibre, Velocity, FibreLength, SarcomereLength, CSA = ReadFile(File = file)
#         Positive_Force, Peak_force = VelocityAnalysis(Data = Data, Filename = os.path.basename(file), CSA=CSA, Graph = True)

#         OrganizedData = {
#             'Subject': Subject,
#             'Muscle': Muscle,
#             'Velocity': Velocity,
#             'Fibre Length': FibreLength,
#             'Fibre': Fibre,
#             'Sarcomere Length': SarcomereLength,
#             'CSA': CSA,
#             'Time of force redevelopment': Positive_Force-80000,
#             'Isometric Force': Peak_force / CSA}
        
#         final_results.append(OrganizedData)
    
#     for test in final_results:
#         print(test)

#     # return Final_results
#     return final_results 

# #call the main () function to get final_results
# final_results = main()

# # your code for generating the final results
# for File in final_results:
#      excel_results = pd.DataFrame(columns = [key for key in File.keys()])

# for idx, results in enumerate(final_results):
#     for column, data in results.items():
#         excel_results.loc[idx, column] = data

# excel_file_path = filedialog.asksaveasfilename(defaultextension= '.xlsx')

# excel_results.to_excel(excel_file_path, index = False)

# os.system(f'open "{excel_file_path}"')

# if __name__ == "__main__":
#     main()