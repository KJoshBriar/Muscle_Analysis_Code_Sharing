
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
stiffness_time:int = 30 * sample_rate
ktr_start_time:int = int(31.016 * sample_rate)
ktr_end_time:int = int(36 * sample_rate)
relax_stiff_time:int = 60 * sample_rate

class Error(Exception):
    pass

def ReadFile(File: str = None, model: str = None, test: str = None):
    """
    Import each text file and grab characteristic information (depending on model and test)\n 
    - Single fibre
        - Fibre length
        - Sarcomere length
        - Fibre diameter (to calculate CSA)
    
    Also obtain information from filename (depending on model and test)\n
    - Subject ID
    - Muscle
    - Fibre #
    
    Returns:
    - Dataframe containing only real data from text files we need for further analyses
    - Characteristic info obtained from text file and filename
    """
    Subject: str = None
    Fibre: str = None
    Muscle: str = None
    def GetTextData(skiprows: int, columns: dict) -> pd.DataFrame():
        return pd.read_table(
                    File,
                    encoding = 'latin1',
                    header = None,
                    sep='\s+',
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
    print(Filename)
    Subject: str = Filename[0]
    Fibre: int = Filename[1]
    Muscle: str = Filename[2]
    Level: str = Filename[3].split('.')[0]

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
    data_start = Data.index[Data['Time'] == '25000.00'][0]
    Data = pd.DataFrame(Data[data_start:], dtype = float).reset_index(drop = True)
    Baseline_force = Data['Force'][1:500].mean()
    Data['Force'] = Data['Force'] - Baseline_force
    Data['Time'] = Data['Time'].div(1000)
    #Data['Normalized Length'] = Data['Stim'].div(FibreLength)
            
    return Data, Subject, Muscle, Fibre, Level, FibreLength, SarcomereLength, CSA

def ktrAnalysis(Data: pd.DataFrame = None, Filename: str = None, CSA: float = None, Graph: bool = False) -> tuple[float, float, float, pd.Series, pd.Series, float]:
    """
    model ktr \n
    Returns:
     - Stiffness (i.e., \u0394 force / \u0394 normalized length)
     - ktr
     - ktr goodness of fit
     - Estimated x and y data from modeling (can be graphed with real data to visually inspect fit)
     - Average force over the final 500 ms of test
    """

    def StiffnessAnalysis(Data: pd.DataFrame, stiffness_time: int = stiffness_time):
        StiffnessWindow = range(int(stiffness_time) - 100, int(stiffness_time) + 200)
        ForceWindow = range(int(stiffness_time) - 5001, int(stiffness_time) - 1)
        dF = (Data['Force'][StiffnessWindow]).max() - (Data['Force'][ForceWindow]).mean()
        #dLo = (Data['Normalized Length'][StiffnessWindow]).max() - (Data['Normalized Length'][ForceWindow]).mean()
        dLo = (Data['Stim'][StiffnessWindow]).max() - (Data['Stim'][ForceWindow]).mean()
        Strain = ((Data['Stim'][StiffnessWindow]).max() - (Data['Stim'][ForceWindow]).mean())/(Data['Stim'][ForceWindow].mean())
        Stiffness = dF/dLo
        Modulus = (dF/CSA)/Strain


        return Stiffness, Strain, Modulus
    
   # def RStiffnessAnalysis(Data: pd.DataFrame, relax_stiff_time: int = relax_stiff_time):
    #    RStiffnessWindow = range(int(relax_stiff_time) - 100, int(relax_stiff_time) + 200)
     #   RForceWindow = range(int(relax_stiff_time) - 2001, int(relax_stiff_time) - 1)
      #  RdF = (Data['Force'][RStiffnessWindow]).max() - (Data['Force'][RForceWindow]).mean()
       # RdLo = (Data['Stim'][RStiffnessWindow]).max() - (Data['Stim'][RForceWindow]).mean()
        #RStrain = ((Data['Stim'][RStiffnessWindow]).max() - (Data['Stim'][RForceWindow]).mean())/(Data['Stim'][RForceWindow].mean())
        #RStiffness = RdF/RdLo
        #RModulus = (RdF/CSA)/RStrain

        #return RStiffness, RStrain, RModulus

    def ktr_model(x, a, kt, c):
        return a * (1-np.exp(-kt*x)) + c

    def generate_Initial_Parameters(x_data: pd.Series = None, y_data: pd.Series = None):
        # min and max used for bounds
        def sumOfSquaredError(parameterTuple):
            # do not print warnings by genetic algorithm
            warnings.filterwarnings("ignore")
            val = ktr_model(x_data, *parameterTuple)
            return(np.sum((y_data - val) ** 2.0))

        maxX = max(x_data)
        minX = min(x_data)
        maxY = max(y_data)
        minY = min(y_data)
        
        # find max force during ktr window. To account for files where force drops substantially
        # we compare whether the last 500 ms force is lower than 90% of peak force
        # if it is, we use 90% of peak force as our new maximal force to fix curve to
        Max_Force_param: float = maxY - minY
        # Max_Force_param: float = maxY 

        # Force at ktr start
        Force_at_T0: float = y_data[0] 
        
        parameterBounds = []
        # search bounds for a (force when at plateau)
        parameterBounds.append([maxY, maxY])
        # search bounds for kt (range of values software uses to find ktr)
        parameterBounds.append([0, 30])
        # searh bounds for c (force at t=0)
        parameterBounds.append([Force_at_T0, Force_at_T0])

        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
        
        return result.x

    Stiffness, Strain, Modulus = StiffnessAnalysis(Data = Data)
    RStiffness, RStrain, RModulus = StiffnessAnalysis(Data = Data, stiffness_time = relax_stiff_time)
    Peak_force = Data['Force'][200000:210000].mean()
    
    # Fit ktr only to subset of data
    # ktr_start_time = following restretch (i.e., 10160)
    # ktr_end_time = 60000
    modelData = pd.DataFrame(Data[['Time', 'Force']][ktr_start_time:ktr_end_time]).reset_index(drop = False)

    # Find min force value after restretch occurs
    # Becomes real start to modelData
    min_force_index = modelData[['Force']].idxmin()
    def ktr_fitting(offset):
        ktr_start: int = min_force_index[0]+(offset)
        #print(offset)
        # Cutoff end of raw data to avoid any potential bath moves/movements at end of test
        # Would negatively affect curve fitting
        ktr_end: int = ktr_start + 30000

        # Put time and force data into numpy arrays for curve fitting
        x_data = np.array(modelData['Time'][ktr_start:ktr_end])
        time_offset = x_data[0]

        norm_x_data = np.subtract(x_data, time_offset)
        x_data = norm_x_data
        y_data = np.array(modelData['Force'][ktr_start:ktr_end])
        # Find initial parameters for curve fitting
        ktr_Parameters = generate_Initial_Parameters(x_data, y_data)

        # maxfev = number of iterations code will attempt to find optimal curve fit
        maxfev:int = 1000
        try:
            fittedParameters, pcov = curve_fit(ktr_model, x_data, y_data, ktr_Parameters, maxfev=maxfev)
        except:
            try:
                maxfev = 5000
                fittedParameters, pcov = curve_fit(ktr_model, x_data, y_data, ktr_Parameters, maxfev=maxfev)
            except:
                print(Error(f"ktr parameters were not fit after {maxfev} iterations for file: {Filename}. Added to 'Files to Check'"))
                return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        # Generate model predictions

        modelPredictions = ktr_model(x_data, *fittedParameters)
        # Calculate error
        Max_Force_error: float = np.sqrt(pcov[0, 0])
        ktr_error: float = np.sqrt(pcov[1, 1])
        Force_at_T0_error: float = np.sqrt(pcov[2, 2])
        ktr: float = fittedParameters[1]
        absError = modelPredictions - y_data
        SE: float = np.square(absError)  # squared errors
        MSE: float = np.mean(SE)  # mean squared errors
        RMSE: float = np.sqrt(MSE)  # Root Mean Squared Error, RMSE
        GoodnessFit: float = 1.0 - (np.var(absError) / np.var(y_data))
        return ktr, GoodnessFit, time_offset, x_data, fittedParameters
    
    def handle_dem_mouse_clickies(): #Code to define interactable graph
        click_x=round(coords[0][0]*10000)
        #print (click_x-250000)
        click_y=round(coords[0][1])

        ax.plot(Data['Time'][click_x-250000:ktr_end_time], Data['Force'][click_x-250000:ktr_end_time], color='red') #plot graph data from the pre-selected model data of the KTR time
        ax.axvline(x=Data['Time'][click_x-250000], ymin = 0, ymax = 1, color = 'Green', linestyle = '-', linewidth = 2) #show cursor line
        fig.canvas.draw()
        plt.show()

        return click_x

    ktr, GoodnessFit, time_offset, x_data, fittedParameters=ktr_fitting(80) #add a print
    if ktr <90: #if KTR is below ** then it runs this section to allow you to manually adjust the data pull for it
        fig=plt.figure()
        ax=fig.add_subplot()

        ax.plot(modelData['Time'], modelData['Force'], color='black')
        coords=plt.ginput(n=1, show_clicks=True, timeout=9999) #plots the data and then allows you to click a point
        if coords:
            ktr_delay=handle_dem_mouse_clickies() #given a mouse click, we pull the function from above for the interactable graph to get a defined start time to the KTR data
            ktr, GoodnessFit, time_offset, x_data, fittedParameters=ktr_fitting(ktr_delay-560160) #uses the manually selected ktr time to run the KTR fitting
            #add a print
    
    Xmodel = np.linspace(min(x_data), max(x_data), 100)
    Ymodel = ktr_model(Xmodel, *fittedParameters)

    ktrForce = Data['Force'][-5000:].mean()

    if Graph == True:
        fig = plt.figure()
        
        plt.plot(Data['Time'], Data['Force'], color = 'black', label = 'Raw')
        plt.plot(Data['Time'][200000:210000], Data['Force'][200000:210000], color = 'blue', label = 'Peak')
        plt.plot(Data['Time'][1:500], Data['Force'][1:500], color = 'Yellow', label = 'Baseline')
        plt.plot(Xmodel + time_offset, Ymodel, color = 'red', label = 'Fit')
        plt.ylabel('Force (mN)')
        plt.xlabel('Time (s)')
        plt.text(
            x = 0.5, y = 0.1,
            s = f'ktr = {ktr:.3f}\n'
                f'Goodness of fit = {GoodnessFit * 100:.2f}%\n' 
                f'Active Specific Force = {(Peak_force/CSA):.2f}',
            transform = plt.gca().transAxes,
            horizontalalignment = 'center',
            verticalalignment = 'center')
        plt.title(Filename)
        plt.legend()
        plt.show()

    return Stiffness, Strain, Modulus, ktr, GoodnessFit, Xmodel, Ymodel, ktrForce, Peak_force, RStiffness, RStrain, RModulus, ktr_delay

def main():
    #Directory = filedialog.askopenfilenames()    
    Directory = '/Users/keaton.briar/Desktop/Single Fibre Data/Puncture Study/S07/LR'

    AllFiles = []
    for root, subdirs, files in os.walk(Directory):
        for file in files:
            if 'Store' in file:
                continue
            AllFiles.append(os.path.join(root, file))
    AllFiles = sorted(AllFiles)


    final_results = []
    for file in AllFiles:

        Data, Subject, Muscle, Fibre, Level, FibreLength, SarcomereLength, CSA = ReadFile(File = file)
        Stiffness, Strain, Modulus, ktr, GoodnessFit, Xmodel, Ymodel, ktrForce, Peak_force, RStiffness, RStrain, RModulus, ktr_delay = ktrAnalysis(Data = Data, Filename = os.path.basename(file), CSA=CSA, Graph = True)

        OrganizedData = {
            'Subject': Subject,
            'Muscle': Muscle,
            'Fibre Length': FibreLength,
            'Fibre': Fibre,
            'Level': Level,
            'Sarcomere Length': SarcomereLength,
            'CSA': CSA,
            'Stiffness (pCa 4.5)': Stiffness,
            'Stiffness in Relaxing': RStiffness,
            'Modulus (pCa 4.5)': Modulus,
            'Modulus in Relaxing': RModulus,
            'Normalized Modulus': Modulus - RModulus, 
            'ktr (pCa 4.5)': ktr,
            'Goodness of fit (pCa 4.5)': GoodnessFit,
            'Active Specific Force': Peak_force / CSA,
            'KTR Delay': ktr_delay - ktr_start_time,
            'KTR Delay True': ((ktr_delay - ktr_start_time)-250000)/10}
            #'ktr Force (pCa 4.5)': ktrForce,
            #'ktr Specific Force (pCa 4.5)': ktrForce / CSA}
        
        final_results.append(OrganizedData)
    
    for test in final_results:
        print(test)

    # return Final_results
    return final_results 

#call the main () functio to get final_results
final_results = main()

# your code for generating the final results
for File in final_results:
     excel_results = pd.DataFrame(columns = [key for key in File.keys()])

for idx, results in enumerate(final_results):
    for column, data in results.items():
        excel_results.loc[idx, column] = data

excel_file_path = '/Users/keaton.briar/Desktop/Single Fibre Data/Puncture Study/S07/SO7LR.xlsx'

excel_results.to_excel(excel_file_path, index = False)


if __name__ == "__main__":
    main()