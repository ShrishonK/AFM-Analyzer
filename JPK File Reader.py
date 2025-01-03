from pyfmreader import loadfile
from pyfmrheo.routines.TingFit import doTingFit
import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy import integrate
import os

numFiles = 0
ContactPoint = []
TurnaroundPoint = []
FileName = []
AreaCurve = []

# This function makes sure that it goes through a main folder and looks through each sub-folder and reads the data
def process_folder_and_files(base_folder):
    if not os.path.isdir(base_folder):
        print(f"The path {base_folder} is not a valid directory.")
        return
    
    all_curves = []

    # Iterate through subfolders and files
    for root, dirs, files in os.walk(base_folder):
        folder_curves = []
        for file in files:
            if not file.endswith(".jpk-force-map"):
                continue
            
            file_path = os.path.join(root, file)
            curve_data = mainFunc(file_path)
            folder_curves.append(curve_data)
            all_curves.append(curve_data)

        if folder_curves:  # If we have data for this folder
            plt.figure(figsize=(10, 5))
            
            # Plot all individual curves in black
            for curve in folder_curves:
                plt.plot(curve['x'], curve['y'], color='black', alpha=0.1, linewidth=1)
            
            # Find the length of each curve
            curve_lengths = [len(curve['x']) for curve in folder_curves]
            
            # Find the minimum length among all curves
            min_length = min(curve_lengths)
            
            # Truncate all curves to the minimum length
            x_values = np.array([curve['x'][:min_length] for curve in folder_curves])
            y_values = np.array([curve['y'][:min_length] for curve in folder_curves])
            
            # Calculate average x and y values
            avg_x = np.mean(x_values, axis=0)
            avg_y = np.mean(y_values, axis=0)
            
            # Plot average curve
            plt.plot(avg_x, avg_y, 'r-', linewidth=3, label='Average')
            
            plt.xlabel('Indentation [nm]', fontsize=15)
            plt.ylabel('Force [pN]', fontsize=15)
            plt.title(f'Force Curves for {os.path.basename(root)}')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(root, 'combined_plot.png'))
            plt.close()

def mainFunc(filePath):
    # Define global parameters for plotting
    plt.rcParams["figure.figsize"] = (10,5)

    # Define path of file to process
    file_path = filePath
    new_path = os.path.dirname(file_path)
    filename = file_path.split('\\')[-1]
    FileName.append(filename)

    # Load File
    file = loadfile(file_path)
    filemetadata = file.filemetadata
    print(filemetadata['file_type'])

    # Get some of the file metadata
    closed_loop = filemetadata['z_closed_loop']
    file_deflection_sensitivity = filemetadata['defl_sens_nmbyV'] #nm/V
    file_spring_constant = filemetadata['spring_const_Nbym'] #N/m
    height_channel = filemetadata['height_channel_key']

    deflection_sensitivity = file_deflection_sensitivity / 1e9 #m/V
    spring_constant = file_spring_constant

    print(f"Closed loop: {closed_loop}")
    print(f"Height channel: {height_channel}")
    print(f"Deflection Sens.: {deflection_sensitivity} m/V")
    print(f"Spring Constant: {spring_constant} N/m")

    # Define parameters to perform the HertzFit
    maxnoncontact = 1e-6 #um
    pts_downsample = 300

    param_dict = {
        'height_channel': height_channel,   # Channel where to find the height data
        'def_sens': deflection_sensitivity, # Deflection sensitivity in m/V
        'k': spring_constant,               # Spring constant in N/m
        'contact_model': 'paraboloid',      # Geometry of the indenter: paraboloidal, conical, pyramidal
        'tip_param': 5e-06,                 # Tip raidus in meters or tip angle in degrees
        'curve_seg': 'extend',              # Segement to perform the fit: extend or retract
        'correct_tilt': False,              # Perform tilt correction
        'tilt_min_offset': 1e-08,           # Maximum range where to perform the tilt correction in meters
        'tilt_max_offset': 1e-06,           # Minimum range where to perform the tilt correction in meters
        'poisson': 0.5,                     # Poisson's ratio
        'poc_method': 'RoV',                # Method to find the contact point: RoV or RegulaFalsi
        'poc_win': 4e-07,                   # Window size for the RoV method
        'max_ind': 0.0,                     # Maximum indentation range for fit in meters
        'min_ind': 0.0,                     # Minimum indentation range for fit in meters
        'max_force': 0.0,                   # Maximum force range for fit in Newtons
        'min_force': 0.0,                   # Minimum force range for fit in Newtons
        'fit_range_type': 'full',           # Fit data range: full, indentation or force
        'vdragcorr': False,                 # Compute viscous drag correction from baseline separation.
        'polyordr': 2,                      # Order of Polynomial to fit to the baselines to compute viscous drag
        'rampspeed': 0.0,                   # Ramp speed in m/s
        'compute_v_flag': False,            # Compute ramp speed
        't0': 1,                            # Scaling t0
        'd0': 0.0,                          # Initial point of contact
        'slope': 0.0,                       # Initial slope
        'auto_init_E0': True,               # Estimate automatically the initial value of the Young's Modulus
        'E0': 1000,                         # Initial Young's Modulus value
        'tc': 0.0,                          # Initial contact time
        'auto_init_betaE': True,            # Estimate automatically the fluidity exponent beta
        'fluid_exp': 0.2,                   # Initial fluidity exponent
        'f0': 0.0,                          # Initial F0 value
        'vdrag': 2.5e-06,                   # Viscous drag value in N/ms
        'model_type': 'analytical',         # Viscoelastic model type
        'smoothing_win': 5,                 # Number of points to use for smoothing the signal
        'contact_offset': maxnoncontact,    # Baseline offset for the Hertz Fit
        'fit_line': False,                  # Fit line to the baseline
        'downsample_flag': True,            # Downsample the signal for Hertz Fit
        'pts_downsample': pts_downsample,   # Number of points to downsample
        'offset_type':'percentage',         # How to correct for baseline offset: percentage or value
        'max_offset':.3,                    # Max percentage to compute offset
        'min_offset':0                      # Min percentage to compute offset
    }



    # Select curve by index
    curve_idx = 0
    force_curve = file.getcurve(curve_idx)

    # Preprocess curve
    force_curve.preprocess_force_curve(param_dict['def_sens'], param_dict['height_channel'])

    # JPK files require the height signal to be shifted
    if filemetadata['file_type'] in ('jpk-force', 'jpk-force-map', 'jpk-qi-data'):
        force_curve.shift_height()

    # Run fit
    ting_result, hertz_result = doTingFit(force_curve, param_dict)

    # Check Hertz Result values
    #hertz_result.fit_report()
    # Check Ting Result values
    #ting_result.fit_report()



    # Prepare data for plotting
    segs = force_curve.get_segments()
    ext_data = segs[0][1]
    ret_data = segs[1][1]
    idx_tc = (np.abs(ext_data.indentation - 0)).argmin()
    t0 = ext_data.time[-1]
    indentation = np.r_[ext_data.indentation, ret_data.indentation]
    time = np.r_[ext_data.time, ret_data.time + t0]
    force = np.r_[ext_data.force, ret_data.force]
    fit_mask = indentation > (-1 * maxnoncontact)
    tc = time[idx_tc]
    ind_fit = indentation[fit_mask]
    force_fit = force[fit_mask]
    force_fit = force_fit - force_fit[0]
    time_fit = time[fit_mask]
    tc_fit = tc-time_fit[0]
    time_fit = time_fit - time_fit[0] - tc_fit

    downfactor= len(time_fit) // pts_downsample
    idxDown = list(range(0, len(time_fit), downfactor))

    idx_tc = (np.abs(time_fit[idxDown] - ting_result.tc)).argmin()
    d0 = ind_fit[idxDown][idx_tc]

    # Need to find the Turnaround point
    maxVal = max((ind_fit[idxDown] - d0) * 1e9)
    print(f"Turnaround Point: {maxVal}")
    TurnaroundPoint.append(maxVal)

    
    
    areaArrayX = []
    areaArrayY = []
    areaIndex = []
    count = 0

    # Filter out all values that don't lie within the x ranges
    for i in ind_fit[idxDown]:
        if ((i - d0)* 1e9) >= hertz_result.delta0 and ((i - d0)* 1e9) <= maxVal:
            areaArrayX.append((i - d0)* 1e9)
            areaIndex.append(count)
        count+=1

    for j in areaIndex:
        areaArrayY.append(force_fit[idxDown][j] * 1e12)

    ContactPoint.append(hertz_result.delta0)
    x = areaArrayX
    y = areaArrayY

    goingCurveX = []
    goingCurveY = []
    returningCurveX = []
    returningCurveY = []
    indexVal = 0

    for k in range(0, len(x)):
        if x[k] == maxVal:
            goingCurveX.append(x[k])
            goingCurveY.append(y[k])
            indexVal = k
            break
        else:
            goingCurveX.append(x[k])
            goingCurveY.append(y[k])

    for m in range(indexVal, len(x)):
        returningCurveX.append(x[m])
        returningCurveY.append(y[m])
    
    # Finds the area
    area = integrate.trapezoid(y, x)  # Use trapezoidal rule for integration
    print(f"Enclosed area under the curve: {area}")
    AreaCurve.append(area)


    coefficientsGoing = np.polyfit(goingCurveX, goingCurveY, 1)
    coeffiicentsReturning = np.polyfit(returningCurveX, returningCurveY, 1)
    polynomialGoing = np.poly1d(coefficientsGoing)
    polynomialReturning = np.poly1d(coeffiicentsReturning)

    plt.clf()

    # Create points for the line
    x_line = np.linspace(min(goingCurveX), maxVal, 100)
    y_line = polynomialGoing(x_line)



    # Plot
    plt.plot(
        (ind_fit[idxDown] - d0) * 1e9,
        force_fit[idxDown] * 1e12,
        color='black',
        linewidth=3,
        label='Experimental Data'
    )

    plt.plot(
        goingCurveX,
        goingCurveY,
        color='#599cff',
        linewidth=3,
        label='Going Curve'
    )
    plt.plot(x_line, y_line, color='red', label=f'Line of best fit: y = {coefficientsGoing[0]:.2f}x + {coefficientsGoing[1]:.2f}')

    plt.plot(
        returningCurveX,
        returningCurveY,
        color='#ed6b7a',
        linewidth=3,
        label='Returning Curve'
    )

    plt.fill_between(x, y, 0, color='lightblue', label='Area') 


    

    plt.axvline(x=hertz_result.delta0, color="green", linestyle="--", linewidth=2, label="Contact Point")
    plt.axvline(x=maxVal, color="green", linestyle="--", linewidth=2, label="Turn Point")
    plt.xlabel('Indentation [nm]', fontsize = 15)
    plt.ylabel('Force [pN]', fontsize = 15)
    plt.legend()
    plt.savefig(f'{new_path}\\{filename}.png')
    #plt.show() # Comment out this line to make the CSV file immediately without having to close each plot
    plt.close()
    
    return {
        'x': (ind_fit[idxDown] - d0) * 1e9,
        'y': force_fit[idxDown] * 1e12,
        'going_x': goingCurveX,
        'going_y': goingCurveY,
        'returning_x': returningCurveX,
        'returning_y': returningCurveY,
        'contact_point': hertz_result.delta0,
        'turn_point': maxVal,
        'filename': filename,
        'new_path': new_path
    }
    

process_folder_and_files("Enter file path here")

# Replace the file path below with your local file path
with open("Enter file path here\\samplename.csv", 'w', newline = '') as csvfile:
    my_writer = csv.writer(csvfile)
    my_writer.writerow(['Filename', 'Contact Point', 'Turnaround Point', 'Area Curve'])
    for i in range(len(FileName)):
        my_writer.writerow([FileName[i], ContactPoint[i], TurnaroundPoint[i], AreaCurve[i]])
