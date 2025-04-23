import numpy as np
import pandas as pd


def calculate_incidence_deflection(Beta1,Beta2, Sigma, toc_max, foil_type): #Beta 1, solidity c/S, max thickness over chord ratio, type of foil (DCA/NACA-65)
    """Calculates the incidence and deflection based on the velocity triangles and blade foil shape.
    """

    
    # Load data
    I0_data = pd.read_csv('Blade_Angle_Data/I0_10.csv')
    n_data = pd.read_csv('Blade_Angle_Data/n-coefficient.csv')
    Ki_data = pd.read_csv('Blade_Angle_Data/Kit-coefficient.csv')
    Delta0_data = pd.read_csv('Blade_Angle_Data/Delta0_10.csv')
    Kd_data = pd.read_csv('Blade_Angle_Data/Kdeltat-coefficient.csv')
    m_data = pd.read_csv('Blade_Angle_Data/m-coefficient.csv')
    b_data = pd.read_csv('Blade_Angle_Data/b-coefficient.csv')
    
    I010 = interpolate_Multicurve(Beta1, Sigma, I0_data)
    n = interpolate_Multicurve(Beta1, Sigma, n_data)
    Ki_t = interpolate_single_curve(toc_max, Ki_data)
    Kd_t = interpolate_single_curve(toc_max, Kd_data)
    Delta010 = interpolate_Multicurve(Beta1, Sigma, n_data)
    m = interpolate_single_curve(Beta1, m_data)
    b = interpolate_single_curve(Beta1, b_data)


    if foil_type == "DCA":
        Ki_sh = 0.7
        Kd_sh = 0.75
        
    elif foil_type == "NACA-65":
        Ki_sh = 1.1
        Kd_sh = 1.1
    else:
        raise ValueError("No foil_type was provided! Function ::: calculate_incidence_deflection(Beta1, Sigma, toc_max, foil_type) where foil_type is 'DCA' or 'NACA-65' ")
    I0 = Ki_sh * Ki_t * I010
    Delta0 = Kd_sh * Kd_t * Delta010
    theta = ( Beta2 - Beta1 + Delta0 - I0)/(1 - (m/Sigma**b)+n)
    I = I0 + n*theta
    Delta = Delta0 + (m/Sigma**b)*theta
    return(I, Delta)
    
    

def interpolate_Multicurve(x, Sigma, file):
    """
    Interpolates y-value for a given Beta1 and Sigma.
    'file' is a DataFrame where each column is named 'Sigma{value}' and rows are Beta1 values.
    """


    sigma_columns = [col for col in file.columns if col.startswith("Sigma")]

    # Create mapping of column names and sigma values
    sigma_values = [float(col.replace("Sigma", "")) for col in sigma_columns]
    sigma_map = dict(zip(sigma_values, sigma_columns)) #map is used to go from sigma to relevant column name. This will then be used to call the relevant columns


    sigma_sorted = sorted(sigma_map.keys()) #list of all sigma values available

    # Check for exact match
    if Sigma in sigma_sorted:
        ID = sigma_map[Sigma]#name of column corresponding to value Sigma
        x_val, y_val = get_xy_for_label(file,ID)
        return np.interp(x, x_val, y_val)

    # Find bounding sigmas
    lower_idx = np.searchsorted(sigma_sorted, Sigma) - 1
    upper_idx = lower_idx + 1
    SigmaL = sigma_sorted[lower_idx]
    SigmaU = sigma_sorted[upper_idx]


    if lower_idx < 0 or upper_idx >= len(sigma_sorted):
        raise ValueError("Sigma value out of interpolation range.")
    else:
            #lower calculations for interpolation
            IDL = sigma_map[SigmaL]
            x_valL, y_valL = get_xy_for_label(file,IDL)
            yL = np.interp(x, x_valL, y_valL)

            #upper calculation
            IDU = sigma_map[SigmaU]
            x_valU, y_valU = get_xy_for_label(file,IDU)
            yU = np.interp(x, x_valU, y_valU)

            #interpolation
            y = yL + (Sigma-SigmaL)*(yU-yL)/(SigmaU-SigmaL)
            return y


def interpolate_single_curve(x, file):
    # can just find the first label then run get_xy_for_label(df,label) then interpolate
    All_columns = [col for col in file.columns]#create list of labels
    if len(All_columns)==2:
        x_val, y_val = get_xy_for_label(file,All_columns[0])#get first and only name
        return np.interp(x, x_val, y_val)
    else:
        raise ValueError("You are asking for single interpolation of a file with multiple curves!")

# Function to get list of x and y for sigma
def get_xy_for_label(df,label):
    # Find the column index for the specified Sigma
    dfn = df.drop(index=0).reset_index(drop=True)
    col_index = dfn.columns.get_loc(label)

    #print(col_index)
    # Get the X and Y column names
    x_col = dfn.columns[col_index]
    y_col = dfn.columns[col_index + 1]

    # Extract and convert values to float
    x_values = list(dfn[x_col].astype(float).reset_index(drop=True))
    y_values = list(dfn[y_col].astype(float).reset_index(drop=True))

    return x_values, y_values


#Run to test
print(calculate_incidence_deflection(Beta1=50,Beta2=70, Sigma=1.01, toc_max=0.05, foil_type='DCA'))

