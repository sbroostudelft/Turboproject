## Old stuff
import numpy as np
import pandas as pd


def calculate_incidence_deflection(Beta1, Sigma, toc_max, foil_type): #Beta 1, solidity c/S, max thickness over chord ratio, type of foil (DCA/NACA-65)
    """Calculates the incidence and deflection based on the velocity triangles and blade foil shape.
    """

    
    # Load data
    I0_data = pd.read_csv('I0_10.csv')
    n_data = pd.read_csv('n-coefficient.csv')
    Ki_data = pd.read_csv('Kit-coefficient.csv')
    Delta0_data = pd.read_csv('Delta0_10.csv')
    Kd_data = pd.read_csv('Kdeltat-coefficient.csv')
    m_data = pd.read_csv('m-coefficient.csv')
    b_data = pd.read_csv('b-coefficient.csv')

    I0 = interpolate_Multicurve(Beta1, Sigma, I0_data)
    return I0

def interpolate_Multicurve(Beta1, Sigma, file):
    """
    Interpolates y-value for a given Beta1 and Sigma.
    'file' is a DataFrame where each column is named 'Sigma{value}' and rows are Beta1 values.
    """

    # Only keep columns with 'Sigma' in the name
    sigma_columns = [col for col in file.columns if col.startswith("Sigma")]
    Y_columns = [col for col in file.columns if col.startswith("Unnamed")]
    
    # Create mapping of column names and sigma values
    sigma_values = [float(col.replace("Sigma", "")) for col in sigma_columns]
    sigma_map = dict(zip(sigma_values, sigma_columns)) #map is used to go from sigma to relevant column name. This will then be used to call the relevant columns
    print(sigma_map , "this was sigma map")

    # Sort the sigmas
    sigma_sorted = sorted(sigma_map.keys())

    # Check for exact match
    if Sigma in sigma_sorted:
        col = sigma_map[Sigma]
        print(col)
        print(file[col].values)
        print(file.index.values)
        return np.interp(Beta1, file.index.values, file[col].values)

    # Find bounding sigmas
    lower_idx = np.searchsorted(sigma_sorted, Sigma) - 1
    upper_idx = lower_idx + 1

    if lower_idx < 0 or upper_idx >= len(sigma_sorted):
        raise ValueError("Sigma value out of interpolation range.")



    sigma_low = sigma_sorted[lower_idx]
    sigma_high = sigma_sorted[upper_idx]

    col_low = sigma_map[sigma_low]
    col_high = sigma_map[sigma_high]
    print("I'm still standing!")

    print(file.index.values)
    print(file[col_low].values)
    print(file[col_high].values)
    # Interpolate for Beta1 at both bounding sigmas
    y_low = np.interp(Beta1, file.index.values, file[col_low].values)
    y_high = np.interp(Beta1, file.index.values, file[col_high].values)

    # Interpolate between the two y-values
    y_interp = np.interp(Sigma, [sigma_low, sigma_high], [y_low, y_high])

    return y_interp


#Run to test
print(calculate_incidence_deflection(Beta1=60, Sigma=1.00, toc_max=0, foil_type=0))
