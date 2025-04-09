import pandas as pd

# Load the CSV file
file_path = "I0_10.csv"
df = pd.read_csv(file_path)

# Drop the first row which contains repeated 'X' and 'Y' labels
df_clean = df.drop(index=0).reset_index(drop=True)

# Define function to get x and y values for a given Sigma
def get_xy_for_sigma(df, sigma_label):
    # Find the column index for the specified Sigma
    col_index = df.columns.get_loc(sigma_label)
    
    # Get the X and Y column names
    x_col = df.columns[col_index]
    y_col = df.columns[col_index + 1]
    
    # Extract and convert values to float
    x_values = df[x_col].astype(float).reset_index(drop=True)
    y_values = df[y_col].astype(float).reset_index(drop=True)
    
    return x_values, y_values

# Example usage
sigma = "Sigma1.2"
x_vals, y_vals = get_xy_for_sigma(df_clean, sigma)

# Display first few values
print(x_vals.head())
print(y_vals.head())
