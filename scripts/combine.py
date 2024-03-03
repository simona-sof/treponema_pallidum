import pandas as pd

# COMBINE THE THREE LAT-LONG files
# Read the three CSV files of lat_long data from each metadata into DataFrames
df1 = pd.read_csv("config/lat_longs_Nichols.csv", header=None)
df2 = pd.read_csv("config/lat_longs_SS14.csv", header=None)
df3 = pd.read_csv("config/lat_longs_TPE.csv", header=None)

# Concatenate the DataFrames vertically
combined_df = pd.concat([df1, df2, df3], ignore_index=True)

# Drop duplicates based on latitude and longitude columns
combined_df.drop_duplicates()

# Write the combined DataFrame to a new CSV file
combined_df[3] = 'location'
desired_order = [3, 0, 1, 2]
combined_df = combined_df.iloc[:, desired_order]
combined_df.to_csv("config/lat_longs.tsv", index=False, sep="\t", header=False)