import sys
import pandas as pd

# Get file and columns to pivot
file_dir, col1, col2, values = sys.argv[1:]

# Read DataFrame and pivot and write to file
df = pd.read_csv(file_dir)
pivoted = df.pivot(index=col1, columns=col2, values=values)
pivoted.to_csv('pivoted_output.txt', index=True)
print('Output written to pivoted_output.txt')
