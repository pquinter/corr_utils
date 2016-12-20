"""
Pivot a table
Takes three arguments from stdin
file_dir: path to file to pivot, must be comma separated with header
col1: name of column to be used as index
col2: name of column to be used as header
values: name of column with values of intersection col1,col2

Run from terminal as:
python pivot.py file_dir col1 col2 values

Writes pivotted table to pivoted_output.txt

Author: Porfirio Quintero-Cadena
"""

import sys
import pandas as pd

# Get file and columns to pivot
file_dir, col1, col2, values = sys.argv[1:]

# Read DataFrame and pivot and write to file
df = pd.read_csv(file_dir)
pivoted = df.pivot(index=col1, columns=col2, values=values)
pivoted.to_csv('pivoted_output.txt', index=True)
print('Output written to pivoted_output.txt')
