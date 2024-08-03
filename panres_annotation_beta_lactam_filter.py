import pandas as pd
import sys

# Read the input file path and output file path from Snakemake
input_file = sys.argv[1]
output_file = sys.argv[2]

# Read the TSV file
data = pd.read_csv(input_file, sep='\t', comment='#')

# Filter and process the data
filtered_data = (data[data['variable'] == 'class']
                 .groupby('gene')
                 .filter(lambda x: len(x['value'].unique()) == 1 and x['value'].iloc[0] == 'beta_lactam')
                 .drop_duplicates(subset=['gene'])
                 .loc[:, ['gene']])
  
# Ensure all genes have the version ending "_v1.0.1"  
filtered_data['gene'] = filtered_data['gene'].apply(lambda x: x if x.endswith('_v1.0.1') else x + '_v1.0.1')


# Save the filtered data to the output file
filtered_data.to_csv(output_file, sep='\t', index=False, header=False)
