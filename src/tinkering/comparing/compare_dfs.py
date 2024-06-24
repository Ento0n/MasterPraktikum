import pandas as pd

path_df1 = "output/domains_2_chunks"
path_df2 = "output/domains_1_chunks"

df1 = pd.read_csv(path_df1)
df2 = pd.read_csv(path_df2)

comparison = df1.compare(df2)
print(comparison)
print(len(comparison))
