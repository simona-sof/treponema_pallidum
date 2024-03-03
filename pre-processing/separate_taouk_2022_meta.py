import pandas as pd

def lineage(df, lineage):
    return df.loc[(df['Major_lineage']==lineage) & (df['ENA/SSR_Accession']!='-')]

def getDate(df):
    date = df['Collection_date']
    return date

def separate_meta(df):
    new_df = pd.DataFrame()
    new_df['strain'] = df['Sample_ID']
    new_df['name'] = df['Major_lineage']
    new_df['accession'] = df['ENA/SSR_Accession']
    new_df['date'] = getDate(df)
    new_df['region'] = df['Continent']
    new_df['location'] = df['Geo_country'].apply(lambda x: x.lower() if isinstance(x, str) else x)
    new_df['authors'] = df['Study']
    new_df['authors'] = new_df['authors'].replace('This study', 'Taouk_2022')
    new_df['sample_type'] = df['Sample_type']
    new_df['reads_or_assemblies'] = df['Reads_or_assembly'].apply(lambda x: x.lower() if isinstance(x, str) else x)
    return new_df

def main():
    df = pd.read_excel("pre-processing/taouk_2022.xlsx", header=3)
    SS14_df = lineage(df, 'SS14')
    Nichols_df = lineage(df, 'Nichols')

    SS14 = separate_meta(SS14_df)
    Nichols = separate_meta(Nichols_df)
    SS14.to_csv('pre-processing/meta/meta_taouk_2022_SS14.csv', index=False)
    Nichols.to_csv('pre-processing/meta/meta_taouk_2022_Nichols.csv', index=False)

if __name__ == "__main__":
    main()