import pandas as pd

def lineage(df, lineage):
    return df.loc[(df['TPA_Lineage']==lineage) & (df['SRR/ENA_Accession']!='-')]

def getDate(df):
    date = df['Sample_Year']
    date = date.apply(lambda x: f"{x}-XX-XX")
    return date

def separate_meta(df):
    new_df = pd.DataFrame()
    new_df['strain'] = df['Sample_Name']
    new_df['name'] = df['TPA_Lineage']
    new_df['accession'] = df['SRR/ENA_Accession']
    new_df['date'] = getDate(df)
    new_df['region'] = df['Continent']
    new_df['location'] = df['Geo_Country'].apply(lambda x: x.lower() if isinstance(x, str) else x)
    new_df['authors'] = df['Citation']
    new_df['authors'] = new_df['authors'].replace('This_Study', 'Beale_2021')
    new_df['sample_type'] = df['Sample_Type']
    new_df['reads_or_assemblies'] = df['Reads_or_assemblies']
    return new_df

def main():
    df = pd.read_excel("pre-processing/beale_2021.xlsx", header=0)
    SS14_df = lineage(df, 'SS14')
    Nichols_df = lineage(df, 'Nichols')

    # SS14_df['SRA.ENA_Accession'].to_csv('accessions/SS14_beale_2019_accessions.txt',index=False, header=False)
    # Nichols_df['SRA.ENA_Accession'].to_csv('accessions/Nichols_beale_2019_accessions.txt',index=False, header=False)
    SS14 = separate_meta(SS14_df)
    Nichols = separate_meta(Nichols_df)
    SS14.to_csv('pre-processing/meta/meta_beale_2021_SS14.csv', index=False)
    Nichols.to_csv('pre-processing/meta/meta_beale_2021_Nichols.csv', index=False)

if __name__ == "__main__":
    main()