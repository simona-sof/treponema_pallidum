import pandas as pd

def lineage(df, lineage):
    return df.loc[(df['Lineage']==lineage) & (df['SRA.ENA_Accession']!='-')]

def getDate(df):
    date = df['Sample_Date']
    date = date.apply(lambda x: f"{x}-XX-XX")
    return date

def separate_meta(df):
    new_df = pd.DataFrame()
    new_df['strain'] = df['SampleShort']
    new_df['name'] = df['Lineage']
    new_df['accession'] = df['SRA.ENA_Accession']
    new_df['date'] = getDate(df)
    new_df['region'] = df['Continent']
    new_df['location'] = df['GeoCountry'].apply(lambda x: x.lower() if isinstance(x, str) else x)
    new_df['authors'] = df['Citation']
    new_df['authors'] = new_df['authors'].replace('This_Study', 'Beale_2019')
    new_df['sample_type'] = df['SampleType']
    new_df['reads_or_assemblies'] = df['Analysed_reads_or_assembly']
    return new_df

def main():
    df = pd.read_excel("pre-processing/beale_2019.xlsx", header=1)
    SS14_df = lineage(df, 'SS14')
    Nichols_df = lineage(df, 'Nichols')
    # SS14_df['SRA.ENA_Accession'].to_csv('accessions/SS14_beale_2019_accessions.txt',index=False, header=False)
    # Nichols_df['SRA.ENA_Accession'].to_csv('accessions/Nichols_beale_2019_accessions.txt',index=False, header=False)
    SS14 = separate_meta(SS14_df)
    Nichols = separate_meta(Nichols_df)
    SS14.to_csv('pre-processing/meta/meta_beale_2019_SS14.csv', index=False)
    Nichols.to_csv('pre-processing/meta/meta_beale_2019_Nichols.csv', index=False)

if __name__ == "__main__":
    main()