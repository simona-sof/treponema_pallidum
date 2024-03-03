import pandas as pd

def lineage(df):
    return df.loc[df['ENA_Accession']!='-']

def getDate(df):
    date = df['Sample_Year']
    date = date.apply(lambda x: f"{x}-XX-XX")
    return date

def separate_meta(df):
    new_df = pd.DataFrame()
    new_df['strain'] = df['Sample']
    new_df['name'] = 'TPE'
    new_df['accession'] = df['ENA_Accession']
    new_df['date'] = getDate(df)
    # new_df['region'] = df['Continent']
    new_df['location'] = df['Country'].apply(lambda x: x.lower() if isinstance(x, str) else x)
    new_df['authors'] = df['Study']
    new_df['authors'] = new_df['authors'].replace('This_Study', 'Beale_2020')
    new_df['sample_type'] = 'Clinical_Swab'
    new_df['reads_or_assemblies'] = df['Reads_or_sim-reads']
    return new_df

def main():
    df = pd.read_excel("pre-processing/beale_2020.xlsx", header=0)
    SS14_df = lineage(df)
    TPE = separate_meta(SS14_df)
    print(TPE)
    TPE.to_csv('pre-processing/meta/meta_beale_2020_TPE.csv', index=False)

if __name__ == "__main__":
    main()