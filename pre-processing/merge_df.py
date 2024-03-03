import pandas as pd

def main():
    df_SS14_2019 = pd.read_csv("pre-processing/meta/meta_beale_2019_SS14.csv")
    df_SS14_2021 = pd.read_csv("pre-processing/meta/meta_beale_2021_SS14.csv")
    df_SS14_2022 = pd.read_csv("pre-processing/meta/meta_taouk_2022_SS14.csv")
    df_SS14_2023 = pd.read_csv("pre-processing/meta/meta_beale_2023_SS14.csv")
    ss14_df = pd.concat([df_SS14_2021, df_SS14_2022, df_SS14_2023, df_SS14_2019])
    ss14_df = ss14_df.drop_duplicates(subset='strain')
    ss14_df = ss14_df.drop_duplicates(subset='accession')
    ss14_df['accession'] = ss14_df['accession'].replace('ERS1884626', 'ERR2756209')
    ss14_df['accession'] = ss14_df['accession'].replace('ERS1884591', 'ERR2756180')
    

    df_Nichols_2019 = pd.read_csv("pre-processing/meta/meta_beale_2019_Nichols.csv")
    df_Nichols_2021 = pd.read_csv("pre-processing/meta/meta_beale_2021_Nichols.csv")
    df_Nichols_2022 = pd.read_csv("pre-processing/meta/meta_taouk_2022_Nichols.csv")
    df_Nichols_2023 = pd.read_csv("pre-processing/meta/meta_beale_2023_Nichols.csv")
    nichols_df = pd.concat([df_Nichols_2021, df_Nichols_2022, df_Nichols_2023,df_Nichols_2019])
    nichols_df = nichols_df.drop_duplicates(subset='strain')    
    nichols_df = nichols_df.drop_duplicates(subset='accession')

    # ss14_df['accession'].to_csv('accessions/SS14_accessions.txt',index=False, header=False)
    ss14_df.loc[ss14_df['reads_or_assemblies']=='reads', ['accession']].to_csv('accessions/SS14_accessions.txt',index=False, header=False)
    nichols_df.loc[nichols_df['reads_or_assemblies']=='reads', ['accession']].to_csv('accessions/Nichols_accessions.txt',index=False, header=False)
    ss14_df.loc[ss14_df['reads_or_assemblies']=='reads', ['strain']].to_csv('clusters/SS14.txt',index=False, header=False)
    nichols_df.loc[nichols_df['reads_or_assemblies']=='reads', ['strain']].to_csv('clusters/Nichols.txt',index=False, header=False)

    ss14_df.to_csv('pre-processing/meta/SS14_beale_2019_2021_2023_taouk_2022_meta.csv',index=False)
    nichols_df.to_csv('pre-processing/meta/Nichols_beale_2019_2021_2023_taouk_2022_meta.csv',index=False)


if __name__ == "__main__":
    main()