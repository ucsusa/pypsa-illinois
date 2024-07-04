import pandas as pd
from unyt import g, kg, tonne, pound, kWh, MWh, GWh


if __name__ == "__main__":

    dfs = pd.read_html("https://www.eia.gov/tools/faqs/faq.php?id=74&t=11")
    df = dfs[0].head(3)
    df = df.replace({'Natural gas':'Natural Gas'})
    df.columns = df.columns.droplevel()
    df = df.set_index(df.iloc[:,0])
    df = df[['pounds per kWh']]

    df = df.iloc[:,0].astype('float').apply(lambda x: (x*(pound/kWh)).to(tonne/MWh)).to_frame()

    df.index.name = ''
    df.columns = ['tCO2 per MWh']
    df.to_csv(snakemake.output.emissions)