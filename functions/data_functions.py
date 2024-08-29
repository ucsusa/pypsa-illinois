import pandas as pd


def load_heatrates():
    
    heatrate_url = "https://www.eia.gov/electricity/annual/html/epa_08_01.html"
    heatrates = pd.read_html(heatrate_url)[1].set_index("Year")
    
    return heatrates