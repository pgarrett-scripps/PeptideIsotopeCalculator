import streamlit as st


from util import (
    get_single_app_input,
    add_dumb_mobile_buffer,
    validate_input,
    construct_isotope_df,
    construct_figure,
    get_query_params_url,
    shorten_url,
)

st.set_page_config(layout="centered", page_title="IsoCalc", page_icon="📊")

st.logo(
    image="yates_lab.png",
    size='large',
    link="https://github.com/pgarrett-scripps/peptacular/blob/opt/peptacular_logo.png?raw=true",
    icon_image="spectra.png",
)
with st.sidebar:
    st.title("IsoCalc", text_alignment="center")
    st.markdown(
        """
        Calculate and visualize the isotopic distribution for a peptide, formula, or neutral mass. 
        Peptide and Formula inputs must be 
        [proforma2.0 compliant](https://peptacular.readthedocs.io/en/latest/modules/getting_started.html#proforma-notation).
        Neutral Mass will be converted to a composition using the averagine peptide composition.
        """,
        text_alignment="center",
    )

    # Get all input parameters from the user
    st.subheader("Options", divider="grey")
    params = get_single_app_input()
    add_dumb_mobile_buffer()



title_c, _, button_c = st.columns([2, 1, 1])
help_msg = "This page's URL automatically updates with your input and can be shared with others. You can optionally use the Generate TinyURL button to create a shortened URL."
title_c.header("Results", help=help_msg)
validate_input(params)
df = construct_isotope_df(params)
fig = construct_figure(df, params)

import pandas as pd
summary_df = pd.DataFrame({
    "Sequence": [params.sequence],
    "Neutral Mass (Da)": [round(params.neutral_mass, 5)],
    "m/z": [round(params.mz, 5)],
    "Composition": [params.chemical_formula]
})
st.dataframe(summary_df, hide_index=True, use_container_width=True)

st.plotly_chart(fig)

# Show isotope table
# st.title("Isotopic Distribution Table")

height = min(int(35.3 * (len(df) + 1)), 1000)
st.dataframe(
    df,
    width='stretch',
    hide_index=True,
    height=height,
    column_order=[
        "neutral_mass",
        "mz",
        "abundance",
    ],
    column_config={
        "neutral_mass": st.column_config.NumberColumn(
            "Neutral Mass",
            help="Neutral mass of the isotope.",
            width="small",
            format="%.4f",
        ),
        "mz": st.column_config.NumberColumn(
            "m/z",
            help="Mass to charge ratio of the isotope.",
            width="small",
            format="%.4f",
        ),
        "abundance": st.column_config.NumberColumn(
            "Abundance",
            help="Abundance of the isotope.",
            width="small",
            format="%.4f",
        ),
    },
)

# download df
st.download_button(
    label="Download DataFrame as CSV",
    data=df.to_csv(index=False),
    file_name="isotopic_distribution.csv",
    mime="text/csv",
    type="secondary",
    width='stretch',
    on_click="ignore",
    help="Download the isotopic distribution table as a CSV file.",
)

if button_c.button("Generate TinyURL", key="generate_tinyurl", type="primary"):
    url_params = {k: st.query_params.get_all(k) for k in st.query_params.keys()}
    page_url = f"{st.context.url}{get_query_params_url(url_params)}"
    short_url = shorten_url(page_url)

    @st.dialog(title="Share your results")
    def url_dialog(url: str):
        st.write(f"Shortened URL: {url}")

    url_dialog(short_url)

st.divider()

col1, col2 = st.columns(2)
with col1:
    st.markdown("[**Iso-Calc**](https://github.com/pgarrett-scripps/PeptideIsotopeCalculator)")
    st.image("https://zenodo.org/badge/779470286.svg", width=150)
with col2:
    st.markdown("[**Peptacular**](https://github.com/pgarrett-scripps/peptacular)")
    st.image("https://zenodo.org/badge/591504879.svg", width=150)
