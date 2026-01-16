import streamlit as st
import pandas as pd

from util import (
    get_multi_app_input,
    add_dumb_mobile_buffer,
    validate_input,
    construct_isotope_df,
    construct_multi_isotope_figure,
    get_query_params_url,
    shorten_url,
)

st.set_page_config(page_title="Multi-IsoCalc", page_icon="📊")
st.logo(
    image="yates_lab.png",
    size='large',
    link="https://github.com/pgarrett-scripps/peptacular/blob/opt/peptacular_logo.png?raw=true",
    icon_image="spectra.png",
)
with st.sidebar:
    st.title("Multi-IsoCalc", text_alignment="center")
    st.markdown(
        """
        Calculate and visualize the isotopic distributions of multiple peptides, formulas, or neutral masses. 
        Peptide and Formula inputs must be 
        [proforma2.0 compliant](https://peptacular.readthedocs.io/en/latest/modules/getting_started.html#proforma-notation).
        Neutral Mass will be converted to a composition using the averagine peptide.
        """,
        text_alignment="center",
    )

    st.subheader("Options", divider="grey")
    params = get_multi_app_input()
    add_dumb_mobile_buffer()


title_c, _, button_c = st.columns([2, 1, 1])
help_msg = "This page's URL automatically updates with your input and can be shared with others. You can optionally use the Generate TinyURL button to create a shortened URL."
title_c.header("Results", help=help_msg)

df = pd.DataFrame()
try:
    for single_iso_param in params.single_iso_inputs:
        validate_input(single_iso_param)
        single_df = construct_isotope_df(single_iso_param)
        single_df["sequence"] = single_iso_param.sequence
        df = pd.concat([df, single_df], ignore_index=True)

    
    df["relative_abundance"] = None
    max_abundance = df["abundance"].max()
    sum_abundance = df["abundance"].sum()
    for single_iso_param in params.single_iso_inputs:
        # Normalize the abundances to sum to 100%
        seq_flag = df["sequence"] == single_iso_param.sequence
        if params.is_intensity_sum:
            # Normalize to sum to 100%
            df.loc[seq_flag, "relative_abundance"] = (
                df.loc[seq_flag, "abundance"] / sum_abundance
            )
        else:
            # Normalize to the most abundant peak
            df.loc[seq_flag, "relative_abundance"] = (
                df.loc[seq_flag, "abundance"] / max_abundance
            )

    # Create a summary dataframe for all parameters
    summary_data: list[dict[str, object]] = []
    for single_iso_param in params.single_iso_inputs:
        summary_data.append({
            "Sequence": single_iso_param.sequence,
            "Neutral Mass (Da)": round(single_iso_param.neutral_mass, 5),
            "m/z": round(single_iso_param.mz, 5),
            "Composition": single_iso_param.chemical_formula
        })

    summary_df = pd.DataFrame(summary_data)
    st.dataframe(summary_df, hide_index=True, width='stretch') # type: ignore

    fig = construct_multi_isotope_figure(
        df, line_width=params.line_width, is_log=params.is_log
    )
    st.plotly_chart(fig) # type: ignore

    # Show isotope table
    # st.title("Isotopic Distribution Table")

    height = min(int(35.3 * (len(df) + 1)), 1000)
    st.dataframe( # type: ignore
        df,
        width='stretch',
        hide_index=True,
        height=height,
        column_order=[
            "sequence",
            "neutral_mass",
            "mz",
            "abundance",
            "relative_abundance",
        ],
        column_config={
            "sequence": st.column_config.TextColumn(
                "Sequence",
                help="Sequence of the peptide.",
                width="small",
            ),
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
                format="%.2f",
            ),
            "relative_abundance": st.column_config.NumberColumn(
                "Relative Abundance",
                help="Abundance of the isotope.",
                width="small",
                # make percentage
                format="percent",
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
        url_params: dict[str, list[str]] = {k: st.query_params.get_all(k) for k in st.query_params.keys()}
        page_url = f"{st.context.url}{get_query_params_url(url_params)}"
        short_url = shorten_url(page_url)

        @st.dialog(title="Share your results")
        def url_dialog(url: str):
            st.write(f"Shortened URL: {url}")

        url_dialog(short_url)

    st.divider()

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("[**Multi-Iso-Calc**](https://github.com/pgarrett-scripps/PeptideIsotopeCalculator)")
        st.image("https://zenodo.org/badge/779470286.svg", width=150)
    with col2:
        st.markdown("[**Peptacular**](https://github.com/pgarrett-scripps/peptacular)")
        st.image("https://zenodo.org/badge/591504879.svg", width=150)


except Exception as e:
    st.error(f"Error validating input: {e}")

