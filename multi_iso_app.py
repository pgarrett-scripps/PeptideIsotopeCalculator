import streamlit as st
import peptacular as pt
import pandas as pd
import plotly.graph_objects as go
import pandas as pd
from streamlit_js_eval import get_page_location

from util import (get_multi_app_input, 
                  add_dumb_mobile_buffer, 
                  validate_input, 
                  construct_isotope_df,
                  construct_multi_isotope_figure,
                    get_query_params_url,
                    shorten_url)

st.set_page_config(
    page_title="Multi-IsoCalc", page_icon="📊"
)

with st.sidebar:
    st.title("Multi-Isotopic Distribution Calculator 📊")

    st.caption(
        "Calculate and visualize the isotopic distributions of multiple peptides, formulas, or neutral masses. Peptide and Formula inputs "
        "must be proforma2.0 compliant. Neutral Mass is in Daltons (Da). Distributions are calculated prior "
        "to applying charge state or ambiguous masses (so charge adducts and delta mass mods are ignored)."
    )

    st.caption(
        "Made using [peptacular](https://github.com/pgarrett-scripps/peptacular): [![DOI](https://zenodo.org/badge/591504879.svg)](https://doi.org/10.5281/zenodo.15054278)"
    )

    st.subheader("Input Table")

    params = get_multi_app_input()
    add_dumb_mobile_buffer()


top_window, bottom_window = st.container(), st.container()

with bottom_window:
    page_loc = get_page_location()

with top_window:


    title_c, _, button_c = st.columns([2, 1, 1])
    help_msg = "This page's URL automatically updates with your input and can be shared with others. You can optionally use the Generate TinyURL button to create a shortened URL."
    title_c.header("Results", help=help_msg)

    df = pd.DataFrame()
    for single_iso_param in params.single_iso_inputs:
        validate_input(single_iso_param)
        single_df = construct_isotope_df(single_iso_param)
        single_df['sequence'] = single_iso_param.sequence
        df = pd.concat([df, single_df], ignore_index=True)

    # Find maximum abundance across all m/z for relative scaling
    grouped_df = df.groupby('mass_to_charge_ratio')
    total_by_mz = grouped_df['abundance'].sum()
    max_total_abundance = total_by_mz.max()

    # Create a single table with multiple rows for all parameters
    table_html = """
    <table style="width:100%; margin-bottom:20px; border-collapse:collapse; border-radius:4px; overflow:hidden; box-shadow:0 2px 3px rgba(0,0,0,0.1);">
        <thead>
            <tr style="background-color:#f2f2f2;">
                <th style="padding:12px 15px; text-align:left; border-bottom:2px solid #dddddd;">Sequence</th>
                <th style="padding:12px 15px; text-align:left; border-bottom:2px solid #dddddd;">Neutral Mass (Da)</th>
                <th style="padding:12px 15px; text-align:left; border-bottom:2px solid #dddddd;">m/z</th>
                <th style="padding:12px 15px; text-align:left; border-bottom:2px solid #dddddd;">Composition</th>
            </tr>
        </thead>
        <tbody>
    """

    # Add a row for each parameter
    for single_iso_param in params.single_iso_inputs:
        table_html += f"""
            <tr>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{single_iso_param.sequence}</td>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{round(single_iso_param.neutral_mass, 5)}</td>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{round(single_iso_param.mz, 5)}</td>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{single_iso_param.chemical_formula}</td>
            </tr>
        """

    # Close the table
    table_html += """
        </tbody>
    </table>
    """

    # Display the complete table once
    st.html(table_html)

    fig = construct_multi_isotope_figure(df)
    st.plotly_chart(fig)

    # Show isotope table
    st.title("Isotopic Distribution Table")

    # fix realtive intensity (should be realtive to the sum of all rows with the same mass)
    # reordering the columns
    # Group by mass_to_charge_ratio to identify overlapping isotopes



    height = min(int(35.2 * (len(df) + 1)), 1000)
    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        height=height,
        column_order=[
            "sequence",
            "neutral_mass",
            "mass_to_charge_ratio",
            "abundance",
            "global_realtive_abundance",
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
            "mass_to_charge_ratio": st.column_config.NumberColumn(
                "Mass to Charge Ratio",
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
                "global_realtive_abundance": st.column_config.NumberColumn(
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
        type='primary',
        use_container_width=True,
        on_click='ignore',
        help="Download the isotopic distribution table as a CSV file.",

    )


    if page_loc and 'origin' in page_loc:
        url_origin = page_loc['origin']
        if button_c.button("Generate TinyURL", key="generate_tinyurl", type="primary"):
            url_params = {k: st.query_params.get_all(k) for k in st.query_params.keys()}
            page_url = f"{url_origin}{get_query_params_url(url_params)}"
            short_url = shorten_url(page_url)


            @st.dialog(title="Share your results")
            def url_dialog(url):
                st.write(f"Shortened URL: {url}")

            url_dialog(short_url)