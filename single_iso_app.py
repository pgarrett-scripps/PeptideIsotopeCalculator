import streamlit as st
from streamlit_js_eval import get_page_location


from util import (get_single_app_input, 
                  add_dumb_mobile_buffer, 
                  validate_input, 
                  construct_isotope_df,
                  construct_figure,
                  get_query_params_url,
                  shorten_url)

st.set_page_config(
    layout="centered", page_title="IsoCalc", page_icon="📊"
)

with st.sidebar:

    # Inject custom CSS to set the width of the sidebar
    st.markdown(
        """
        <style>
            section[data-testid="stSidebar"] {
                width: 600px !important; # Set the width to your desired value
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(f"""
        <div style='text-align: center; padding: 15px; top-margin: 0px'>
            <h3 style='margin: 0; font-size: 1.5em; color: #333;'>IsoCalc 📊</h3>
            <p style='font-size: 1.0em; line-height: 1.6; color: #555;'>
                Calculate and visualize the isotopic distribution for a peptide, formula, or neutral mass. 
                Peptide and Formula inputs must be 
                <a href="https://peptacular.readthedocs.io/en/latest/modules/getting_started.html#proforma-notation" 
                target="_blank" style='color: #007BFF; text-decoration: none;'>proforma2.0 compliant</a>.
                Neutral Mass is in Daltons (Da). Distributions are calculated prior to applying charge state 
                or ambiguous masses (so charge adducts and delta mass mods are ignored). Powered by 
                <a href="https://github.com/pgarrett-scripps/peptacular" target="_blank" style='color: #007BFF; text-decoration: none;'>
                    <strong>Peptacular</strong>
                </a>.
            </p>
        </div>
    """, unsafe_allow_html=True)


    # Get all input parameters from the user
    params = get_single_app_input()
    add_dumb_mobile_buffer()


top_window, bottom_window = st.container(), st.container()

with bottom_window:
    page_loc = get_page_location()

with top_window:

    title_c, _, button_c = st.columns([2, 1, 1])
    help_msg = "This page's URL automatically updates with your input and can be shared with others. You can optionally use the Generate TinyURL button to create a shortened URL."
    title_c.header("Results", help=help_msg)
    validate_input(params)
    df = construct_isotope_df(params)
    fig = construct_figure(df, params)

    st.html(f"""
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
            <tr>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{params.sequence}</td>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{round(params.neutral_mass, 5)}</td>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{round(params.mz, 5)}</td>
                <td style="padding:12px 15px; border-bottom:1px solid #dddddd;">{params.chemical_formula}</td>
            </tr>
        </tbody>
    </table>
    """)

    st.plotly_chart(fig)

    # Show isotope table
    st.title("Isotopic Distribution Table")

    height = min(int(35.3 * (len(df) + 1)), 1000)
    st.dataframe(
        df,
        use_container_width=True,
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
        type='secondary',
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


    st.divider()

    st.markdown(f"""
    <div style='display: flex; justify-content: space-between; align-items: center; padding: 15px 0; border-top: 0px solid #ddd;'>
        <div style='text-align: left; font-size: 1.1em; color: #555;'>
            <a href="https://github.com/pgarrett-scripps/PeptideIsotopeCalculator" target="_blank" 
                style='text-decoration: none; color: #007BFF; font-weight: bold;'>
                Iso-Calc
            </a>
            <a href="https://doi.org/10.5281/zenodo.15170926" target="_blank" style="margin-left: 12px;">
                <img src="https://zenodo.org/badge/779470286.svg" alt="DOI" 
                        style="vertical-align: middle; height: 20px;">
            </a>
        </div>
        <div style='text-align: right; font-size: 1.1em; color: #555;'>
            <a href="https://github.com/pgarrett-scripps/peptacular" target="_blank" 
                style='text-decoration: none; color: #007BFF; font-weight: bold;'>
                Peptacular
            </a>
            <a href="https://doi.org/10.5281/zenodo.15054278" target="_blank" style="margin-left: 12px;">
                <img src="https://zenodo.org/badge/591504879.svg" alt="DOI" 
                        style="vertical-align: middle; height: 20px;">
            </a>
        </div>
    </div>
""", unsafe_allow_html=True)
