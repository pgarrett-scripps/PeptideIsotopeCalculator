import streamlit as st
import streamlit_permalink as stp
import peptacular as pt
import pandas as pd
import plotly.graph_objects as go


st.set_page_config(
    layout="centered", page_title="Isotopic Distribution Calculator", page_icon="📊"
)

with st.sidebar:
    st.title("Isotopic Distribution Calculator 📊")

    st.caption(
        "Calculate the isotopic distribution of a peptide, formula, or neutral mass. Peptide and Formula inputs "
        "must be proforma2.0 compliant. Neutral Mass is in Daltons (Da). Distributions are calculated prior "
        "to applying charge state (so proton isotopic abundance is not considered)."
    )

    st.caption(
        "Made using [peptacular](https://github.com/pgarrett-scripps/peptacular): [![DOI](https://zenodo.org/badge/591504879.svg)](https://doi.org/10.5281/zenodo.15054278)"
    )

    input_type = stp.radio(
        label="Input Type",
        options=["Peptide", "Formula", "Mass"],
        index=0,
        horizontal=True,
    )

    if input_type == "Peptide":
        sequence_input = stp.text_input(
            label="Sequence",
            value="PEPTIDE",
            help="Enter the peptide sequence (e.g. PEPTIDE)",
        )

        sequence_input = sequence_input.replace(" ", "")
        sequence_input = sequence_input.replace("\n", "")
        sequence_input = sequence_input.strip()

        if not sequence_input:
            st.warning("Please enter a sequence.")
            st.stop()

        empty_caption = st.empty()
        ion_type = stp.radio(
            label="Ion Type",
            options=list("abcxyzpi"),
            index=6,
            horizontal=True,
            help="Ion type to calculate, p = precursor",
        )

        try:
            sequence_mass = pt.mass(sequence_input, ion_type=ion_type)
        except ValueError as e:
            st.error(e)
            st.stop()

        try:
            composition, delta_mass = pt.comp_mass(sequence_input, ion_type)
        except ValueError as e:
            st.error(e)
            st.stop()

        if delta_mass:
            st.warning(
                "Ambiguous Modification! Using averagine Formula to estimate delta mass composition."
            )

            delta_mass_comp = pt.estimate_element_counts(
                delta_mass, pt.get_isotope_mods(sequence_input)
            )
            combined_comp = {
                k: composition.get(k, 0) + delta_mass_comp.get(k, 0)
                for k in set(composition) | set(delta_mass_comp)
            }
            composition = {k: round(v, 2) for k, v in combined_comp.items()}

        formula = pt.write_chem_formula(composition, hill_order=True)
        empty_caption.caption(
            f"Neutral Mass: {sequence_mass:.4f} (Da), Chem Formula: {formula}"
        )

    elif input_type == "Formula":
        sequence_input = stp.text_input(
            label="Formula",
            value="C6H12O6",
            help="Enter the chemical formula (e.g. C6H12O6). Isotopes must be enclosed in "
            "brackets (e.g. C[13]H12O6)",
        )

        empty_caption = st.empty()
        try:
            sequence_mass = pt.chem_mass(sequence_input)
        except ValueError as e:
            st.error(e)
            st.stop()

        st.caption(f"Neutral Mass: {sequence_mass:.4f} Da")
        composition = pt.parse_chem_formula(sequence_input)

        if len(composition) >= 10:
            st.warning(
                "The formula is too complex (Limit of 10 unique elements). Please check the input and try again."
            )
            st.stop()

        formula = pt.write_chem_formula(composition)

    elif input_type == "Mass":
        sequence_input = stp.number_input(
            label="Neutral (Monoisotopic) Mass",
            value=1000.0,
            max_value=100_001.0,
            min_value=0.0,
            step=100.0,
            help="Enter the neutral mass (e.g. 1000 Da)",
        )
        sequence_mass = sequence_input
        try:
            composition = pt.estimate_comp(sequence_mass)
        except ValueError as e:
            st.error(e)
            st.stop()

        composition = {k: round(v, 2) for k, v in composition.items()}
        st.caption(
            f"Chem Formula: {pt.write_chem_formula(composition, hill_order=True)}"
        )

    charge = stp.number_input(
        label="Charge",
        value=1,
        min_value=0,
        max_value=100,
        step=1,
        help="Enter the charge state (e.g. 2+)",
    )

    c1, c2 = st.columns(2)
    with c1:
        max_isotopes = stp.number_input(
            label="Max Number of Isotopes",
            value=100,
            min_value=1,
            max_value=200,
            step=1,
            help="Maximum number of isotopes to display",
        )

    with c2:
        min_abundance_threshold = stp.selectbox(
            label="Min Relative Abundance",
            options=[0.00001, 0.0001, 0.001, 0.01, 0.1],
            index=2,
            help="Minimum abundance threshold to display",
        )

    use_neutron = stp.toggle(
        label="Neutron Offset",
        value=True,
        help="Use neutron offsets for isotopic distribution instead of mass of each isotope. "
        "This gives isotopic peaks at Neutral Mass +/- N * Neutron Mass",
    )

    distribution_resolution = 2
    if not use_neutron:
        distribution_resolution = stp.slider(
            label="Resolution",
            value=2,
            min_value=0,
            max_value=10,
            step=1,
            help="Resolution of the distribution (Round to nearest X decimal places)",
        )
    else:
        neutron_value = stp.selectbox(
            label="Neutron Mass",
            options=["Neutron", "C13-C12", "Averagine (Peptide)", "Custom"],
            index=2,
            help="Neutron mass to use for isotopic distribution calculation",
        )

        if neutron_value == "Neutron":
            neutron_value = pt.NEUTRON_MASS
        elif neutron_value == "C13-C12":
            neutron_value = 1.00335
        elif neutron_value == "Averagine (Peptide)":
            neutron_value = 1.002856
        elif neutron_value == "Custom":
            neutron_value = stp.number_input(
                label="Custom Neutron Mass",
                value=pt.NEUTRON_MASS,
                format="%.6f",
                help="Custom neutron mass to use for isotopic distribution calculation",
            )

        st.caption(f"Neutron Mass: {neutron_value:.6f} Da")

    if sequence_mass >= 100_000:
        st.error("The mass is too high. Please check the input and try again.")
        st.stop()

    # blank space for mobile cause streamlit is weird
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)

if input_type == "Peptide":

    isotopes = pt.isotopic_distribution(
        composition,
        max_isotopes,
        min_abundance_threshold,
        distribution_resolution,
        use_neutron,
    )

elif input_type == "Formula":

    isotopes = pt.isotopic_distribution(
        composition,
        max_isotopes,
        min_abundance_threshold,
        distribution_resolution,
        use_neutron,
    )

elif input_type == "Mass":

    isotopes = pt.estimate_isotopic_distribution(
        sequence_input,
        max_isotopes,
        min_abundance_threshold,
        distribution_resolution,
        use_neutron,
    )
if len(isotopes) == 0:
    st.error("No isotopes found.")
    st.stop()

if use_neutron:
    # add mass to the isotopes
    isotopes = [
        (sequence_mass + mass * neutron_value, abundance)
        for mass, abundance in isotopes
    ]

df = pd.DataFrame(isotopes, columns=["neutral_mass", "relative_abundance"])

# filter the df
df = df[df["relative_abundance"] > min_abundance_threshold]
df["relative_abundance"] = df["relative_abundance"] * 100

# make percent
df["relative_abundance"] = df["relative_abundance"].round(3)

# sort the df
df = df.sort_values(by="relative_abundance", ascending=False)

# filter df by max isotopes
df = df.head(max_isotopes)

# apply the charge state
df["mass_to_charge_ratio"] = (
    (df["neutral_mass"] + charge * pt.PROTON_MASS) / charge
    if charge > 0
    else df["neutral_mass"]
)

st.title("Results")

# Assuming `df` is your DataFrame with isotopes
fig = go.Figure()

# Add scatter plot for the isotopes markers
fig.add_trace(
    go.Scatter(
        x=df["mass_to_charge_ratio"],
        y=df["relative_abundance"],
        mode="markers",
        marker=dict(size=4),
        name="Isotopes",
    )
)

# Add lines from each isotope marker to the base
for idx, row in df.iterrows():
    fig.add_trace(
        go.Scatter(
            x=[row["mass_to_charge_ratio"], row["mass_to_charge_ratio"]],
            y=[0, row["relative_abundance"]],
            mode="lines",
            line=dict(color="grey", width=1.5),
            showlegend=False,
        )
    )

# Customize the plot for better readability
fig.update_layout(
    title=f'Isotopic Distribution: {sequence_input} {"da" if isinstance(sequence_input, int) else ""}',
    xaxis_title="Mass To Charge Ratio",
    yaxis_title="Relative Abundance",
    margin=dict(l=40, r=40, t=40, b=40),
)

st.plotly_chart(fig)

# center title
# st.markdown("<center><h1>Isotopic Distribution Table</h1></center>", unsafe_allow_html=True)
st.title("Isotopic Distribution Table")
# st.dataframe(df, use_container_width=True, hide_index=True)

# reordering the columns
df = df[["neutral_mass", "mass_to_charge_ratio", "relative_abundance"]]
df["percentage_abundance"] = (
    df["relative_abundance"] / df["relative_abundance"].sum() * 100
)

# format the columns
df["neutral_mass"] = df["neutral_mass"].round(4)
df["mass_to_charge_ratio"] = df["mass_to_charge_ratio"].round(4)
df["relative_abundance"] = df["relative_abundance"].round(3)
df["percentage_abundance"] = df["percentage_abundance"].round(3)


# download
height = min(int(35.2 * (len(df) + 1)), 1000)
st.dataframe(
    df,
    use_container_width=True,
    hide_index=True,
    height=height,
    column_config={
        "neutral_mass": st.column_config.NumberColumn(
            "Neutral Mass", help="Neutral mass of the isotope.", width="small"
        ),
        "mass_to_charge_ratio": st.column_config.NumberColumn(
            "Mass to Charge Ratio",
            help="Mass to charge ratio of the isotope.",
            width="small",
        ),
        "relative_abundance": st.column_config.NumberColumn(
            "Relative Abundance",
            help="Relative abundance of the isotope.",
            width="small",
        ),
        "percentage_abundance": st.column_config.NumberColumn(
            "Percentage Abundance",
            help="Percentage abundance of the isotope.",
            width="small",
        ),
    },
)
