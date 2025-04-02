import streamlit as st
import streamlit_permalink as stp
import peptacular as pt
import pandas as pd
import plotly.graph_objects as go
import pandas as pd

st.set_page_config(
    page_title="Peptide Isotopic Distribution Calculator", page_icon="💻"
)

with st.sidebar:
    st.title("Peptide Isotopic Distribution Calculator 💻")

    st.caption(
        "Calculate the isotopic distribution of one or more peptides. Peptide and Formula inputs "
        "must be proforma2.0 compliant."
    )

    st.caption(
        "Made using [peptacular](https://github.com/pgarrett-scripps/peptacular): [![DOI](https://zenodo.org/badge/591504879.svg)](https://doi.org/10.5281/zenodo.15054278)"
    )

    st.subheader("Peptide Input")
    default_df = pd.DataFrame(
        {
            "sequence": ["PEPTIDE", "PEPTIDEPEPTIDE"],
            "charge": [2, 4],
            "ion_type": ["p", "p"],
            "intensity": [1000, 2000],
            "input_type": ["peptide", "peptide"],
        }
    )

    df = stp.data_editor(
        default_df,
        column_config={
            "sequence": st.column_config.TextColumn(
                "Sequence",
                default="PEPTIDE",
                help="Peptide sequence. Must be proforma2.0 compliant.",
                width="small",
                required=True,
            ),
            "charge": st.column_config.NumberColumn(
                "Charge", default=2, help="Charge state of the peptide.", width="small"
            ),
            "ion_type": st.column_config.SelectboxColumn(
                "Ion Type",
                help="Ion type to calculate the mass of the peptide.",
                default="p",
                options=list("abcxyzpi"),
                width="small",
                required=True,
            ),
            "intensity": st.column_config.NumberColumn(
                "Intensity",
                default=1000,
                help="Intensity of the peptide.",
                width="small",
                required=True,
            ),
            "input_type": st.column_config.SelectboxColumn(
                "Input Type",
                default="peptide",
                help="Input type of the peptide.",
                options=["peptide", "formula", "mass"],
                width="small",
                required=True,
            ),
        },
        hide_index=True,
        use_container_width=True,
        key="peptide_input",
        num_rows="dynamic",
        compress=True,
    )

    df = df.dropna(axis=0, how="any")
    df["charge"] = df["charge"].astype(int)
    df["intensity"] = df["intensity"].astype(float)

    # reset index
    df.reset_index(drop=True, inplace=True)

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
    with c1:
        use_neutron = stp.toggle(
            label="Neutron Offset",
            value=True,
            help="Use neutron offsets for isotopic distribution instead of mass of each isotope. "
            "This gives isotopic peaks at Neutral Mass +/- N * Neutron Mass",
        )

    with c2:
        use_intensity_sum = stp.toggle(
            label="Use Intensity Sum",
            value=False,
            help="Use intensity sum instead of absolute abundance for isotopic distribution.",
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


all_isotopes = []
masses = []
mzs = []
comps = []
for i, row in df.iterrows():
    sequence = row["sequence"]
    ion_type = row["ion_type"]
    charge = row["charge"]
    input_type = row["input_type"]
    intensity = row["intensity"]

    if input_type == "peptide":

        try:
            sequence_mass = pt.mass(sequence, ion_type=ion_type)
        except ValueError as e:
            st.error(e)
            st.stop()

        try:
            sequence_mz = pt.mz(sequence, charge=charge, ion_type=ion_type)
        except ValueError as e:
            st.error(e)
            st.stop()

        try:
            composition, delta_mass = pt.comp_mass(sequence, ion_type)
        except ValueError as e:
            st.error(e)
            st.stop()

        if delta_mass:
            st.warning(
                "Ambiguous Modification! Using averagine Formula to estimate delta mass composition."
            )

            delta_mass_comp = pt.estimate_element_counts(
                delta_mass, pt.get_isotope_mods(sequence)
            )
            combined_comp = {
                k: composition.get(k, 0) + delta_mass_comp.get(k, 0)
                for k in set(composition) | set(delta_mass_comp)
            }
            composition = {k: round(v, 2) for k, v in combined_comp.items()}

        formula = pt.write_chem_formula(composition, hill_order=True)

        isotopes = pt.isotopic_distribution(
            composition,
            max_isotopes,
            min_abundance_threshold,
            distribution_resolution,
            use_neutron,
        )

    elif input_type == "formula":

        try:
            sequence_mass = pt.chem_mass(sequence)
        except ValueError as e:
            st.error(e)
            st.stop()

        try:
            sequence_mz = pt.chem_mz(sequence, charge=charge)
        except ValueError as e:
            st.error(e)
            st.stop()

        composition = pt.parse_chem_formula(sequence)

        if len(composition) >= 10:
            st.warning(
                "The formula is too complex (Limit of 10 unique elements). Please check the input and try again."
            )
            st.stop()

        formula = pt.write_chem_formula(composition)

        isotopes = pt.isotopic_distribution(
            composition,
            max_isotopes,
            min_abundance_threshold,
            distribution_resolution,
            use_neutron,
        )

    elif input_type == "mass":

        sequence_mass = float(sequence)
        try:
            composition = pt.estimate_comp(sequence_mass)
        except ValueError as e:
            st.error(e)
            st.stop()

        sequence_mz = (sequence_mass + charge * pt.constants.PROTON_MASS) / charge

        composition = {k: round(v, 2) for k, v in composition.items()}

        isotopes = pt.estimate_isotopic_distribution(
            sequence_mass,
            max_isotopes,
            min_abundance_threshold,
            distribution_resolution,
            use_neutron,
        )


    print(isotopes)

    if use_neutron:
        # add mass to the isotopes
        isotopes = [
            (sequence_mz + (mass * neutron_value / charge), abundance)
            for mass, abundance in isotopes
        ]
    else:
        # add mass to the isotopes
        isotopes = [
            ((mass + pt.PROTON_MASS * charge) / charge, abundance)
            for mass, abundance in isotopes
        ]

    if use_intensity_sum:
        # makes the isotopes abundancees sum to abundance
        total_abundance = sum([abundance for _, abundance in isotopes])
        isotopes = [
            (mass, (abundance / total_abundance) * intensity) for mass, abundance in isotopes
        ]
    else:
        isotopes = [
            (mass, abundance * intensity) for mass, abundance in isotopes
        ]

    all_isotopes.append(isotopes)
    masses.append(sequence_mass)
    mzs.append(sequence_mz)
    comps.append(composition)

df["neutral_mass"] = masses
df["mass_to_charge_ratio"] = mzs
df["composition"] = comps
df['isotope'] = all_isotopes


st.title("Results")


# st.plotly_chart(fig)

# merge all isotopes, if they have the same mass (rounded to distribution_resolution)
isotope_dict = {}
for isotopes in all_isotopes:
    for mz, abundance in isotopes:
        mz = round(mz, distribution_resolution)
        isotope_dict[mz] = isotope_dict.get(mz, 0) + abundance

# normalize the abundance so that the highest peak is 100%
total_abundance = sum(isotope_dict.values())
max_abundance = max(isotope_dict.values())
relative_isotope_dict = {k: v / max_abundance for k, v in isotope_dict.items()}
absolute_isotope_dict = {k: v for k, v in isotope_dict.items()}

df = pd.DataFrame(
    relative_isotope_dict.items(),
    columns=["mass_to_charge_ratio", "relative_abundance"],
)

# add absolute abundance
df["absolute_abundance"] = df["mass_to_charge_ratio"].map(absolute_isotope_dict)

# filter the df
df = df[df["relative_abundance"] > min_abundance_threshold]
df["relative_abundance"] = df["relative_abundance"] * 100
df["percentage_abundance"] = (
    df["relative_abundance"] / df["relative_abundance"].sum() * 100
)

# make percent
df["relative_abundance"] = df["relative_abundance"].round(3)


# sort the df
df = df.sort_values(by="relative_abundance", ascending=False)

# filter df by max isotopes
df = df.head(max_isotopes)


show_absolute_abundance = stp.toggle(
    label="Show Absolute Abundance",
    value=True,
    help="Show absolute abundance instead of relative abundance.",
)

if show_absolute_abundance:
    column_to_show = "absolute_abundance"
else:
    column_to_show = "relative_abundance"

# Assuming `df` is your DataFrame with isotopes
fig = go.Figure()

# Add scatter plot for the isotopes markers
fig.add_trace(
    go.Scatter(
        x=df["mass_to_charge_ratio"],
        y=df[column_to_show],
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
            y=[0, row[column_to_show]],
            mode="lines",
            line=dict(color="grey", width=1.5),
            showlegend=False,
        )
    )

# Customize the plot for better readability
fig.update_layout(
    title=f"Combined Isotopic Distribution",
    xaxis_title="Mass To Charge Ratio",
    yaxis_title="Abundance",
    margin=dict(l=40, r=40, t=40, b=40),
)

st.plotly_chart(fig)

# center title
# st.markdown("<center><h1>Isotopic Distribution Table</h1></center>", unsafe_allow_html=True)
st.title("Isotopic Distribution Table")
# st.dataframe(df, use_container_width=True, hide_index=True)

# reordering the columns
df = df[["absolute_abundance", "mass_to_charge_ratio", "relative_abundance"]]
df["percentage_abundance"] = (
    df["relative_abundance"] / df["relative_abundance"].sum() * 100
)

# format the columns
df["mass_to_charge_ratio"] = df["mass_to_charge_ratio"].round(4)
df["relative_abundance"] = df["relative_abundance"].round(3)
df["percentage_abundance"] = df["percentage_abundance"].round(3)


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
        "absolute_abundance": st.column_config.NumberColumn(
            "Absolute Abundance",
            help="Absolute abundance of the isotope.",
            width="small",
        ),
    },
)
