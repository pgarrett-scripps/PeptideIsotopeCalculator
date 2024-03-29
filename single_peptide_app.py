import base64

import streamlit as st
import peptacular as pt
import pandas as pd
import plotly.graph_objects as go

with st.sidebar:
    st.title('Isotope Distribution Calculator')

    st.caption('Calculate the isotopic distribution of a peptide, formula, or neutral mass. Peptide and Formula inputs '
               'must be proforma2.0 compliant. Neutral Mass is in Daltons (Da). Distributions are calculated prior '
               'to applying charge state (so proton isotopic abundance is not considered).')

    st.caption('Made with [peptacular](https://pypi.org/project/peptacular/)')

    input_type = st.radio(label='Input Type',
                          options=['Peptide', 'Formula', 'Neutral Mass'],
                          index=0,
                          horizontal=True)

    c1, c2 = st.columns([7, 3])

    charge = c2.number_input(label='Charge',
                             value=1,
                             min_value=1,
                             max_value=100,
                             step=1,
                             help='Enter the charge state (e.g. 2+)')

    if input_type == 'Peptide':
        sequence_input = c1.text_input(label='Sequence',
                                       value='PEPTIDE',
                                       help='Enter the peptide sequence (e.g. PEPTIDE)')

        sequence_input = sequence_input.replace(' ', '')
        sequence_input = sequence_input.replace('\n', '')
        sequence_input = sequence_input.strip()

        if not sequence_input:
            st.warning('Please enter a sequence.')
            st.stop()

        empty_caption = st.empty()
        ion_type = st.radio(label='Ion Type',
                            options=['p', 'b', 'y', 'a', 'c', 'x', 'z', 'i'],
                            index=0,
                            horizontal=True,
                            help='Ion type to calculate, p = precursor')

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
            st.warning('Ambiguous Modification! Using averagine Formula to estimate delta mass composition.')

            delta_mass_comp = pt.estimate_element_counts(delta_mass, pt.get_isotope_mods(sequence_input))
            combined_comp = {k: composition.get(k, 0) + delta_mass_comp.get(k, 0)
                             for k in set(composition) | set(delta_mass_comp)}
            composition = {k: round(v, 2) for k, v in combined_comp.items()}

        formula = pt.write_chem_formula(composition, hill_order=True)
        empty_caption.caption(f'Neutral Mass: {sequence_mass:.4f} (Da), Chem Formula: {formula}')

    elif input_type == 'Formula':
        sequence_input = c1.text_input(label='Formula',
                                       value='C6H12O6',
                                       help='Enter the chemical formula (e.g. C6H12O6). Isotopes must be enclosed in '
                                            'brackets (e.g. C[13]H12O6)')

        empty_caption = st.empty()
        try:
            sequence_mass = pt.chem_mass(sequence_input)
        except ValueError as e:
            st.error(e)
            st.stop()

        st.caption(f'Neutral Mass: {sequence_mass:.4f} Da')
        composition = pt.parse_chem_formula(sequence_input)

        if len(composition) >= 10:
            st.warning(
                'The formula is too complex (Limit of 10 unique elements). Please check the input and try again.')
            st.stop()

        formula = pt.write_chem_formula(composition)

    elif input_type == 'Neutral Mass':
        sequence_input = c1.number_input(label='Neutral (Monoisotopic) Mass',
                                         value=1000.0,
                                         max_value=100_001.0,
                                         min_value=0.0,
                                         step=100.0,
                                         help='Enter the neutral mass (e.g. 1000 Da)')
        sequence_mass = sequence_input
        try:
            composition = pt.estimate_element_counts(sequence_mass)
        except ValueError as e:
            st.error(e)
            st.stop()

        composition = {k: round(v, 2) for k, v in composition.items()}
        st.caption(f'Chem Formula: {pt.write_chem_formula(composition, hill_order=True)}')

    c1, c2 = st.columns(2)
    max_isotopes = c1.number_input(label='Max Number of Isotopes',
                                   value=100,
                                   min_value=1,
                                   max_value=200,
                                   step=1,
                                   help='Maximum number of isotopes to display')

    min_abundance_threshold = c2.selectbox(label='Min Relative Abundance',
                                           options=[0.00001, 0.0001, 0.001, 0.01, 0.1],
                                           index=2,
                                           help='Minimum abundance threshold to display')

    use_neutron = st.toggle(label='Neutron Offset',
                            value=False,
                            help='Use neutron offsets for isotopic distribution instead of mass of each isotope. '
                                 'This gives isotopic peaks at Neutral Mass +/- N * Neutron Mass')

    distribution_resolution = st.slider(label='Resolution',
                                        value=5,
                                        min_value=0,
                                        max_value=10,
                                        step=1,
                                        help='Resolution of the distribution (Round to nearest X decimal places)',
                                        disabled=use_neutron)

    st.caption('Low Resolution: 0-1, Medium Resolution: 2-5, High Resolution: 6-10')

    if sequence_mass >= 100_001:
        st.error('The mass is too high. Please check the input and try again.')
        st.stop()

    # blank space for mobile
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)

if input_type == 'Peptide':

    isotopes = pt.isotopic_distribution(composition,
                                        max_isotopes,
                                        min_abundance_threshold,
                                        distribution_resolution,
                                        use_neutron)

elif input_type == 'Formula':

    isotopes = pt.isotopic_distribution(composition,
                                        max_isotopes,
                                        min_abundance_threshold,
                                        distribution_resolution,
                                        use_neutron)

elif input_type == 'Neutral Mass':

    isotopes = pt.estimate_isotopic_distribution(sequence_input,
                                                 max_isotopes,
                                                 min_abundance_threshold,
                                                 distribution_resolution,
                                                 use_neutron)
if len(isotopes) == 0:
    st.error('No isotopes found.')
    st.stop()

if use_neutron:
    # add mass to the isotopes
    isotopes = [(sequence_mass + mass * pt.constants.NEUTRON_MASS, abundance) for mass, abundance in isotopes]

df = pd.DataFrame(isotopes, columns=['neutral_mass', 'relative_abundance'])

# filter the df
df = df[df['relative_abundance'] > min_abundance_threshold]
df['relative_abundance'] = df['relative_abundance'] * 100

# make percent
df['relative_abundance'] = df['relative_abundance'].round(3)

# sort the df
df = df.sort_values(by='relative_abundance', ascending=False)

# filter df by max isotopes
df = df.head(max_isotopes)

# apply the charge state
df['mass_to_charge_ratio'] = (df['neutral_mass'] + charge * 1.007276466) / charge if charge > 0 else df['neutral_mass']

# Assuming `df` is your DataFrame with isotopes
fig = go.Figure()

# Add scatter plot for the isotopes markers
fig.add_trace(go.Scatter(x=df['mass_to_charge_ratio'], y=df['relative_abundance'], mode='markers',
                         marker=dict(size=4), name='Isotopes'))

# Add lines from each isotope marker to the base
for idx, row in df.iterrows():
    fig.add_trace(
        go.Scatter(x=[row['mass_to_charge_ratio'], row['mass_to_charge_ratio']], y=[0, row['relative_abundance']],
                   mode='lines', line=dict(color='grey', width=1.5),
                   showlegend=False))

# Customize the plot for better readability
fig.update_layout(title=f'Isotopic Distribution: {sequence_input} {"da" if isinstance(sequence_input, int) else ""}',
                  xaxis_title="Mass To Charge Ratio", yaxis_title="Relative Abundance",
                  margin=dict(l=40, r=40, t=40, b=40))

st.plotly_chart(fig)

# center title
st.markdown("<center><h1>Isotopic Distribution Table</h1></center>", unsafe_allow_html=True)
# st.dataframe(df, use_container_width=True, hide_index=True)

# reordering the columns
df = df[['neutral_mass', 'mass_to_charge_ratio', 'relative_abundance']]
df['percentage_abundance'] = df['relative_abundance'] / df['relative_abundance'].sum() * 100

# format the columns
df['neutral_mass'] = df['neutral_mass'].round(4)
df['mass_to_charge_ratio'] = df['mass_to_charge_ratio'].round(4)
df['relative_abundance'] = df['relative_abundance'].round(3)
df['percentage_abundance'] = df['percentage_abundance'].round(3)

# center
st.markdown("<style>table{color: #000000; text-align: center;}</style>", unsafe_allow_html=True)


def get_table_download_link(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="isotopic_distribution.csv">Download Table</a>'
    return href


# download
st.markdown(get_table_download_link(df), unsafe_allow_html=True)

filter_by = st.toggle(label='Sort by Mass', value=False)

if filter_by:
    df = df.sort_values(by='neutral_mass')
else:
    df = df.sort_values(by='relative_abundance', ascending=False)

st.markdown(df.style.hide(axis="index").to_html(), unsafe_allow_html=True)
