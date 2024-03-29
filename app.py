import base64
from typing import Dict

import streamlit as st
import peptacular as pt
import pandas as pd
import plotly.graph_objects as go

if 'peptide_count' not in st.session_state:
    st.session_state['peptide_count'] = 1

with st.sidebar:
    st.title('Isotope Distribution Calculator')

    st.caption('Calculate the isotopic distribution of a peptide, formula, or neutral mass. Peptide and Formula inputs '
               'must be proforma2.0 compliant. Neutral Mass is in Daltons (Da). Distributions are calculated prior '
               'to applying charge state (so proton isotopic abundance is not considered).')

    st.caption('Made with [peptacular](https://pypi.org/project/peptacular/)')

    c1, c2 = st.columns(2)
    if c1.button('Add Peptide', key='add_peptide', use_container_width=True):
        st.session_state['peptide_count'] += 1

    if c2.button('Remove Peptide', key='remove_peptide', use_container_width=True):

        if st.session_state['peptide_count'] > 1:
            st.session_state['peptide_count'] -= 1
        else:
            st.toast('At least one peptide is required')

    for i in range(st.session_state['peptide_count']):
        with st.expander(f'Peptide {i+1}', expanded=True):
            sequence_input = st.text_input(label='Sequence',
                                           value='PEPTIDE',
                                           help='Enter the peptide sequence (e.g. PEPTIDE)',
                                           key=f'sequence_input_{i}')
            c1, c2, c3 = st.columns(3)
            charge = c1.number_input(label='Charge',
                                        value=2,
                                        min_value=0,
                                        step=1,
                                        help='Charge state of the peptide',
                                        key=f'charge_{i}')

            ion_types = c2.selectbox(label='Ion Type',
                                        options=['p', 'b', 'y', 'a', 'c', 'z', 'x', ],
                                        index=0,
                                        help='Ion type to calculate the mass of the peptide',
                                        key=f'ion_types_{i}')

            intensity = c3.number_input(label='Intensity',
                                        value=1000,
                                        min_value=0,
                                        step=1000,
                                        help='Intensity of the peptide',
                                        key=f'intensity_{i}')


    # Get the values
    sequences = [st.session_state[f'sequence_input_{i}'] for i in range(st.session_state['peptide_count'])]
    charges = [st.session_state[f'charge_{i}'] for i in range(st.session_state['peptide_count'])]
    ion_types = [st.session_state[f'ion_types_{i}'] for i in range(st.session_state['peptide_count'])]
    intensities = [st.session_state[f'intensity_{i}'] for i in range(st.session_state['peptide_count'])]

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
                                        disabled=False)

    st.caption('Low Resolution: 0-1, Medium Resolution: 2-5, High Resolution: 6-10')


data = []
for sequence, charge, ion_types, intensity in zip(sequences, charges, ion_types, intensities):
    data.append({'sequence': sequence, 'charge': int(charge), 'ion_type': ion_types, 'intensity': float(intensity)})

df = pd.DataFrame(data)


def get_chem_composition(sequence, ion_type, charge) -> Dict[str, float]:
    try:
        composition, delta_mass = pt.comp_mass(sequence, ion_type, charge)
    except ValueError as e:
        st.error(e)
        st.stop()

    if delta_mass:
        st.warning('Ambiguous Modification! Using averagine Formula to es timate delta mass composition.')

        delta_mass_comp = pt.estimate_element_counts(delta_mass, pt.get_isotope_mods(sequence_input))
        combined_comp = {k: composition.get(k, 0) + delta_mass_comp.get(k, 0)
                         for k in set(composition) | set(delta_mass_comp)}
        composition = {k: round(v, 2) for k, v in combined_comp.items()}

    return composition


df['neutral_mass'] = df.apply(lambda x: pt.mass(x.sequence, ion_type=x.ion_type), axis=1)
df['mass_to_charge_ratio'] = df.apply(lambda x: pt.mz(x.sequence, charge=x.charge, ion_type=x.ion_type), axis=1)
df['composition'] = df.apply(lambda x: get_chem_composition(x.sequence, x.ion_type, x.charge), axis=1)

st.dataframe(df, use_container_width=True)

all_isotopes = []
for _, row in df.iterrows():
    isotopes = pt.isotopic_distribution(row['composition'],
                                        max_isotopes,
                                        min_abundance_threshold,
                                        distribution_resolution,
                                        use_neutron)
    sequence_mass = pt.mass(row['sequence'], ion_type=row['ion_type'], charge=row['charge'])

    if use_neutron:
        # add mass to the isotopes
        isotopes = [(sequence_mass + mass * pt.constants.NEUTRON_MASS, abundance) for mass, abundance in isotopes]

    # convert to m/z
    isotopes = [(mass/row['charge'],abundance) for mass, abundance in isotopes]

    # multiply by intensity
    isotopes = [(mass, abundance * row['intensity']) for mass, abundance in isotopes]

    all_isotopes.append(isotopes)

isotope_dict2 = []
for i, isotopes in enumerate(all_isotopes):
    for mz, abundance in isotopes:
        mz = round(mz, distribution_resolution)
        isotope_dict2.append({'mass_to_charge_ratio': mz, 'relative_abundance': abundance, 'sequence': sequences[i]})

df2 = pd.DataFrame(isotope_dict2)
st.dataframe(df2, use_container_width=True)

# make a bar chart
# make a bar chart
fig = go.Figure()

# Add bar chart which colors the bars by sequence (stacking mode)
for sequence in df2['sequence'].unique():
    df_temp = df2[df2['sequence'] == sequence]
    fig.add_trace(go.Bar(x=df_temp['mass_to_charge_ratio'], y=df_temp['relative_abundance'], name=sequence))

# Customize the plot for better readability
# Add barmode='stack' to stack the bars
fig.update_layout(title='Isotopic Distribution',
                  xaxis_title="Mass To Charge Ratio", yaxis_title="Relative Abundance",
                  margin=dict(l=40, r=40, t=40, b=40),
                  barmode='stack')  # This is where you set the bars to stack

st.plotly_chart(fig)

# merge all isotopes, if they have the same mass (rounded to distribution_resolution)
isotope_dict = {}
for isotopes in all_isotopes:
    for mz, abundance in isotopes:
        mz = round(mz, distribution_resolution)
        isotope_dict[mz] = isotope_dict.get(mz, 0) + abundance

# normalize the abundance so that the highest peak is 100%
total_abundance = sum(isotope_dict.values())
max_abundance = max(isotope_dict.values())
isotope_dict = {k: v / max_abundance for k, v in isotope_dict.items()}


df = pd.DataFrame(isotope_dict.items(), columns=['mass_to_charge_ratio', 'relative_abundance'])

# filter the df
df = df[df['relative_abundance'] > min_abundance_threshold]
df['relative_abundance'] = df['relative_abundance'] * 100

# make percent
df['relative_abundance'] = df['relative_abundance'].round(3)

# sort the df
df = df.sort_values(by='relative_abundance', ascending=False)

# filter df by max isotopes
df = df.head(max_isotopes)


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
df = df[['mass_to_charge_ratio', 'relative_abundance']]
df['percentage_abundance'] = df['relative_abundance'] / df['relative_abundance'].sum() * 100

# format the columns
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

filter_by = st.toggle(label='Sort by M/Z', value=False)

if filter_by:
    df = df.sort_values(by='mass_to_charge_ratio')
else:
    df = df.sort_values(by='relative_abundance', ascending=False)

st.markdown(df.style.hide(axis="index").to_html(), unsafe_allow_html=True)


