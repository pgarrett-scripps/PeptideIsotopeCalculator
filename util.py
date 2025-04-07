import streamlit as st
import streamlit_permalink as stp
import peptacular as pt
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import requests

from dataclasses import dataclass
from typing import Literal, Optional
from urllib.parse import quote_plus

from constants import *


@dataclass
class SingleIsoInput:
    """
    Class to hold the input values for a single isotope.
    """

    input_type: Literal["Peptide", "Formula", "Mass"]
    sequence: str
    ion_type: Optional[Literal["a", "b", "c", "x", "y", "z", "p", "i"]]
    charge: int
    max_isotopes: int
    min_abundance_threshold: float
    use_neutron: bool
    distribution_resolution: int
    neutron_value: float
    intensity: float = 1
    is_intensity_sum: bool = False

    @property
    def is_peptide(self) -> bool:
        """
        Check if the input type is peptide.
        """
        return self.input_type == "Peptide"
    
    @property
    def is_formula(self) -> bool:
        """
        Check if the input type is formula.
        """
        return self.input_type == "Formula"
    
    @property
    def is_mass(self) -> bool:
        """
        Check if the input type is mass.
        """
        return self.input_type == "Mass"
    
    @property
    def composition(self) -> dict:
        """
        Get the composition of the input.
        """
        if self.is_peptide:
            return pt.comp_mass(self.sequence, self.ion_type)[0]        
        if self.is_formula:
            return pt.parse_chem_formula(self.sequence)
        if self.is_mass:
            return pt.estimate_comp(float(self.sequence))
        
        raise ValueError("Invalid input type")
    
    @property
    def chemical_formula(self) -> str:
        """
        Get the chemical formula of the input.
        """
        if self.is_peptide:
            return pt.write_chem_formula(self.composition, hill_order=True, precision=1)
        if self.is_formula:
            return pt.write_chem_formula(self.composition, precision=1)
        if self.is_mass:
            return pt.write_chem_formula(self.composition, hill_order=True, precision=1)
        
        raise ValueError("Invalid input type")
    
    @property
    def neutral_mass(self) -> float:
        """
        Get the neutral mass of the input.
        """
        if self.is_peptide:
            return pt.mass(self.sequence, ion_type=self.ion_type)
        if self.is_formula:
            return pt.chem_mass(self.sequence)
        if self.is_mass:
            return float(self.sequence)
        
        raise ValueError("Invalid input type")
    
    @property
    def mass(self) -> float:
        """
        Get the mass of the input.
        """
        if self.is_peptide:
            return pt.mass(self.sequence, ion_type=self.ion_type, charge=self.charge)
        if self.is_formula:
            return pt.chem_mass(self.sequence) + self.charge * pt.PROTON_MASS
        if self.is_mass:
            return float(self.sequence) + self.charge * pt.PROTON_MASS
        
        raise ValueError("Invalid input type")


    @property
    def mz(self) -> float:
        """
        Get the m/z of the input.
        """
        if self.is_peptide:
            return (self.neutral_mass + self.charge * pt.PROTON_MASS) / self.charge
        if self.is_formula or self.is_mass:
            return (self.neutral_mass + self.charge * pt.PROTON_MASS) / self.charge
        
        raise ValueError("Invalid input type")
    
    
    @property
    def isotopes(self):
        # Calculate isotopic distributions based on input type
        if self.is_peptide or self.is_formula:
            isotopes = pt.isotopic_distribution(
                self.composition,
                self.max_isotopes,
                self.min_abundance_threshold,
                self.distribution_resolution,
                self.use_neutron,
                neutron_mass=self.neutron_value,
                output_masses_for_neutron_offset=True,
                distribution_abundance=self.intensity,
                is_abundance_sum=self.is_intensity_sum,
            )
        elif self.is_mass:
            isotopes = pt.estimate_isotopic_distribution(
                float(self.sequence),
                self.max_isotopes,
                self.min_abundance_threshold,
                self.distribution_resolution,
                self.use_neutron,
                neutron_mass=self.neutron_value,
                output_masses_for_neutron_offset=True,
                distribution_abundance=self.intensity,
                is_abundance_sum=self.is_intensity_sum,
            )
        else:
            raise ValueError("Invalid input type")
        
        # add delta mass
        if self.delta_mass != 0:
            isotopes = [
                (mass + self.delta_mass, abundance)
                for mass, abundance in isotopes
            ]

        return isotopes


    @property
    def delta_mass(self) -> float:
        if self.is_peptide:
            _, delta_mass = pt.comp_mass(self.sequence, self.ion_type)
            return delta_mass
        if self.is_formula:
            return 0
        if self.is_mass:
            return 0
        
        raise ValueError("Invalid input type")

def validate_input(self):

    
    if self.is_peptide and not self.sequence:
        st.warning("Please enter a sequence.")
        st.stop()
    
    if len(self.isotopes) == 0:
        st.error("No isotopes found.")
        st.stop()

    if self.mass >= MAX_MASS:
        st.error("The mass is too high. Please check the input and try again.")
        st.stop()

    if len(self.composition) >= 10:
        st.error(
            "The formula is too complex (Limit of 10 unique elements). Please check the input and try again."
        )
        st.stop()

    if self.delta_mass != 0:
        st.warning(
            "Ambiguous Modification! Using averagine Formula to estimate delta mass composition."
        )

    # check if composition is has floats
    if self.is_peptide or self.is_formula:
        for k, v in self.composition.items():
            if isinstance(v, float):
                st.warning(
                    f"Composition has floats: {k}: {v}. This will be rounded to: {round(v)} before calculating the isotopic distribution."
                )


@dataclass
class MultiIsoInput:
    single_iso_inputs: list[SingleIsoInput]

    @property
    def is_intensity_sum(self) -> bool:
        """
        Check if the input type is intensity sum.
        """
        return self.single_iso_inputs[0].is_intensity_sum

            
def get_input_settings() -> tuple:
    
    c1, c2 = st.columns(2)
    with c1:
        max_isotopes = stp.number_input(
            label="Max Number of Isotopes",
            value=DEFAULT_ISOTOPES,
            min_value=MIN_ISOTOPES,
            max_value=MAX_ISOTOPES,
            step=ISOTOPES_STEP,
            key="max_isotopes",
            help="Maximum number of isotopes to display",
        )

        use_neutron = stp.toggle(
            label="Neutron Offset",
            value=DEFAULT_USE_NEUTRON,
            key="use_neutron",
            help="Use neutron offsets for isotopic distribution instead of mass of each isotope. "
            "This gives isotopic peaks at Neutral Mass +/- N * Neutron Mass",
        )

    with c2:
        min_abundance_threshold = stp.selectbox(
            label="Min Relative Abundance",
            options=ABUNDANCE_OPTIONS,
            index=DEFAULT_ABUNDANCE_INDEX,
            help="Minimum abundance threshold to display",
        )

        is_intensity_sum = stp.toggle(
                label="Intensity Sum",
                value=DEFAULT_IS_INTENSITY_SUM,
                key="is_intensity_sum",
                help="If selected, then the sum of all isotopic peaks should equal to 100%. If false. the largest peak will be set to 100%.",
            )
        
    distribution_resolution = DEFAULT_DISTRIBUTION_RESOLUTION
    neutron_value = None
    if not use_neutron:
        distribution_resolution = stp.slider(
            label="Resolution",
            value=DEFAULT_DISTRIBUTION_RESOLUTION,
            min_value=MIN_DISTRIBUTION_RESOLUTION,
            max_value=MAX_DISTRIBUTION_RESOLUTION,
            step=DISTRIBUTION_RESOLUTION_STEP,
            help="Resolution of the distribution (Round to nearest X decimal places)",
        )
        
    neutron_option = stp.selectbox(
        label="Neutron Mass",
        options=NEUTRON_MASS_OPTIONS,
        index=DEFAULT_NEUTRON_MASS_INDEX,
        help="Neutron mass to use for isotopic distribution calculation",
    )

    if neutron_option == "Neutron":
        neutron_value = NEUTRON_MASS_VALUE
    elif neutron_option == "C13-C12":
        neutron_value = C13_C12_DIFF
    elif neutron_option == "Averagine (Peptide)":
        neutron_value = AVERAGINE_MASS
    elif neutron_option == "Custom":
        neutron_value = stp.number_input(
            label="Custom Neutron Mass",
            value=NEUTRON_MASS_VALUE,
            format="%.6f",
            help="Custom neutron mass to use for isotopic distribution calculation",
        )

    st.caption(f"Neutron Mass: {neutron_value:.6f} Da")

    return max_isotopes, min_abundance_threshold, use_neutron, distribution_resolution, is_intensity_sum, neutron_value

def get_single_app_input() -> SingleIsoInput:
    """
    Get a single app input from the user and return a SingleIsoInput object.
    No validation is performed on the input.
    """
    input_type = stp.radio(
        label="Input Type",
        options=SEQUENCE_INPUT_OPTIONS,
        index=DEFAULT_INPUT_INDEX,
        horizontal=True,
        key="input_type",
    )

    sequence_input = ""
    ion_type = None
    neutron_value = 1.002856  # Default to Averagine

    # Handle different input types
    if input_type == "Peptide":
        sequence_input = stp.text_input(
            label="Sequence",
            value=DEFAULT_SEQUENCE,
            help="Enter the peptide sequence (e.g. PEPTIDE)",
            key="sequence_input",
        )

        ion_type = stp.radio(
            label="Ion Type",
            options=ION_TYPES,
            index=DEFAULT_ION_TYPE_INDEX,
            horizontal=True,
            help="Ion type to calculate, p = precursor",
            key="ion_type",
        )

    elif input_type == "Formula":
        sequence_input = stp.text_input(
            label="Formula",
            value=DEFAULT_FORMULA,
            help="Enter the chemical formula (e.g. C6H12O6). Isotopes must be enclosed in "
            "brackets (e.g. C[13]H12O6)",
            key="formula_input",
        )

    elif input_type == "Mass":
        sequence_input = str(stp.number_input(
            label="Neutral (Monoisotopic) Mass",
            value=DEFAULT_MASS,
            max_value=MAX_MASS,
            min_value=MIN_MASS,
            step=MASS_STEP,
            key="mass_input",
            format="%.4f",
            help="Enter the neutral mass (e.g. 1000 Da)",
        ))

    # Common inputs for all types
    charge = stp.number_input(
        label="Charge",
        value=DEFAULT_CHARGE,
        min_value=MIN_CHARGE,
        max_value=MAX_CHARGE,
        step=CHARGE_STEP,
        key="charge",
        help="Enter the charge state (e.g. 2+)",
    )

    # Get other input settings
    max_isotopes, min_abundance_threshold, use_neutron, distribution_resolution, is_intensity_sum, neutron_value = get_input_settings()


    return SingleIsoInput(
        input_type=input_type,
        sequence=sequence_input,
        ion_type=ion_type,
        charge=charge,
        max_isotopes=max_isotopes,
        min_abundance_threshold=min_abundance_threshold,
        use_neutron=use_neutron,
        distribution_resolution=distribution_resolution,
        neutron_value=neutron_value,
        intensity=1,
        is_intensity_sum=is_intensity_sum,
    )


def get_multi_app_input() -> MultiIsoInput:
    """
    Get a multi app input from the user and return a MultiIsoInput object.
    No validation is performed on the input.
    """
    
    default_df = pd.DataFrame(
        {
            "sequence": ["PEPTIDE", "PEPTIDEPEPTIDE"],
            "charge": [2, 4],
            "ion_type": ["p", "p"],
            "intensity": [1000, 2000],
            "input_type": ["Peptide", "Peptide"],
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
                options=["Peptide", "Formula", "Mass"],
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

    max_isotopes, min_abundance_threshold, use_neutron, distribution_resolution, is_intensity_sum, neutron_value = get_input_settings()

    # create a list of SingleIsoInput objects
    single_iso_inputs = []
    for _, row in df.iterrows():
        single_iso_inputs.append(
            SingleIsoInput(
                input_type=row["input_type"],
                sequence=row["sequence"],
                ion_type=row["ion_type"],
                charge=row["charge"],
                max_isotopes=max_isotopes,
                min_abundance_threshold=min_abundance_threshold,
                use_neutron=use_neutron,
                distribution_resolution=distribution_resolution,
                neutron_value=neutron_value,  # Default to Averagine
                intensity=row["intensity"],
                is_intensity_sum=is_intensity_sum,
            )
        )

    return MultiIsoInput(single_iso_inputs=single_iso_inputs)



def add_dumb_mobile_buffer():
    # blank space for mobile cause streamlit is weird
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)
    


def construct_isotope_df(params: SingleIsoInput) -> pd.DataFrame:
    """
    Construct a DataFrame of isotopes from the parameters.
    """
    df = pd.DataFrame(params.isotopes, columns=["neutral_mass", "abundance"])
    
    # filter the df
    df = df[df["abundance"] > params.min_abundance_threshold]

    # sort the df
    df = df.sort_values(by="abundance", ascending=False)

    # apply the charge state
    df["mass_to_charge_ratio"] = (
        (df["neutral_mass"] + params.charge * pt.PROTON_MASS) / params.charge
    )

    return df

def construct_figure(df: pd.DataFrame, params: SingleIsoInput) -> go.Figure:
    """
    Construct a Plotly figure for the isotopes.
    """
    # Assuming `df` is your DataFrame with isotopes
    fig = go.Figure()

    # Add lines from each isotope marker to the base
    for idx, row in df.iterrows():
        fig.add_trace(
            go.Scatter(
                x=[row['mass_to_charge_ratio'], row['mass_to_charge_ratio']],
                y=[0, row["abundance"]*100],
                mode="lines",
                line=dict(color="grey", width=3),
                showlegend=False,
            )
        )

    
    # Add scatter plot for the isotopes markers
    fig.add_trace(
        go.Scatter(
            x=df['mass_to_charge_ratio'],
            y=df['abundance']*100,
            mode="markers",
            marker=dict(size=6),
            name="Isotopes",
        )
    )


    # Customize the plot for better readability
    fig.update_layout(
        title=f'Isotopic Distribution: {params.sequence} {"da" if params.input_type == "Mass" else ""}',
        xaxis_title="Mass To Charge Ratio",
        yaxis_title="Realtive Abundance (%)",
        margin=dict(l=40, r=40, t=40, b=40),
    )

    return fig


def construct_multi_isotope_figure(df: pd.DataFrame) -> go.Figure:
    """
    Construct a Plotly figure for multiple isotope distributions where each sequence
    is represented by a different color, and isotopes with the same m/z are stacked.
    
    Args:
        df: DataFrame containing isotope data with 'sequence', 'mass_to_charge_ratio',
            and 'abundance' columns.
        
    Returns:
        Plotly Figure object with colored and stacked isotope distributions.
    """
    # Create a new figure
    fig = go.Figure()
    
    # Get unique sequences for coloring
    sequences = df['sequence'].unique()
    
    # Create a color map for the sequences
    colormap = plt.cm.get_cmap('tab10', len(sequences))
    sequence_colors = {seq: f'rgba({int(255*colormap(i)[0])}, {int(255*colormap(i)[1])}, {int(255*colormap(i)[2])}, 0.8)'
                      for i, seq in enumerate(sequences)}
    
    # Group by m/z values to stack isotopes with the same m/z
    grouped_df = df.groupby('mass_to_charge_ratio')
    
    # Track which sequences have already been added to the legend
    legend_shown = set()
    
    # Find maximum abundance sum for y-axis scaling
    max_total_abundance = 0
    
    # Process each m/z group
    for mz, group in grouped_df:
        # Sort by sequence to ensure consistent stacking order
        group = group.sort_values('sequence')
        
        # Stack heights for this m/z
        current_height = 0
        
        for _, row in group.iterrows():
            sequence = row['sequence']
            abundance = row['abundance']
            
            # Add stacked bar for this sequence
            fig.add_trace(
                go.Scatter(
                    x=[mz, mz],
                    y=[current_height, current_height + abundance],
                    mode="lines",
                    line=dict(color=sequence_colors[sequence], width=3),
                    name=sequence,
                    legendgroup=sequence,
                    showlegend=(sequence not in legend_shown)  # Show legend only once per sequence
                )
            )
            
            # Add marker at the top of each bar segment
            fig.add_trace(
                go.Scatter(
                    x=[mz],
                    y=[current_height + abundance],
                    mode="markers",
                    marker=dict(
                        size=6,
                        color=sequence_colors[sequence]
                    ),
                    legendgroup=sequence,
                    showlegend=False,
                )
            )
            
            # Update current height for stacking
            current_height += abundance
            
            # Mark this sequence as shown in the legend
            legend_shown.add(sequence)
        
        # Update the maximum total abundance
        max_total_abundance = max(max_total_abundance, current_height)
    
    # Customize the plot for better readability
    fig.update_layout(
        title="Isotopic Distribution",
        xaxis_title="Mass To Charge Ratio",
        yaxis=dict(
            title="Abundance",
            range=[0, max_total_abundance * 1.05],
            showgrid=True,
        ),
        margin=dict(l=40, r=40, t=40, b=40),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right", 
            x=0.99
        )
    )
    
    return fig


def listify(o=None):
    if o is None:
        res = []
    elif isinstance(o, list):
        res = o
    elif isinstance(o, str):
        res = [o]
    else:
        res = [o]
    return res


def shorten_url(url: str) -> str:
    """Shorten a URL using TinyURL."""
    api_url = f"http://tinyurl.com/api-create.php?url={url}"
    
    try:
        response = requests.get(api_url)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        return f"Error: {e}"
    

def get_query_params_url(params_dict):
    """
    Create url params from alist of parameters and a dictionary with values.

    Args:
        params_list (str) :
            A list of parameters to get the value of from `params_dict`
        parmas_dict (dict) :
            A dict with values for the `parmas_list .
        **kwargs :
            Extra keyword args to add to the url
    """
    return "?" + "&".join(
        [
            f"{key}={quote_plus(str(value))}"
            for key, values in params_dict.items()
            for value in listify(values)
        ]
    )