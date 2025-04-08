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
from help_messages import *


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
    line_width: int = 3
    is_log: bool = False

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
                chemical_formula=self.composition,
                max_isotopes=self.max_isotopes,
                min_abundance_threshold=self.min_abundance_threshold,
                distribution_resolution=self.distribution_resolution,
                use_neutron_count=self.use_neutron,
                conv_min_abundance_threshold=None,
                distribution_abundance=self.intensity,
                is_abundance_sum=self.is_intensity_sum,
                output_masses_for_neutron_offset=True,
                neutron_mass=self.neutron_value,

            )
        elif self.is_mass:
            isotopes = pt.estimate_isotopic_distribution(
                neutral_mass=float(self.sequence),
                max_isotopes=self.max_isotopes,
                min_abundance_threshold=self.min_abundance_threshold,
                distribution_resolution=self.distribution_resolution,
                use_neutron_count=self.use_neutron,
                conv_min_abundance_threshold=None,
                distribution_abundance=self.intensity,
                is_abundance_sum=self.is_intensity_sum,
                output_masses_for_neutron_offset=True,
                neutron_mass=self.neutron_value,
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
    
    @property
    def line_width(self) -> int:
        """
        Get the line width of the input.
        """
        return self.single_iso_inputs[0].line_width
    
    @property
    def is_log(self) -> bool:
        """
        Check if the input type is log.
        """
        return self.single_iso_inputs[0].is_log

            
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
            help=MAX_ISOTOPES_HELP,
        )

        use_neutron = stp.toggle(
            label="Neutron Offset",
            value=DEFAULT_USE_NEUTRON,
            key="use_neutron",
            help=NEUTRON_OFFSET_HELP,
        )

    with c2:
        min_abundance_threshold = stp.selectbox(
            label="Min Relative Abundance",
            options=ABUNDANCE_OPTIONS,
            index=DEFAULT_ABUNDANCE_INDEX,
            help=MIN_ABUNDANCE_HELP,
        )

        is_intensity_sum = stp.toggle(
                label="Intensity Sum",
                value=DEFAULT_IS_INTENSITY_SUM,
                key="is_intensity_sum",
                help=INTENSITY_SUM_HELP)
        
    distribution_resolution = DEFAULT_DISTRIBUTION_RESOLUTION
    neutron_value = None
    if not use_neutron:
        distribution_resolution = stp.slider(
            label="Resolution",
            value=DEFAULT_DISTRIBUTION_RESOLUTION,
            min_value=MIN_DISTRIBUTION_RESOLUTION,
            max_value=MAX_DISTRIBUTION_RESOLUTION,
            step=DISTRIBUTION_RESOLUTION_STEP,
            help=RESOLUTION_HELP,
            key="distribution_resolution",
        )
        
    neutron_option = stp.selectbox(
        label="Neutron Mass",
        options=NEUTRON_MASS_OPTIONS,
        index=DEFAULT_NEUTRON_MASS_INDEX,
        help=NEUTRON_MASS_HELP,
        key="neutron_mass_option",
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
            help=CUSTOM_NEUTRON_MASS_HELP,
            key="custom_neutron_mass",
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
        help=INPUT_TYPE_HELP,
    )

    sequence_input = ""
    ion_type = None
    neutron_value = 1.002856  # Default to Averagine

    # Handle different input types
    if input_type == "Peptide":
        sequence_input = stp.text_input(
            label="Sequence",
            value=DEFAULT_SEQUENCE,
            placeholder=DEFAULT_SEQUENCE,
            help=SEQUENCE_HELP,
            key="sequence_input",
            
        )

        ion_type = stp.radio(
            label="Ion Type",
            options=ION_TYPES,
            index=DEFAULT_ION_TYPE_INDEX,
            horizontal=True,
            help=ION_TYPE_HELP,
            key="ion_type",
        )

    elif input_type == "Formula":
        sequence_input = stp.text_input(
            label="Chemical Formula",
            value=DEFAULT_FORMULA,
            help=FORMULA_HELP,
            placeholder=DEFAULT_FORMULA,
            key="formula_input",
        )

    elif input_type == "Mass":
        sequence_input = str(stp.number_input(
            label="Neutral Mass",
            value=DEFAULT_MASS,
            max_value=MAX_MASS,
            min_value=MIN_MASS,
            step=MASS_STEP,
            key="mass_input",
            format="%.4f",
            help=MASS_HELP,
        ))

    c1, c2 = st.columns(2)
    with c1:
        # Common inputs for all types
        charge = stp.number_input(
            label="Charge",
            value=DEFAULT_CHARGE,
            min_value=MIN_CHARGE,
            max_value=MAX_CHARGE,
            step=CHARGE_STEP,
            key="charge",
            help=CHARGE_HELP,
        )

    with c2:
        intensity = stp.number_input(
            label="intensity",
            value=100.0,
            min_value=0.0,
            key="intensity",
            help=INTENSITY_HELP,
        )


    with st.expander("Advanced Options", expanded=False):
        # Get other input settings
        max_isotopes, min_abundance_threshold, use_neutron, distribution_resolution, is_intensity_sum, neutron_value = get_input_settings()
        line_width = stp.number_input(
            "Line Width",
            min_value=1,
            max_value=10,
            value=3,
            step=1,
            help="Set the line width for the plot.",
            key="line_width")
        
        is_log = stp.toggle(
            label="Log Scale",
            value=False,
            key="is_log",
            help="Use log scale for the y-axis.",
        )


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
        intensity=intensity,
        is_intensity_sum=is_intensity_sum,
        line_width=line_width,
        is_log=is_log
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
                help=SEQUENCE_HELP,
                width="small",
                required=True,
            ),
            "charge": st.column_config.NumberColumn(
                "Charge", default=2, help=CHARGE_HELP, width="small"
            ),
            "ion_type": st.column_config.SelectboxColumn(
                "Ion Type",
                help=ION_TYPE_HELP,
                default="p",
                options=list("abcxyzpi"),
                width="small",
                required=True,
            ),
            "intensity": st.column_config.NumberColumn(
                "Intensity",
                default=1000,
                help=INTENSITY_HELP,
                width="small",
                required=True,
            ),
            "input_type": st.column_config.SelectboxColumn(
                "Input Type",
                default="peptide",
                help=INPUT_TYPE_HELP,
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

    with st.expander("Advanced Options", expanded=False):
        # Get other input settings
        max_isotopes, min_abundance_threshold, use_neutron, distribution_resolution, is_intensity_sum, neutron_value = get_input_settings()

        line_width = stp.number_input(
            "Line Width",
            min_value=1,
            max_value=10,
            value=3,
            step=1,
            help="Set the line width for the plot.",
            key="line_width",
        )

        is_log = stp.toggle(
            label="Log Scale",
            value=False,
            key="is_log",
            help="Use log scale for the y-axis.",
        )

    

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
                line_width=line_width,
                is_log=is_log,
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

    # sort the df
    df = df.sort_values(by="abundance", ascending=False)

    # realtive abundance
    if params.is_intensity_sum:
        df["reative_abundance"] = df["abundance"] / df["abundance"].sum()
    else:
        df["reative_abundance"] = df["abundance"] / df["abundance"].max()

    # apply the charge state
    df["mz"] = (
        (df["neutral_mass"] + params.charge * pt.PROTON_MASS) / params.charge
    )

    return df

def construct_figure(df: pd.DataFrame, params: SingleIsoInput) -> go.Figure:
    """
    Construct a Plotly figure for the isotopes with both relative and absolute abundance axes.
    """
    # Create figure with secondary y-axis
    fig = go.Figure()

    # Add lines from each isotope marker to the base
    for idx, row in df.iterrows():
        fig.add_trace(
            go.Scatter(
                x=[row['mz'], row['mz']],
                y=[0, row["reative_abundance"]*100],
                mode="lines",
                line=dict(color="grey", width=params.line_width),
                showlegend=False,
                yaxis="y"
            )
        )
    
    # Add scatter plot for the isotopes markers with relative abundance
    fig.add_trace(
        go.Scatter(
            x=df['mz'],
            y=df['reative_abundance']*100,
            mode="markers",
            marker=dict(size=params.line_width*2, color="grey"),
            name="Relative Abundance (%)",
            yaxis="y",
            showlegend=False,

        )
    )
    
    # Add scatter plot for absolute abundance values on secondary axis
    fig.add_trace(
        go.Scatter(
            x=df['mz'],
            y=df['abundance'],
            mode="markers",
            marker=dict(size=params.line_width*2, color="grey", opacity=0.6),
            name="Abundance",
            yaxis="y2",
            #dont show legend
            showlegend=False,

        )
    )

    y2_color = "rgba(77, 150, 214, 0.77)"
    y_color = "rgba(55, 55, 55, 0.37)"

    # Customize the plot with dual y-axes
    fig.update_layout(
        xaxis_title="m/z",
        yaxis=dict(
            title="Relative Abundance (%)",
            titlefont=dict(color="grey"),
            tickfont=dict(color="grey"),
            # change grid line color
            gridcolor="rgba(55, 55, 55, 0.17)",
            rangemode="nonnegative",  # Ensure non-negative range starting at 0
            #make log
            type="log" if params.is_log else "linear",

        ),
        yaxis2=dict(
            title="Abundance",
            titlefont=dict(color="rgba(77, 150, 214, 0.87)"),
            tickfont=dict(color="rgba(77, 150, 214, 0.87)"),
            showgrid=False,
            # change grid line colro
            gridcolor="rgba(77, 150, 214, 0.27)",
            anchor="x",
            overlaying="y",
            side="right",
            rangemode="nonnegative",  # Ensure non-negative range starting at 0
            type="log" if params.is_log else "linear",

        ),
        margin=dict(l=40, r=60, t=40, b=40),
        legend=dict(            
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )

    return fig


def construct_multi_isotope_figure(df: pd.DataFrame, line_width=3, is_log=False) -> go.Figure:
    """
    Construct a Plotly figure for multiple isotope distributions where each sequence
    is represented by a different color, and isotopes with the same m/z are stacked.
    
    Args:
        df: DataFrame containing isotope data with 'sequence', 'mass_to_charge_ratio',
            and 'abundance' columns.
        
    Returns:
        Plotly Figure object with colored and stacked isotope distributions with dual y-axes.
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
    grouped_df = df.groupby('mz')
    
    # Track which sequences have already been added to the legend
    legend_shown = set()
    
    # Find maximum abundance sum for y-axis scaling
    max_total_abundance = 0
    
    # Calculate relative abundances (primary y-axis)
    total_abundance_per_mz = {}
    max_abundance_overall = df['abundance'].max()
    
    # First calculate total abundance per m/z group
    for mz, group in grouped_df:
        total_abundance = group['abundance'].sum()
        total_abundance_per_mz[mz] = total_abundance
        max_total_abundance = max(max_total_abundance, total_abundance)
    
    # Process each m/z group
    for mz, group in grouped_df:
        # Sort by sequence to ensure consistent stacking order
        group = group.sort_values('sequence')
        
        # Stack heights for this m/z
        current_height_abs = 0
        current_height_rel = 0
        
        for _, row in group.iterrows():
            sequence = row['sequence']
            abundance = row['abundance']
            
            # Calculate relative abundance (0-100%)
            relative_abundance = row['relative_abundance'] * 100
            
            # Add stacked bar for absolute abundance (secondary y-axis)
            fig.add_trace(
                go.Scatter(
                    x=[mz, mz],
                    y=[current_height_abs, current_height_abs + abundance],
                    mode="lines",
                    line=dict(color=sequence_colors[sequence], width=line_width),
                    name=sequence,
                    legendgroup=sequence,
                    showlegend=(sequence not in legend_shown),  # Show legend only once per sequence
                    yaxis="y2"
                )
            )
            
            # Add marker at the top of each absolute abundance bar segment
            fig.add_trace(
                go.Scatter(
                    x=[mz],
                    y=[current_height_abs + abundance],
                    mode="markers",
                    marker=dict(
                        size=line_width*2,
                        color=sequence_colors[sequence]
                    ),
                    legendgroup=sequence,
                    showlegend=False,
                    yaxis="y2"
                )
            )
            
            # Add stacked bar for relative abundance (primary y-axis)
            fig.add_trace(
                go.Scatter(
                    x=[mz, mz],
                    y=[current_height_rel, current_height_rel + relative_abundance],
                    mode="lines",
                    line=dict(color=sequence_colors[sequence], width=0, dash="dot"),
                    legendgroup=sequence,
                    showlegend=False,
                    yaxis="y"
                )
            )
            
            # Update current heights for stacking
            current_height_abs += abundance
            current_height_rel += relative_abundance
            
            # Mark this sequence as shown in the legend
            legend_shown.add(sequence)
    
    y_color = "rgba(55, 55, 55, 0.77)"
    y2_color = "rgba(77, 150, 214, 0.77)"
    
    # Customize the plot with dual y-axes
    fig.update_layout(
        title="Isotopic Distribution",
        xaxis_title="m/z",
        yaxis=dict(
            title="Relative Abundance (%)",
            titlefont=dict(color=y_color),
            tickfont=dict(color=y_color),
            gridcolor="rgba(55, 55, 55, 0.17)",
            rangemode="nonnegative",  # Ensure non-negative range starting at 0
            zeroline=True,           # Show zero line
            type="log" if is_log else "linear",

        ),
        yaxis2=dict(
            title="Absolute Abundance",
            titlefont=dict(color=y2_color),
            tickfont=dict(color=y2_color),
            anchor="x",
            overlaying="y",
            side="right",
            showgrid=False,
            gridcolor="rgba(77, 150, 214, 0.27)",
            rangemode="nonnegative",  # Ensure non-negative range starting at 0
            zeroline=True,           # Show zero line
            # is log
            type="log" if is_log else "linear",
        ),
        margin=dict(l=40, r=60, t=40, b=40),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    # Force both axes to include zero in their range
    fig.update_yaxes(rangemode="nonnegative", constrain="domain")
    
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