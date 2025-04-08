"""
Help message constants for the Peptide Isotope Calculator.
This file contains all help messages used throughout the application.
"""

# Input Type Help Messages
INPUT_TYPE_HELP = "Select the type of input you want to use: a peptide sequence, a chemical formula, or neutral mass value."

# Peptide Sequence Input Help
SEQUENCE_HELP = "Enter a peptide sequence using standard one-letter amino acid codes (e.g., PEPTIDE). Modifications are supported using ProForma 2.0 notation."

# Formula Input Help
FORMULA_HELP = "Enter a chemical formula using standard notation (e.g., C6H12O6). For isotopes, use bracket notation (e.g., [13C6]H12O6)."

# Mass Input Help
MASS_HELP = "Enter the neutral mass in Daltons. The calculator will estimate a chemical composition using averagine peptide composition."

# Ion Type Help
ION_TYPE_HELP = "Select the fragment ion type to calculate:\n- p: Precursor ion (intact peptide)\n- a/b/c: N-terminal fragments\n- x/y/z: C-terminal fragments\n- i: Immonium ion"

# Charge Help
CHARGE_HELP = "Enter the charge state of the ion (e.g., 2+). This affects the calculated m/z values but will not affect the isotopic distribution calculation."

# Intensity Help
INTENSITY_HELP = "Enter the abundance/intensity value for scaling your results. Higher values represent more abundant species."

# Max Isotopes Help
MAX_ISOTOPES_HELP = "Maximum number of isotope peaks to calculate and display. Higher values may be needed for larger molecules or higher charge states."

# Min Abundance Threshold Help
MIN_ABUNDANCE_HELP = "Minimum relative abundance threshold (%) for displaying isotope peaks. Peaks below this threshold will be excluded."

# Neutron Offset Help
NEUTRON_OFFSET_HELP = "When enabled, displays isotopic distribution as neutral mass ± N × neutron mass. When disabled, calculates the exact mass of each isotopologue."

# Resolution Help
RESOLUTION_HELP = "Controls the decimal precision of the calculated masses. Higher values provide more precise m/z values."

# Intensity Sum Help
INTENSITY_SUM_HELP = "When enabled, the sum of all isotopic peaks equals 100%. When disabled, the most abundant peak is set to 100% (normalized to the base peak)."

# Neutron Mass Help
NEUTRON_MASS_HELP = """Select the mass difference to use between isotope peaks:
- Neutron: Standard neutron mass (1.00866 Da)
- C13-C12: Mass difference between carbon isotopes (1.00335 Da)
- Averagine (Peptide): Average mass shift in peptides (1.00286 Da)
- Custom: Specify your own value"""

# Custom Neutron Mass Help
CUSTOM_NEUTRON_MASS_HELP = (
    "Enter a custom mass difference (in Da) to use between isotope peaks."
)

# Data Editor Help Messages
SEQUENCE_EDITOR_HELP = (
    "Enter peptide sequence. Must follow ProForma 2.0 notation for modifications."
)
CHARGE_EDITOR_HELP = "Specify the charge state of the peptide/molecule."
ION_TYPE_EDITOR_HELP = "Ion type for mass calculation: p (precursor), a/b/c (N-terminal), x/y/z (C-terminal), i (immonium)."
INTENSITY_EDITOR_HELP = "Relative abundance/intensity of this molecule in the plot."
INPUT_TYPE_EDITOR_HELP = (
    "Select whether this entry is a peptide sequence, formula, or mass value."
)

# Plot Help Messages
MZ_AXIS_HELP = "Mass-to-charge ratio (m/z) in Daltons/charge."
RELATIVE_ABUNDANCE_HELP = "Relative abundance of each isotope peak as a percentage."
ABSOLUTE_ABUNDANCE_HELP = "Absolute abundance values for each isotope peak."
