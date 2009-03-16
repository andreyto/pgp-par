"""
Global variables for mass-spec analysis
"""

IsotopeWeights = {}

# Keys are ion type names, values are the corresponding ion instances
AllIonDict = {}

# Masses for amino acids (keys: 1-letter peptide abbreviations, like "G" or "D")
AminoMass = {}
AminoMassRight = {}

# List of all amino acid (left) masses:
AminoMasses = []

# Dictionary of post-translational modifications.  Keys are modification
# names (in lower-case).
PTMods = {}

# Truncated, 3- or 4-character keys:
PTModByShortName = {}
PTModList = []

AminoAcids = {} # key: single-letter abbreviation ("A" -> Alanine)
FixedMods = {"C":57.0518} # The protecting group on C is enabled, by default!

#List of  ModificationTypeObject (see Utils.py) the user defines as invivo or invitro
InVivoMods = []
InVitroMods = []

