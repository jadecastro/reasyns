# This is a specification definition file for the LTLMoP toolkit.
# Format details are described at the beginning of each section below.


======== SETTINGS ========

Actions: # List of action propositions and their state (enabled = 1, disabled = 0)

CompileOptions:
convexify: True
parser: structured
symbolic: False
use_region_bit_encoding: False
synthesizer: jtlv
fastslow: True
decompose: True

CurrentConfigName:
dubins

Customs: # List of custom propositions

RegionFile: # Relative path of region description file
dubins.regions

Sensors: # List of sensor propositions and their state (enabled = 1, disabled = 0)

======== SPECIFICATION ========

RegionMapping: # Mapping between region names and their decomposed counterparts
A = p1
B = p2

Spec: # Specification in structured English
visit A
