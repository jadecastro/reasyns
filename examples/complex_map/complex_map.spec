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
complex_map

Customs: # List of custom propositions

RegionFile: # Relative path of region description file
complex_map.regions

Sensors: # List of sensor propositions and their state (enabled = 1, disabled = 0)
crowded_B, 1
crowded_E, 1


======== SPECIFICATION ========

RegionMapping: # Mapping between region names and their decomposed counterparts
A = p1
C = p3
B = p2
E = p5
D = p4
G = p7
F = p6
H = p8

Spec: # Specification in structured English
visit A

