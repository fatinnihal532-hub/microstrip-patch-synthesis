"""
Laminate catalogue used in the design-space study.

Values are the manufacturers' nominal datasheet figures (dielectric
constant and loss tangent near 10 GHz) and are treated as frequency
independent over the narrow band of a single patch. Real laminates have
lot-to-lot and frequency dispersion of a few percent in Dk; this study
uses the nominal values as given, not measured samples.
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Substrate:
    name: str
    eps_r: float
    tan_delta: float
    # Common catalogue thicknesses (mm). A board can only be bought in
    # these; confirm availability with the supplier before relying on one.
    standard_thicknesses_mm: tuple = ()


RT_DUROID_5880 = Substrate("RT/duroid 5880", 2.20, 0.0009,
                           (0.127, 0.254, 0.381, 0.508, 0.787, 1.575, 3.175))
RO4003C = Substrate("RO4003C", 3.55, 0.0027,       # Rogers' recommended *design* Dk
                    (0.203, 0.305, 0.406, 0.508, 0.813, 1.524))
FR4 = Substrate("FR-4", 4.40, 0.020,               # typical generic FR-4 figures
                (0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.2))
RT_DUROID_6010 = Substrate("RT/duroid 6010.2LM", 10.2, 0.0023,
                           (0.254, 0.635, 1.270, 1.905, 2.540))

CATALOGUE = (RT_DUROID_5880, RO4003C, FR4, RT_DUROID_6010)
