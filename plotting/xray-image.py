#! /bin/python3

import yt
from yt import derived_field

# from yt.fields.species_fields import add


@derived_field(name=("gas", "temperature"), units="K", sampling_type="cell")
def _temperature(field, data):
    return data.ds.arr(data["idefix-vtk", "TMP"], "K")


ds = yt.load(
    "p2r003-trc/data.0016.vtk",
    definitions_header="definitions.h",
    default_species_fields="ionized",
)

xray_fields = yt.add_xray_emissivity_field(
    ds, 0.4, 5.0, metallicity=0.15, table_type="apec"
)

for field in xray_fields:
    sp = yt.ProjectionPlot(ds, "x", field, origin="native")
    sp.save()
