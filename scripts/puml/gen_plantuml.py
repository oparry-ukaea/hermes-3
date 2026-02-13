import glob
import hpp2plantuml
import os
import plantuml

from hpp2plantuml.hpp2plantuml import MEMBER_PROP_MAP as access_type_map
from hpp2plantuml.hpp2plantuml import LINK_TYPE_MAP as link_type_map

this_dir = os.path.dirname(__file__)
repo_root = os.path.dirname(os.path.dirname(this_dir))
adas_output_base = "adas_reactions_classes"
amjuel_output_base = "amjuel_reactions_classes"
HCX_output_base = "HCX_reactions_classes"


def is_mem_var(l, mem_var_name):
    tmp = l.strip()
    if tmp and tmp[0] in access_type_map.values():
        return tmp[1:].startswith(mem_var_name)
    return False


def is_static(l):
    tmp = l.strip()
    static_pref = "{static}"
    return tmp and tmp[1 : 1 + len(static_pref)] == static_pref


def is_access(l, access_type):
    try:
        access_str = access_type_map[access_type]
    except KeyError:
        print(f"'{access_type}' not recognised as a valid access type")
        raise
    tmp = l.strip()
    return tmp.startswith(access_str)


def is_skipped_relationship(l, rel):
    tmp = l.split()
    if len(tmp) == 3 and tmp[1] in link_type_map.values():
        return tmp[0].endswith(rel) or tmp[2].endswith(rel)
    return False


def edit_diag(
    diag_str,
    skip_access=[],
    skip_member_vars=[],
    skip_relationships=[],
    skip_static=False,
    replace_strs={},
    preamble=[],
):
    sep = "\n"
    lines = diag_str.split(sep)
    # Filter out some lines
    access_excluded = lambda l: True in [is_access(l, access) for access in skip_access]
    mem_var_excluded = lambda l: True in [
        is_mem_var(l, mem_var) for mem_var in skip_member_vars
    ]
    relationship_excluded = lambda l: True in [
        is_skipped_relationship(l, rel) for rel in skip_relationships
    ]
    static_excluded = lambda l: is_static(l) if skip_static else False
    lines = [
        replace_strs.get(l, l)
        for l in lines
        if not access_excluded(l)
        and not mem_var_excluded(l)
        and not relationship_excluded(l)
        and not static_excluded(l)
    ]
    # for access in skip_access:
    #     lines = [line for line in lines if not is_access(line, access)]

    # Insert the preamble at the start
    start_idx = lines.index("@startuml")
    tmp_lines = lines[: start_idx + 1]
    tmp_lines.extend(preamble)
    tmp_lines.extend(lines[start_idx + 1 :])
    lines = tmp_lines

    diag_str = sep.join(lines)
    for find, rep in replace_strs.items():
        diag_str = diag_str.replace(find, rep)
    return diag_str


def gen_diag_str(header_paths, **diagram_kwargs):
    # Construct diagram object
    diagram_kwargs = {}
    diag = hpp2plantuml.Diagram(**diagram_kwargs)
    # Populate diagram with headers
    diag.create_from_file_list(header_paths)
    # Get diagram string
    return diag.render()


def get_this_dir_path():
    return os.path.realpath(os.path.dirname(__file__))


def write_puml(diag_str, output_path):
    with open(output_path, "wt") as fh:
        fh.write(diag_str)


def find_headers(patterns):
    header_paths = []
    for pattern in patterns:
        header_paths.extend(glob.glob(f"{repo_root}/include/{pattern}*.hxx"))
    return header_paths


def gen_reactions_diagrams(**edit_kws):
    adas_header_paths = find_headers(["adas"])
    gen_diagram_from_headers(
        adas_header_paths,
        adas_output_base,
        preamble=["skinparam wrapWidth 150"],
        **edit_kws,
    )
    amjuel_header_paths = find_headers(["amjuel", "reaction"])
    gen_diagram_from_headers(amjuel_header_paths, amjuel_output_base, **edit_kws)
    HCX_header_paths = find_headers(["hydrogen_charge_exchange", "reaction"])
    gen_diagram_from_headers(HCX_header_paths, HCX_output_base, **edit_kws)


def gen_diagram_from_headers(header_paths, output_base, **edit_kws):
    """Generate diagram string from header files"""
    diagram_kwargs = dict()
    diag_str = gen_diag_str(header_paths, **diagram_kwargs)
    diag_str = edit_diag(diag_str, **edit_kws)
    # Write the diagram string to file
    # Using the string directly with pl.processes doesn't seem to work...
    puml_path = os.path.join(get_this_dir_path(), f"{output_base}.puml")
    write_puml(diag_str, puml_path)

    gen_diagram_from_puml(puml_path, output_base)


def gen_diagram_from_puml(puml_path, output_base):
    # Generate the png from the diagram file
    pl = plantuml.PlantUML("http://www.plantuml.com/plantuml/img/")
    png_path = os.path.join(get_this_dir_path(), f"{output_base}.png")
    pl.processes_file(puml_path, outfile=png_path)
    print(f"Wrote {png_path}")


# gen_reactions_diagrams()

# gen_reactions_diagrams(
#     skip_member_vars=["calculate_rates", "electron_reaction"],
#     replace_strs={
#         "set_diagnostic_fields(Field3D& reaction_rate, Field3D& momentum_exchange, Field3D& energy_exchange, Field3D& energy_loss)": "set_diagnostic_fields(...)"
#     },
#     preamble=["skinparam wrapWidth 75"],
# )

for base in ["reactions_refactor_amjuel1", "reactions_refactor_amjuel2", "reactions_refactor_amjuel_subclasses"]:
# for base in [adas_output_base,amjuel_output_base,HCX_output_base]:
    gen_diagram_from_puml(
        f"{this_dir}/{base}.puml",
        f"{this_dir}/{base}",
    )
