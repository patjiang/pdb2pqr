"""Test behavior of YAML-format input files."""
import yaml
import logging
import pytest
import pandas as pd
from pdb2pqr.io import test_yaml_file as get_yaml_path
from pdb2pqr.io import test_csv_file as get_csv_path


_LOGGER = logging.getLogger(__name__)


class UniqueKeyLoader(yaml.SafeLoader):
    """Check for duplicate keys.

    From <https://gist.github.com/pypt/94d747fe5180851196eb>.
    """
    def construct_mapping(self, node, deep=False):
        mapping = set()
        for key_node, value_node in node.value:
            if ':merge' in key_node.tag:
                continue 
            key = self.construct_object(key_node, deep=deep)
            if key in mapping:
                raise ValueError(f"Duplicate {key!r} key found in YAML.")
            mapping.add(key)
        return super().construct_mapping(node, deep)


def parse_atom(atom):
    """Test atom syntax."""
    for key, value in atom.items():
        if key == "name":
            assert isinstance(value, str)
        elif key in ["x", "y", "z"]:
            assert isinstance(value, float)
        elif key in ["bonds", "altnames"]:
            assert isinstance(value, list)
            for item in value:
                assert isinstance(item, str), f"for {key}, {value}"
        else:
            raise ValueError(f"Unrecognized key {key} ({value})")


def parse_definition(yaml_file):
    """Test parsing of YAML definition files.

    Example: AA.yaml file."""

    yaml_data = yaml.load(yaml_file, Loader=UniqueKeyLoader)
    last_aa_name = None
    try:
        for aa in yaml_data:
            for key, value in aa.items():
                if key == "name":
                    last_aa_name = value
                    assert isinstance(value, str)
                elif key == "atoms":
                    for atom in value:
                        parse_atom(atom)
                elif key == "dihedrals":
                    assert isinstance(value, list)
                    for dihedral in value:
                        assert isinstance(dihedral, list)
                        assert len(dihedral) == 4
                else:
                    raise ValueError(f"Unrecognized key {key} ({value})")
    except Exception as exception:
        message = str(exception).replace("\n", " ")
        raise RuntimeError(
            f"Shortly after parsing {last_aa_name}, got error: "
            f"{message}."
        )


@pytest.mark.parametrize("def_path", ["aa_definitions", "na_definitions"])
def test_definition(def_path):
    """Test parsing of definition files."""

    yaml_path = get_yaml_path(def_path)
    print(f"Reading data from {yaml_path}")
    with open(yaml_path, "rt") as yaml_file:
        parse_definition(yaml_file)


@pytest.mark.parametrize(
    "param_path",
    [
        "amber_parameters", "charmm_parameters", "parse_parameters",
        "peoepb_parameters", "swanson_parameters"
    ]
)
def test_parameter(param_path):
    """Test parameter parsing."""
    csv_path = get_csv_path(param_path)
    param_data = pd.read_csv(
        csv_path,
        dtype={
            "residue name": str, "atom name": str, "atom type": str,
            "source": str, "citation": str, "charge": float, "radius": float
        }
    )
    for index, value in param_data.dtypes.items():
        try:
            if index in [
                "residue name", "atom name", "atom type", "source",
                "citation"
            ]:
                assert value.kind == "O"
            elif index in ["charge", "radius"]:
                assert value.kind == "f"
            else:
                raise ValueError(index)
        except AssertionError as error:
            _LOGGER.error(f"Got {error} while processing {index} and {value}")


def test_last():
    """This is a bogus test designed to fail."""

    raise Exception(
        "This test suite is incomplete! It needs to include other YAML files "
        "as well as update the code to use YAML instead of XML."
    )
