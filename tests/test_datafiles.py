"""Test behavior of YAML-format input files."""
import yaml
import logging
import pytest
import re
import pandas as pd
from pdb2pqr.io import test_yaml_file as get_yaml_path
from pdb2pqr.io import test_csv_file as get_csv_path


_LOGGER = logging.getLogger(__name__)
OPT_TYPES = {"Flip", "Carboxylic", "Alcoholic", "Water", "Generic"}


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


def parse_definition_atom(atom):
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
    """Test parsing of YAML definition files."""

    yaml_data = yaml.load(yaml_file, Loader=UniqueKeyLoader)
    last_residue_name = None
    try:
        for residue in yaml_data:
            for key, value in residue.items():
                if key == "name":
                    last_residue_name = value
                    assert isinstance(value, str)
                elif key == "atoms":
                    for atom in value:
                        parse_definition_atom(atom)
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
            f"Shortly after parsing {last_residue_name}, got error: "
            f"{message}."
        )


def parse_name_atom(atom):
    """Parse atom entry from names file."""
    for key, value in atom.items():
        if key == "original name":
            re.compile(value)
        elif key == "new name":
            assert isinstance(value, str)
        else:
            raise ValueError(f"Unknown entry: {key}, {value}")


def parse_name(yaml_file):
    """Test parsing of YAML renaming files."""
    yaml_data = yaml.load(yaml_file, Loader=UniqueKeyLoader)
    for pattern in yaml_data:
        for key, value in pattern.items():
            if key == "pattern":
                re.compile(value)
            elif key == "atoms":
                for atom in value:
                    parse_name_atom(atom)
            elif key in ["new residue name", "exclude"]:
                assert isinstance(value, str)
            else:
                raise ValueError(f"Unknown entry: {key}: {value}")


def parse_patch_atom(atom):
    """Test parsing of patch file atoms."""
    for key, value in atom.items():
        if key == "atom name":
            assert isinstance(value, str)
        elif key in ["x", "y", "z"]:
            assert isinstance(value, float)
        elif key in ["altnames", "bonds"]:
            assert isinstance(value, list)
            for bond in value:
                assert isinstance(bond, str)
        elif key == "dihedrals":
            for atom_list in value:
                assert len(atom_list) == 4
                for atom in atom_list:
                    assert isinstance(atom, str)
        else:
            raise ValueError(f"Invalid entry: {key} {value}")


@pytest.mark.parametrize("path", ["patches.yaml"])
def test_patches(path):
    """Test parsing of YAML patch files."""
    yaml_path = get_yaml_path(path)
    with open(yaml_path, "rt") as yaml_file:
        yaml_data = yaml.load(yaml_file, Loader=UniqueKeyLoader)
        for patch in yaml_data:
            for key, value in patch.items():
                # TODO - pattern is not currently a proper regex
                if key in ["name", "pattern", "new name"]:
                    assert isinstance(key, str)
                elif key in ["add atoms", "remove atoms"]:
                    for atom in value:
                        parse_patch_atom(atom)
                else:
                    raise ValueError(f"Unrecognized key: {key}: {value}")


@pytest.mark.parametrize("path", ["hydrogen_optimization.yaml"])
def test_hydrogen(path):
    """Test the hydrogen optimization file at path."""
    yaml_path = get_yaml_path(path)
    with open(yaml_path, "rt") as yaml_file:
        yaml_data = yaml.load(yaml_file, Loader=UniqueKeyLoader)
        for residue in yaml_data:
            for key, value in residue.items():
                if key == "name":
                    assert isinstance(key, str)
                elif key == "type":
                    if value not in OPT_TYPES:
                        raise ValueError(f"Unknown optimization type: {value}")
                elif key == "angle":
                    assert isinstance(value, list)
                    assert len(value) == 4
                elif key == "atoms":
                    for atom in value:
                        for atom_key, atom_value in atom.items():
                            if atom_key in ["name", "bond"]:
                                assert isinstance(atom_value, str)
                            else:
                                raise ValueError(
                                    f"Unknown atom "
                                    f"attribute: {atom_key}"
                                )
                else:
                    raise ValueError(f"Unknown entry: {key}: {value}")


@pytest.mark.parametrize(
    "name_path",
    [
        "amber_names", "charmm_names", "parse_names", "peoepb_names",
        "swanson_names", "tyl06_names"
    ]
)
def test_name(name_path):
    """Test parsing of renaming files."""

    yaml_path = get_yaml_path(name_path)
    print(f"Reading data from {yaml_path}")
    with open(yaml_path, "rt") as yaml_file:
        parse_name(yaml_file)


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
        "peoepb_parameters", "swanson_parameters", "tyl06_parameters"
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


def test_zzz():
    """This is a bogus test designed to fail."""

    raise Exception(
        "This test suite is incomplete! It needs to include other YAML files "
        "as well as update the code to use YAML instead of XML."
    )
