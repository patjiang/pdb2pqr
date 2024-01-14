"""Test behavior of YAML-format input files."""
import yaml
import logging
from pdb2pqr.io import test_yaml_file as get_yaml_path


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


def test_definitions():
    """Test parsing of definition files."""

    for def_path in ["aa_definitions", "na_definitions"]:
        yaml_path = get_yaml_path(def_path)
        print(f"Reading data from {yaml_path}")
        with open(yaml_path, "rt") as yaml_file:
            parse_definition(yaml_file)


def last_test():
    """This is a bogus test designed to fail."""

    raise Exception(
        "This test suite is incomplete! It needs to include other YAML files as well as "
        "update the code to use YAML instead of XML."
    )
