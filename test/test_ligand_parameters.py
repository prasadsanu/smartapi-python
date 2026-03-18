"""
Tests for SmartApi.ligand_parameters
=====================================
All tests use temporary in-memory files and do not require Gaussian,
AmberTools, or LigParGen to be installed.
"""

import os
import sys
import tempfile
import textwrap
import unittest
import warnings

root_directory = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir)
)
sys.path.insert(0, root_directory)

from SmartApi.ligand_parameters import (
    KCAL_MOL_TO_K,
    KJ_MOL_TO_K,
    ANGSTROM_PER_NM,
    ChargeReadError,
    GeometryExtractionError,
    LigandParameterError,
    LJParameterError,
    generate_parameters,
    parse_gaussian_gjf,
    parse_gaussian_log,
    parse_geometry,
    read_ligpargen_parameters,
    read_resp_charges,
    validate_parameters,
    write_parameter_file,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_tmp(content: str, suffix: str) -> str:
    """Write *content* to a temp file with *suffix* and return the path."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "w") as fh:
        fh.write(textwrap.dedent(content))
    return path


# Minimal but realistic Gaussian log snippet (methanol, 6 atoms)
_GAUSSIAN_LOG = """\
 Gaussian 16:  AMD64-G16RevC.01 3-Jul-2019
 #P HF/6-31G* Opt

 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        0.0000      0.0000      0.0000
      2          8           0        1.4300      0.0000      0.0000
      3          1           0       -0.3700      1.0280      0.0000
      4          1           0       -0.3700     -0.5140     -0.8900
      5          1           0       -0.3700     -0.5140      0.8900
      6          1           0        1.7200      0.9500      0.0000
 ---------------------------------------------------------------------
 Distance matrix (angstroms):
"""

# Second Standard-orientation block (geometry step 2) – should be returned
_GAUSSIAN_LOG_TWO_STEPS = _GAUSSIAN_LOG + """\
 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        0.0010      0.0010      0.0010
      2          8           0        1.4310      0.0010      0.0010
      3          1           0       -0.3690      1.0290      0.0010
      4          1           0       -0.3690     -0.5130     -0.8890
      5          1           0       -0.3690     -0.5130      0.8910
      6          1           0        1.7210      0.9510      0.0010
 ---------------------------------------------------------------------
 Distance matrix (angstroms):
"""

# A minimal Gaussian .gjf for the same molecule
_GAUSSIAN_GJF = """\
 %mem=4GB
 #P HF/6-31G* Opt ESP

 Methanol

 0 1
 C   0.0000   0.0000   0.0000
 O   1.4300   0.0000   0.0000
 H  -0.3700   1.0280   0.0000
 H  -0.3700  -0.5140  -0.8900
 H  -0.3700  -0.5140   0.8900
 H   1.7200   0.9500   0.0000

"""

# Minimal antechamber MOL2 charge section (6 atoms)
_MOL2_CHARGES = """\
@<TRIPOS>MOLECULE
methanol
 6 5 0 0 0
SMALL
RESP Charges

@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 c3      1  LIG      -0.269700
      2 O1          1.4300    0.0000    0.0000 oh      1  LIG      -0.654900
      3 H1         -0.3700    1.0280    0.0000 hc      1  LIG       0.085100
      4 H2         -0.3700   -0.5140   -0.8900 hc      1  LIG       0.085100
      5 H3         -0.3700   -0.5140    0.8900 hc      1  LIG       0.085100
      6 H4          1.7200    0.9500    0.0000 ho      1  LIG       0.404400
@<TRIPOS>BOND
"""

# Plain-text charge file (one float per line)
_PLAIN_CHARGES = """\
-0.2697
-0.6549
 0.0851
 0.0851
 0.0851
 0.4044
"""

# GROMACS .itp atomtypes section (LigParGen output)
_ITP_PARAMS = """\
[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1 3 yes 0.5 0.5

[ atomtypes ]
; name  at.num  mass  charge  ptype  sigma(nm)  epsilon(kJ/mol)
 opls_135   6   12.011   0.000   A   0.35000E-00   0.27614E+00
 opls_154   8   15.999   0.000   A   0.30650E-00   0.71132E+00
 opls_140   1    1.008   0.000   A   0.25000E-00   0.12552E+00

[ moleculetype ]
"""

# OPLS .prm file (sigma in Å, epsilon in kcal/mol)
_PRM_PARAMS = """\
! OPLS-AA parameters
! type  sigma(Ang)  epsilon(kcal/mol)
opls_135  3.5000  0.0660
opls_154  3.0650  0.1700
opls_140  2.5000  0.0300
"""


# ---------------------------------------------------------------------------
# Unit-conversion constant tests
# ---------------------------------------------------------------------------

class TestUnitConstants(unittest.TestCase):
    def test_kcal_mol_to_K(self):
        self.assertAlmostEqual(KCAL_MOL_TO_K, 503.22, places=2)

    def test_kJ_mol_to_K(self):
        self.assertAlmostEqual(KJ_MOL_TO_K, 120.272, places=3)

    def test_angstrom_per_nm(self):
        self.assertAlmostEqual(ANGSTROM_PER_NM, 10.0, places=5)


# ---------------------------------------------------------------------------
# Geometry parsing – Gaussian log
# ---------------------------------------------------------------------------

class TestParseGaussianLog(unittest.TestCase):

    def setUp(self):
        self.log_path = _write_tmp(_GAUSSIAN_LOG, ".log")
        self.two_step_path = _write_tmp(_GAUSSIAN_LOG_TWO_STEPS, ".log")

    def tearDown(self):
        for p in (self.log_path, self.two_step_path):
            if os.path.isfile(p):
                os.remove(p)

    def test_returns_six_atoms(self):
        atoms = parse_gaussian_log(self.log_path)
        self.assertEqual(len(atoms), 6)

    def test_first_atom_is_carbon(self):
        atoms = parse_gaussian_log(self.log_path)
        self.assertEqual(atoms[0]["element"], "C")

    def test_second_atom_is_oxygen(self):
        atoms = parse_gaussian_log(self.log_path)
        self.assertEqual(atoms[1]["element"], "O")

    def test_coordinates_correct(self):
        atoms = parse_gaussian_log(self.log_path)
        self.assertAlmostEqual(atoms[0]["x"], 0.0)
        self.assertAlmostEqual(atoms[1]["x"], 1.43)

    def test_atom_ids_sequential(self):
        atoms = parse_gaussian_log(self.log_path)
        ids = [a["atom_id"] for a in atoms]
        self.assertEqual(ids, list(range(1, 7)))

    def test_last_block_returned_when_two_steps(self):
        """The second (last) Standard-orientation block must be returned."""
        atoms = parse_gaussian_log(self.two_step_path)
        # The second block shifts every coordinate by +0.001
        self.assertAlmostEqual(atoms[0]["x"], 0.001, places=4)

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            parse_gaussian_log("/nonexistent/path.log")

    def test_empty_file_raises(self):
        path = _write_tmp("", ".log")
        try:
            with self.assertRaises(GeometryExtractionError):
                parse_gaussian_log(path)
        finally:
            os.remove(path)


# ---------------------------------------------------------------------------
# Geometry parsing – Gaussian gjf
# ---------------------------------------------------------------------------

class TestParseGaussianGjf(unittest.TestCase):

    def setUp(self):
        self.gjf_path = _write_tmp(_GAUSSIAN_GJF, ".gjf")

    def tearDown(self):
        if os.path.isfile(self.gjf_path):
            os.remove(self.gjf_path)

    def test_returns_six_atoms(self):
        atoms = parse_gaussian_gjf(self.gjf_path)
        self.assertEqual(len(atoms), 6)

    def test_first_atom_element(self):
        atoms = parse_gaussian_gjf(self.gjf_path)
        self.assertEqual(atoms[0]["element"], "C")

    def test_coordinates(self):
        atoms = parse_gaussian_gjf(self.gjf_path)
        self.assertAlmostEqual(atoms[1]["x"], 1.43)

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            parse_gaussian_gjf("/nonexistent/path.gjf")


# ---------------------------------------------------------------------------
# parse_geometry dispatch
# ---------------------------------------------------------------------------

class TestParseGeometryDispatch(unittest.TestCase):

    def setUp(self):
        self.log_path = _write_tmp(_GAUSSIAN_LOG, ".log")
        self.gjf_path = _write_tmp(_GAUSSIAN_GJF, ".gjf")

    def tearDown(self):
        for p in (self.log_path, self.gjf_path):
            if os.path.isfile(p):
                os.remove(p)

    def test_log_dispatch(self):
        atoms = parse_geometry(self.log_path)
        self.assertEqual(len(atoms), 6)

    def test_gjf_dispatch(self):
        atoms = parse_geometry(self.gjf_path)
        self.assertEqual(len(atoms), 6)


# ---------------------------------------------------------------------------
# RESP charge reader
# ---------------------------------------------------------------------------

class TestReadRespCharges(unittest.TestCase):

    def setUp(self):
        self.mol2_path = _write_tmp(_MOL2_CHARGES, ".mol2")
        self.plain_path = _write_tmp(_PLAIN_CHARGES, ".txt")

    def tearDown(self):
        for p in (self.mol2_path, self.plain_path):
            if os.path.isfile(p):
                os.remove(p)

    def test_mol2_returns_six_charges(self):
        charges = read_resp_charges(self.mol2_path)
        self.assertEqual(len(charges), 6)

    def test_mol2_first_charge(self):
        charges = read_resp_charges(self.mol2_path)
        self.assertAlmostEqual(charges[0], -0.2697, places=4)

    def test_mol2_last_charge(self):
        charges = read_resp_charges(self.mol2_path)
        self.assertAlmostEqual(charges[5], 0.4044, places=4)

    def test_mol2_charge_sum(self):
        charges = read_resp_charges(self.mol2_path)
        self.assertAlmostEqual(sum(charges), -0.2697 - 0.6549 + 3 * 0.0851
                               + 0.4044, places=4)

    def test_plain_returns_six_charges(self):
        charges = read_resp_charges(self.plain_path)
        self.assertEqual(len(charges), 6)

    def test_plain_first_charge(self):
        charges = read_resp_charges(self.plain_path)
        self.assertAlmostEqual(charges[0], -0.2697, places=4)

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            read_resp_charges("/nonexistent/charges.mol2")

    def test_empty_file_raises(self):
        path = _write_tmp("", ".txt")
        try:
            with self.assertRaises(ChargeReadError):
                read_resp_charges(path)
        finally:
            os.remove(path)


# ---------------------------------------------------------------------------
# LigParGen / LJ parameter reader
# ---------------------------------------------------------------------------

class TestReadLigpargenParameters(unittest.TestCase):

    def setUp(self):
        self.itp_path = _write_tmp(_ITP_PARAMS, ".itp")
        self.prm_path = _write_tmp(_PRM_PARAMS, ".prm")

    def tearDown(self):
        for p in (self.itp_path, self.prm_path):
            if os.path.isfile(p):
                os.remove(p)

    # -- ITP --

    def test_itp_returns_three_types(self):
        params = read_ligpargen_parameters(self.itp_path)
        self.assertEqual(len(params), 3)

    def test_itp_has_carbon_type(self):
        params = read_ligpargen_parameters(self.itp_path)
        self.assertIn("opls_135", params)

    def test_itp_sigma_conversion_nm_to_ang(self):
        """sigma_nm * 10 should equal sigma_ang."""
        params = read_ligpargen_parameters(self.itp_path)
        # opls_135 sigma = 0.35000 nm → 3.5000 Å
        self.assertAlmostEqual(params["opls_135"]["sigma_ang"], 3.5, places=4)

    def test_itp_epsilon_conversion_kj_to_K(self):
        """epsilon_kJ * 120.272 should equal epsilon_K."""
        params = read_ligpargen_parameters(self.itp_path)
        expected = 0.27614 * KJ_MOL_TO_K
        self.assertAlmostEqual(
            params["opls_135"]["epsilon_K"], expected, places=2
        )

    # -- PRM --

    def test_prm_returns_three_types(self):
        params = read_ligpargen_parameters(self.prm_path)
        self.assertEqual(len(params), 3)

    def test_prm_epsilon_conversion_kcal_to_K(self):
        """epsilon_kcal * 503.22 should equal epsilon_K."""
        params = read_ligpargen_parameters(self.prm_path)
        expected = 0.0660 * KCAL_MOL_TO_K
        self.assertAlmostEqual(
            params["opls_135"]["epsilon_K"], expected, places=2
        )

    def test_prm_sigma_unchanged(self):
        params = read_ligpargen_parameters(self.prm_path)
        self.assertAlmostEqual(params["opls_135"]["sigma_ang"], 3.5, places=4)

    # -- Error cases --

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            read_ligpargen_parameters("/nonexistent/params.itp")

    def test_unsupported_extension_raises(self):
        path = _write_tmp("dummy", ".xyz")
        try:
            with self.assertRaises(LJParameterError):
                read_ligpargen_parameters(path)
        finally:
            os.remove(path)

    def test_empty_itp_raises(self):
        path = _write_tmp("", ".itp")
        try:
            with self.assertRaises(LJParameterError):
                read_ligpargen_parameters(path)
        finally:
            os.remove(path)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class TestValidateParameters(unittest.TestCase):

    def _make_atoms(self, n=6):
        return [
            {"atom_id": i + 1, "element": "C", "x": 0.0, "y": 0.0, "z": 0.0}
            for i in range(n)
        ]

    def test_valid_passes_silently(self):
        # Near-zero charge sum → should not raise or warn
        charges = [0.0] * 6
        validate_parameters(self._make_atoms(6), charges)

    def test_count_mismatch_raises(self):
        with self.assertRaises(LigandParameterError):
            validate_parameters(self._make_atoms(6), [0.0] * 5)

    def test_atom_type_mismatch_raises(self):
        with self.assertRaises(LigandParameterError):
            validate_parameters(
                self._make_atoms(6), [0.0] * 6, atom_types=["C"] * 5
            )

    def test_large_charge_sum_warns(self):
        charges = [0.5] * 6  # sum = 3.0, deviation from 3 = 0.0 → OK
        # Actually sum = 3.0, nearest int = 3, deviation = 0.0 → no warning
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            validate_parameters(self._make_atoms(6), charges)
        # Deviation from int is 0 → no UserWarning
        user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
        self.assertEqual(len(user_warnings), 0)

    def test_non_integer_charge_sum_warns(self):
        charges = [0.1] * 6  # sum = 0.6, deviation from 1 = 0.4 > 0.01
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            validate_parameters(self._make_atoms(6), charges)
        user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
        self.assertGreater(len(user_warnings), 0)


# ---------------------------------------------------------------------------
# write_parameter_file
# ---------------------------------------------------------------------------

class TestWriteParameterFile(unittest.TestCase):

    def setUp(self):
        self.log_path = _write_tmp(_GAUSSIAN_LOG, ".log")
        self.mol2_path = _write_tmp(_MOL2_CHARGES, ".mol2")
        self.itp_path = _write_tmp(_ITP_PARAMS, ".itp")
        self.out_dir = tempfile.mkdtemp()

    def tearDown(self):
        for p in (self.log_path, self.mol2_path, self.itp_path):
            if os.path.isfile(p):
                os.remove(p)

    def _atoms(self):
        return parse_gaussian_log(self.log_path)

    def _charges(self):
        return read_resp_charges(self.mol2_path)

    def test_file_created(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        write_parameter_file(self._atoms(), self._charges(), output_path=out)
        self.assertTrue(os.path.isfile(out))

    def test_header_present(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        write_parameter_file(self._atoms(), self._charges(), output_path=out)
        with open(out) as fh:
            first_line = fh.readline()
        self.assertIn("AtomID", first_line)
        self.assertIn("Charge", first_line)
        self.assertIn("Sigma", first_line)

    def test_correct_number_of_data_lines(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        atoms = self._atoms()
        write_parameter_file(atoms, self._charges(), output_path=out)
        with open(out) as fh:
            lines = [line for line in fh if line.strip() and not line.startswith("-")]
        # header line + 6 data lines = 7 non-separator lines
        self.assertEqual(len(lines), 7)

    def test_lj_params_written(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        lj = read_ligpargen_parameters(self.itp_path)
        atoms = self._atoms()
        types = ["opls_135", "opls_154", "opls_140",
                 "opls_140", "opls_140", "opls_140"]
        write_parameter_file(
            atoms, self._charges(), lj_params=lj,
            atom_types=types, output_path=out
        )
        with open(out) as fh:
            content = fh.read()
        # sigma for opls_135 = 3.5000 Å
        self.assertIn("3.5000", content)

    def test_count_mismatch_raises(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        with self.assertRaises(LigandParameterError):
            write_parameter_file(self._atoms(), [0.0] * 5, output_path=out)

    def test_no_overlap_in_columns(self):
        """All column values must be separated by at least one space."""
        out = os.path.join(self.out_dir, "parameter.txt")
        write_parameter_file(self._atoms(), self._charges(), output_path=out)
        with open(out) as fh:
            lines = fh.readlines()
        # Skip header and separator lines
        for line in lines[2:]:
            parts = line.split()
            # Should have 6 data columns (no LJ → sigma/epsilon = N/A)
            self.assertGreaterEqual(len(parts), 6)


# ---------------------------------------------------------------------------
# generate_parameters (full workflow)
# ---------------------------------------------------------------------------

class TestGenerateParameters(unittest.TestCase):

    def setUp(self):
        self.log_path = _write_tmp(_GAUSSIAN_LOG, ".log")
        self.gjf_path = _write_tmp(_GAUSSIAN_GJF, ".gjf")
        self.mol2_path = _write_tmp(_MOL2_CHARGES, ".mol2")
        self.plain_path = _write_tmp(_PLAIN_CHARGES, ".txt")
        self.itp_path = _write_tmp(_ITP_PARAMS, ".itp")
        self.out_dir = tempfile.mkdtemp()

    def tearDown(self):
        for p in (
            self.log_path, self.gjf_path, self.mol2_path,
            self.plain_path, self.itp_path,
        ):
            if os.path.isfile(p):
                os.remove(p)

    def test_returns_absolute_path(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        result = generate_parameters(self.log_path, self.mol2_path,
                                     output_path=out)
        self.assertTrue(os.path.isabs(result))

    def test_output_file_exists(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        generate_parameters(self.log_path, self.mol2_path, output_path=out)
        self.assertTrue(os.path.isfile(out))

    def test_with_gjf_input(self):
        out = os.path.join(self.out_dir, "parameter_gjf.txt")
        generate_parameters(self.gjf_path, self.mol2_path, output_path=out)
        self.assertTrue(os.path.isfile(out))

    def test_with_ligpargen_itp(self):
        out = os.path.join(self.out_dir, "parameter_lj.txt")
        types = ["opls_135", "opls_154", "opls_140",
                 "opls_140", "opls_140", "opls_140"]
        generate_parameters(
            self.log_path, self.mol2_path,
            ligpargen_file=self.itp_path,
            atom_types=types,
            output_path=out,
        )
        with open(out) as fh:
            content = fh.read()
        self.assertIn("Sigma", content)
        self.assertNotIn("N/A", content)  # all types should be resolved

    def test_with_plain_charges(self):
        out = os.path.join(self.out_dir, "parameter_plain.txt")
        generate_parameters(self.log_path, self.plain_path, output_path=out)
        self.assertTrue(os.path.isfile(out))

    def test_missing_gaussian_file_raises(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        with self.assertRaises(FileNotFoundError):
            generate_parameters("/no/such/file.log", self.mol2_path,
                                output_path=out)

    def test_missing_charge_file_raises(self):
        out = os.path.join(self.out_dir, "parameter.txt")
        with self.assertRaises(FileNotFoundError):
            generate_parameters(self.log_path, "/no/charges.mol2",
                                output_path=out)

    def test_atom_charge_count_mismatch_raises(self):
        """Charge file with wrong count must raise LigandParameterError."""
        bad_charges = _write_tmp("-0.1\n-0.2\n-0.3\n", ".txt")
        out = os.path.join(self.out_dir, "parameter.txt")
        try:
            with self.assertRaises(LigandParameterError):
                generate_parameters(self.log_path, bad_charges, output_path=out)
        finally:
            os.remove(bad_charges)


if __name__ == "__main__":
    unittest.main(verbosity=2)
