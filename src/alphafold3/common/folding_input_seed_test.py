import sys
from unittest.mock import MagicMock

# Mock necessary modules to avoid import errors
sys.modules['alphafold3.cpp'] = MagicMock()
sys.modules['alphafold3.cpp.cif_dict'] = MagicMock()
sys.modules['alphafold3.cpp.cif_dict'].parse_multi_data_cif.return_value = {}
sys.modules['rdkit'] = MagicMock()
sys.modules['rdkit.Chem'] = MagicMock()
sys.modules['zstandard'] = MagicMock()
sys.modules['jax'] = MagicMock()
sys.modules['jax.numpy'] = MagicMock()
sys.modules['tokamax'] = MagicMock()

import unittest
from alphafold3.common import folding_input

class FoldingInputSeedTest(unittest.TestCase):

  def test_with_multiple_seeds_single_seed(self):
    """Tests with_multiple_seeds when input has a single seed."""
    test_input = folding_input.Input(
        name='test',
        chains=[],
        rng_seeds=[1]
    )
    updated_input = test_input.with_multiple_seeds(5)
    self.assertEqual(list(updated_input.rng_seeds), [1, 2, 3, 4, 5])

  def test_with_multiple_seeds_multiple_seeds_override(self):
    """Tests with_multiple_seeds when input already has multiple seeds (Issue #622)."""
    test_input = folding_input.Input(
        name='test',
        chains=[],
        rng_seeds=[10, 20, 30]
    )
    # This used to raise ValueError, now it should work and use the first seed as base.
    updated_input = test_input.with_multiple_seeds(3)
    self.assertEqual(list(updated_input.rng_seeds), [10, 11, 12])

  def test_with_multiple_seeds_invalid_count(self):
    """Tests that with_multiple_seeds raises ValueError for invalid seed counts."""
    # We can't easily instantiate Input with full validation without more mocks,
    # but we can try basic initialization for logic testing.
    # Note: Input.__post_init__ might fail if we don't mock more things.
    # Let's use a mock for Input if needed, but we want to test the real method.
    
    test_input = folding_input.Input(
        name='test',
        chains=[],
        rng_seeds=[123]
    )
    with self.assertRaisesRegex(ValueError, 'Number of seeds must be greater than 1'):
      test_input.with_multiple_seeds(1)
    with self.assertRaisesRegex(ValueError, 'Number of seeds must be greater than 1'):
      test_input.with_multiple_seeds(0)

if __name__ == '__main__':
  unittest.main()
