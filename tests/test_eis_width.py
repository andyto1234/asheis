import numpy as np
import pytest

from asheis.eis_width.calculation import _instrument_width_array


def test_instrument_width_array_accepts_row_values():
    result = _instrument_width_array([1.0, 2.0, 3.0], (3, 4))

    expected = np.array(
        [
            [1.0, 1.0, 1.0, 1.0],
            [2.0, 2.0, 2.0, 2.0],
            [3.0, 3.0, 3.0, 3.0],
        ]
    )
    np.testing.assert_allclose(result, expected)


def test_instrument_width_array_accepts_column_values():
    result = _instrument_width_array([1.0, 2.0, 3.0, 4.0], (3, 4))

    expected = np.array(
        [
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
            [1.0, 2.0, 3.0, 4.0],
        ]
    )
    np.testing.assert_allclose(result, expected)


def test_instrument_width_array_rejects_incompatible_shape():
    with pytest.raises(ValueError, match="Cannot broadcast slit_width"):
        _instrument_width_array([1.0, 2.0], (3, 4))

