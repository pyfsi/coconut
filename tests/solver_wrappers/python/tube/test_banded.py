import coconut.coupling_components.solver_wrappers.python.banded as bnd

import numpy as np
import unittest


class TestBanded(unittest.TestCase):

    def setUp(self):
        a = np.random.rand(20) * 10
        self.dense_matrix = np.array([
            [a[0], a[1], 0,     0,     0,     0],
            [a[2], a[3], a[4],  0,     0,     0],
            [a[5], a[6], a[7],  a[8],  0,     0],
            [0,    a[9], a[10], a[11], a[12], 0],
            [0,    0,    a[13], a[14], a[15], a[16]],
            [0,    0,    0,     a[17], a[18], a[19]]
        ])
        self.banded_matrix = np.array([
            [0,    a[1], a[4],  a[8],  a[12], a[16]],
            [a[0], a[3], a[7],  a[11], a[15], a[19]],
            [a[2], a[6], a[10], a[14], a[18], 0],
            [a[5], a[9], a[13], a[17], 0,     0]
        ])
        self.au = 1
        self.al = 2

    def test_check_banded(self):
        banded, au, al = bnd.check_banded(self.dense_matrix)
        self.assertFalse(banded)
        banded, au, al = bnd.check_banded(self.banded_matrix)
        self.assertTrue(banded)
        self.assertEqual((au, al), (self.au, self.al))

    def test_to_dense(self):
        d = bnd.to_dense(self.banded_matrix)
        np.testing.assert_array_equal(d, self.dense_matrix)
        d = bnd.to_dense(bnd.transpose_banded(self.banded_matrix))
        np.testing.assert_array_equal(d, self.dense_matrix.T)

    def test_to_banded(self):
        with self.assertRaises(ValueError):
            bnd.to_banded(self.banded_matrix[:, 1:], self.au)
        b = bnd.to_banded(self.dense_matrix, self.au)
        np.testing.assert_array_equal(b, self.banded_matrix[0:2 * self.au + 1])
        b = bnd.to_banded(self.dense_matrix, self.au, self.al)
        np.testing.assert_array_equal(b, self.banded_matrix)
        b = bnd.to_banded(self.dense_matrix.T, self.al, self.au)
        np.testing.assert_array_equal(b, bnd.transpose_banded(self.banded_matrix))

    def test_remove_boundaries(self):
        with self.assertRaises(ValueError):
            bnd.remove_boundaries(self.dense_matrix, 1)
            bnd.remove_boundaries(self.banded_matrix, -1)
        np.testing.assert_array_equal(self.banded_matrix, bnd.remove_boundaries(self.banded_matrix, 0))
        for num in (1, 2):
            b = bnd.remove_boundaries(self.banded_matrix.copy(), num)
            d = self.dense_matrix[num:-num, num:-num]
            np.testing.assert_array_equal(bnd.to_dense(b), d)

    def test_multiply_banded(self):
        with self.assertRaises(ValueError):
            bnd.multiply_banded(self.banded_matrix, self.dense_matrix)
            bnd.multiply_banded(self.dense_matrix, self.banded_matrix)
            bnd.multiply_banded(self.banded_matrix[:-1, :], self.banded_matrix)
        d2 = self.dense_matrix @ self.dense_matrix
        b2 = bnd.multiply_banded(self.banded_matrix, self.banded_matrix)
        np.testing.assert_array_almost_equal(bnd.to_dense(b2), d2)
        b2 = bnd.multiply_banded(self.banded_matrix, self.banded_matrix, output_dense=True)
        np.testing.assert_array_almost_equal(b2, d2)
        bt = bnd.transpose_banded(self.banded_matrix)
        dtdt = self.dense_matrix.T @ self.dense_matrix.T
        btbt = bnd.multiply_banded(bt, bt, output_dense=True)
        np.testing.assert_array_almost_equal(btbt, dtdt)
        ddt = self.dense_matrix @ self.dense_matrix.T
        bbt = bnd.multiply_banded(self.banded_matrix, bt, output_dense=True)
        np.testing.assert_array_almost_equal(bbt, ddt)
        dtd = self.dense_matrix.T @ self.dense_matrix
        btb = bnd.multiply_banded(bt, self.banded_matrix, output_dense=True)
        np.testing.assert_array_almost_equal(btb, dtd)

    def test_multiply_banded_vector(self):
        a = np.random.rand(6) * 10
        self.vector = np.array([a[0], a[1],  a[2],  a[3], a[4], a[5]])
        with self.assertRaises(ValueError):
            bnd.multiply_banded_vector(self.dense_matrix, self.vector)
            bnd.multiply_banded_vector(self.banded_matrix, self.banded_matrix)
            bnd.multiply_banded_vector(self.banded_matrix, 1)
            bnd.multiply_banded_vector(self.banded_matrix, np.array([a[0], a[1]]))
            bnd.multiply_banded_vector(self.banded_matrix[:-1, :], self.vector)
        # with numpy
        b = (self.dense_matrix @ self.vector.reshape(-1, 1)).flatten()
        # 1D vector
        b1 = bnd.multiply_banded_vector(self.banded_matrix, self.vector)
        self.assertTrue(b1.shape == b.shape)
        np.testing.assert_array_almost_equal(b, b1)
        # 2D vector
        b2 = bnd.multiply_banded_vector(self.banded_matrix, self.vector.reshape(-1, 1))
        self.assertTrue(b2.shape == (b.size, 1))
        np.testing.assert_array_almost_equal(b, b2.flatten())

    def test_transpose_banded(self):
        np.testing.assert_array_equal(self.dense_matrix.T, bnd.to_dense(bnd.transpose_banded(self.banded_matrix)))


if __name__ == '__main__':
    unittest.main()
