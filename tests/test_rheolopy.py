import pytest
import numpy as np
from rheolopy.materials import Material
from rheolopy.core import (
    sigma_byerlee,
    eta_dislocation,
    eta_diffusion,
    eta_effective,
    calc_peierls,
    sigma_d,
    compute_dsigma
)

# Use precise gravity as updated in core.py
R = 8.314472
g = 9.80665

def test_sigma_byerlee_compression():
    mat = Material(
        id="test_brittle",
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36
    )
    z = -10000.0
    expected = 3.0 * 2850 * g * abs(z) * (1.0 - 0.36)
    result = sigma_byerlee(mat, z, mode="compression")
    assert np.isclose(result, expected)

def test_sigma_byerlee_extension():
    mat = Material(
        id="test_brittle",
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36
    )
    z = -10000.0
    expected = 0.75 * 2850 * g * abs(z) * (1.0 - 0.36)
    result = sigma_byerlee(mat, z, mode="extension")
    assert np.isclose(result, expected)

def test_sigma_byerlee_invalid_mode():
    mat = Material(id="test_brittle", rho_b=2850, fc_c=3.0, fc_e=0.75, lambda_pore=0.36)
    with pytest.raises(ValueError):
        sigma_byerlee(mat, -10000.0, mode="invalid")

def test_eta_dislocation():
    mat = Material(
        id="test_disloc",
        a_disloc=6.31e-20,
        n=3.0,
        q_disloc=356000,
        d=1e-3,
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36
    )
    eta = eta_dislocation(mat, temp=800.0, strain_rate=1e-15)
    assert np.isfinite(eta)
    assert eta > 0
    assert 1e18 <= eta <= 1e26

def test_eta_dislocation_none_returns_nan():
    mat = Material(id="test_none", a_disloc=None)
    eta = eta_dislocation(mat, temp=800.0, strain_rate=1e-15)
    assert np.isnan(eta)

def test_eta_diffusion():
    mat = Material(
        id="test_diff",
        a_diff=1.5e-15,
        q_diff=300000,
        d=1e-3,
        m=2.5
    )
    eta = eta_diffusion(mat, temp=900.0, strain_rate=1e-15)
    assert np.isfinite(eta)
    assert eta > 0

def test_eta_diffusion_none_returns_nan():
    mat = Material(id="test_none", a_diff=None)
    eta = eta_diffusion(mat, temp=900.0, strain_rate=1e-15)
    assert np.isnan(eta)

def test_eta_effective_harmonic_mean():
    mat = Material(
        id="test_both",
        a_disloc=6.31e-20,
        n=3.0,
        q_disloc=356000,
        a_diff=1.5e-15,
        q_diff=300000,
        d=1e-3,
        m=2.5
    )
    eta_eff, eta_dis, eta_diff = eta_effective(mat, temp=900.0, strain_rate=1e-15)
    assert np.isfinite(eta_eff)
    assert eta_eff <= min(eta_dis, eta_diff)
    expected_eff = 1.0 / (1.0 / eta_dis + 1.0 / eta_diff)
    assert np.isclose(eta_eff, expected_eff)

def test_eta_effective_single_mechanism():
    mat = Material(
        id="test_single",
        a_disloc=6.31e-20,
        n=3.0,
        q_disloc=356000,
        a_diff=None
    )
    eta_eff, eta_dis, eta_diff = eta_effective(mat, temp=900.0, strain_rate=1e-15)
    assert np.isnan(eta_diff)
    assert np.isclose(eta_eff, eta_dis)

def test_calc_peierls_returns_zero_without_params():
    mat = Material(id="test_no_peierls")
    peierls = calc_peierls(mat, temp=800.0, strain_rate=1e-15)
    assert peierls == 0.0

def test_calc_peierls_basic():
    mat = Material(
        id="test_peierls",
        a_disloc=1.1e-16,
        n=3.5,
        q_disloc=530000,
        eps_peierls=5.7e11,
        sigma_peierls=8.5e9,
        stress_pd=200e6
    )
    # At very low temperature (e.g. 600 K), Peierls stress should be active and large
    peierls = calc_peierls(mat, temp=600.0, strain_rate=1e-15)
    assert peierls > 0
    assert np.isfinite(peierls)

def test_sigma_d_with_peierls():
    mat = Material(
        id="test_envelope",
        rho_b=3300,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.0,
        a_disloc=1.1e-16,
        n=3.5,
        q_disloc=530000,
        a_diff=None,
        eps_peierls=5.7e11,
        sigma_peierls=8.5e9,
        stress_pd=200e6
    )
    
    # Very cold, deep: creep is huge, Peierls kicks in
    sigma_peierls_val = calc_peierls(mat, temp=500.0, strain_rate=1e-15)
    eta_eff, _, _ = eta_effective(mat, temp=500.0, strain_rate=1e-15)
    creep = 2 * eta_eff * 1e-15
    
    # Assert that standard creep is indeed above Peierls threshold
    assert creep > 200e6
    assert sigma_peierls_val > 0 and sigma_peierls_val < creep
    
    # Compute overall stress
    is_peierls = sigma_d(mat, z=-40000.0, temp=500.0, strain_rate=1e-15, mode="compression", return_index=True)
    sigma_total = sigma_d(mat, z=-40000.0, temp=500.0, strain_rate=1e-15, mode="compression")
    
    # In some cold settings, Byerlee might still be smaller than Peierls.
    s_byerlee = sigma_byerlee(mat, -40000.0, mode="compression")
    if s_byerlee > sigma_peierls_val:
        assert is_peierls == 1
        assert np.isclose(sigma_total, sigma_peierls_val)
    else:
        assert is_peierls == 0

def test_sigma_d_selects_minimum():
    mat = Material(
        id="test_envelope",
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36,
        a_disloc=6.31e-20,
        n=3.0,
        q_disloc=356000,
        a_diff=1.5e-15,
        q_diff=300000,
        d=1e-3,
        m=2.5
    )
    
    # Shallow, cold: brittle should dominate (smaller than creep)
    sigma_shallow = sigma_d(mat, z=-1000.0, temp=300.0, strain_rate=1e-15, mode="compression")
    sigma_shallow_byerlee = sigma_byerlee(mat, -1000.0, mode="compression")
    assert np.isclose(sigma_shallow, sigma_shallow_byerlee)

    # Deep, hot: creep should dominate (smaller than brittle)
    sigma_deep = sigma_d(mat, z=-50000.0, temp=1200.0, strain_rate=1e-15, mode="compression")
    sigma_deep_byerlee = sigma_byerlee(mat, -50000.0, mode="compression")
    assert sigma_deep < sigma_deep_byerlee

def test_sigma_d_positive_depth_raises():
    mat = Material(id="test_pos")
    with pytest.raises(ValueError):
        sigma_d(mat, z=1000.0, temp=300.0)

def test_sigma_d_zero_depth():
    mat = Material(
        id="test_zero",
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36,
        a_disloc=6.31e-20,
        n=3.0,
        q_disloc=356000
    )
    sigma_surface = sigma_d(mat, z=0.0, temp=300.0, mode="compression")
    assert np.isclose(sigma_surface, 0.0)

def test_compute_dsigma_scalar():
    mat = Material(
        id="test_scalar",
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36
    )
    dsigma_c, dsigma_e = compute_dsigma(mat, z=-1000.0, T=300.0, strain_rate=1e-15)
    assert isinstance(dsigma_c, float)
    assert isinstance(dsigma_e, float)
    assert dsigma_c < 0
    assert dsigma_e > 0

def test_compute_dsigma_array():
    mat = Material(
        id="test_array",
        rho_b=2850,
        fc_c=3.0,
        fc_e=0.75,
        lambda_pore=0.36,
        a_disloc=6.31e-20,
        n=3.0,
        q_disloc=356000,
        d=1e-3
    )
    z = np.array([-1000.0, -2000.0, -3000.0])
    T = np.array([300.0, 400.0, 500.0])
    
    dsigma, depths = compute_dsigma(mat, z, T, strain_rate=1e-15)
    
    assert isinstance(dsigma, np.ndarray)
    assert isinstance(depths, np.ndarray)
    assert len(dsigma) == 2 * len(z)
    assert len(depths) == 2 * len(z)
    # The first half should be compression (negative)
    assert np.all(dsigma[:len(z)] < 0)
    # The second half should be extension (positive)
    assert np.all(dsigma[len(z):] > 0)
