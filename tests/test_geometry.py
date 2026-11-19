"""Tests for the grain geometry system."""
import numpy as np
import pytest
from srm_1d.grain_geometry import (
    build_snapped_geometry, update_cell_geometry, advance_bore_regression,
    advance_endface_regression,
)
from srm_1d.nozzle import Nozzle
from srm_1d.openmotor_adapter import load_ric, convert_geometry, convert_nozzle
from tests._motor_fixtures import (
    hasegawa_propellant_1 as make_hasegawa_propellant_1,
    single_cylinder_geo, example_bates_geo, bates_motor_geo,
    MOTORS_DIR,
)


def make_single_cylinder(D_bore, D_outer, length, N_cells=None):
    """Adapter for v0.5.x test signature.

    The old factory took ``N_cells`` directly; the snapped builder takes
    ``target_propellant_cells``. They map 1:1 for a single inhibited
    cylinder where leading/trailing spacers add ≤2 cells.
    """
    target = N_cells if N_cells is not None else 50
    return single_cylinder_geo(D_bore, D_outer, length,
                               target_propellant_cells=target)


def make_example_bates():
    return example_bates_geo()


def make_bates_motor(D_bore, D_outer, L_segment, N_segments, spacing):
    return bates_motor_geo(D_bore, D_outer, L_segment, N_segments, spacing)


def _hasegawa_a_geo_and_nozzle():
    motor = load_ric(str(MOTORS_DIR / 'hasegawa_a.ric'))
    return convert_geometry(motor['grains']), convert_nozzle(motor['nozzle'])


def make_hasegawa_motor_A_geo():
    return _hasegawa_a_geo_and_nozzle()[0]


def make_hasegawa_motor_A_nozzle():
    return _hasegawa_a_geo_and_nozzle()[1]


def _init_geometry_arrays(geo):
    """Helper: compile geometry and initialize working arrays."""
    prop = make_hasegawa_propellant_1()
    ga = geo.compile_geometry_arrays()
    N = geo.N_cells
    D_port = ga['D_port'].copy()
    regress = ga['regress'].copy()
    A_port, C_burn, D_hyd = np.zeros(N), np.zeros(N), np.zeros(N)
    is_grain = np.zeros(N, dtype=np.bool_)
    endface_msource = np.zeros(N)
    P = np.full(N, 3e6)
    update_cell_geometry(
        regress, D_port, ga['x_centers'], geo.dx, N, ga['N_seg'], ga['cell_D_outer'],
        ga['seg_x_start'], ga['seg_length'],
        ga['seg_fwd_regression'], ga['seg_aft_regression'],
        ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
        ga['cell_segment_id'], P, prop.rho_propellant,
        *prop.tab_arrays(), len(prop.tabs),
        A_port, C_burn, D_hyd, is_grain, endface_msource,
        ga['cell_D_bore_init'], ga['cell_wall_web'],
        ga['cell_segment_type'], ga['cell_fmm_idx'],
        ga['fmm_offset'], ga['fmm_reg_flat'],
        ga['fmm_perim_flat'], ga['fmm_port_flat'],
    )
    return ga, prop, D_port, A_port, C_burn, D_hyd, is_grain, endface_msource, P


class TestFactoryFunctions:
    def test_hasegawa_A_dimensions(self):
        geo = make_hasegawa_motor_A_geo()
        assert geo.D_outer == pytest.approx(0.080)
        assert len(geo.segments) == 1
        assert geo.segments[0].D_bore_fwd == pytest.approx(0.040)

    def test_hasegawa_A_nozzle_throat(self):
        nozzle = make_hasegawa_motor_A_nozzle()
        assert nozzle.D_throat == pytest.approx(0.034)

    def test_bates_segment_count(self):
        geo = make_example_bates()
        assert len(geo.segments) == 4

    def test_bates_motor_length(self):
        geo = make_bates_motor(0.038, 0.070, 0.120, 4, 0.005)
        # 4 propellant segments + 3 internal gaps + a 3-cell leading and
        # 3-cell trailing buffer (MIN_BUFFER_CELLS in
        # build_snapped_geometry; introduced 2026-05-20 to keep the
        # throat boundary resolved by a multi-cell pressure gradient).
        # Buffer cells inherit the snapped dx so the exact contribution
        # depends on cells_target. With target=100 and L_prop=0.48m,
        # dx clamps to the 5mm min-gap, giving 6*5mm = 30mm of extra
        # buffer beyond the 1mm-each requested in the old single-cell
        # convention.
        dx = geo.dx
        expected_segments = 4 * 0.120
        expected_internal_gaps = 3 * dx        # one cell each
        expected_buffers = 2 * 3 * dx          # 3-cell lead + 3-cell trail
        expected = expected_segments + expected_internal_gaps + expected_buffers
        assert geo.L_motor == pytest.approx(expected, abs=1e-9)

    def test_solver_buffers_do_not_trigger_physical_length_warning(self, recwarn):
        motor = load_ric(str(MOTORS_DIR / 'machbusterNew.ric'))
        geo = convert_geometry(motor['grains'], target_propellant_cells=200)
        requested_length = sum(
            grain['properties']['length'] for grain in motor['grains']
        )
        snapped_length = sum(segment.length for segment in geo.segments)

        assert snapped_length == pytest.approx(requested_length, abs=1e-12)
        assert geo.L_motor - snapped_length == pytest.approx(6 * geo.dx)
        assert not recwarn

    def test_physical_length_warning_reports_real_grain_gap_snap(self):
        segments_spec = [{
            'D_bore_fwd': 0.040,
            'length': 0.500,
            'gap_after': 0.080,
        }]

        with pytest.warns(
            UserWarning, match=r'physical grain/gap span by -20\.00 mm',
        ):
            build_snapped_geometry(
                segments_spec, D_outer=0.080, target_propellant_cells=2,
            )

    def test_web_thickness(self):
        geo = make_single_cylinder(0.040, 0.080, 0.500)
        assert geo.web_thickness == pytest.approx(0.020)

    def test_propellant_volume_cylinder(self):
        geo = make_single_cylinder(0.040, 0.080, 0.500)
        V = geo.total_propellant_volume()
        expected = np.pi / 4.0 * (0.080**2 - 0.040**2) * 0.500
        assert V == pytest.approx(expected, rel=1e-6)


class TestOverlapMatching:
    def test_initial_grain_coverage_bates(self):
        """All 8 BATES segment boundaries should have grain coverage."""
        geo = make_example_bates()
        ga, prop, D_port, A_port, C_burn, D_hyd, is_grain, ef, P = \
            _init_geometry_arrays(geo)
        N, dx = geo.N_cells, geo.dx

        # Compute total grain coverage
        total_coverage = 0.0
        for i in range(N):
            k = ga['cell_segment_id'][i]
            if k < 0:
                continue
            x = ga['x_centers'][i]
            x_lo, x_hi = x - 0.5*dx, x + 0.5*dx
            x_fwd = ga['seg_x_start'][k]
            x_aft = ga['seg_x_start'][k] + ga['seg_length'][k]
            overlap = max(0, min(x_hi, x_aft) - max(x_lo, x_fwd))
            total_coverage += overlap

        total_seg_length = sum(s.length for s in geo.segments)
        assert total_coverage == pytest.approx(total_seg_length, rel=1e-6)

    def test_initial_C_burn_consistency(self):
        """Sum of C_burn×dx should equal π×D_bore×total_segment_length."""
        geo = make_example_bates()
        ga, prop, D_port, A_port, C_burn, D_hyd, is_grain, ef, P = \
            _init_geometry_arrays(geo)
        sum_Cburn_dx = np.sum(C_burn * geo.dx)
        D_bore = geo.segments[0].D_bore_fwd
        total_seg_length = sum(s.length for s in geo.segments)
        expected = np.pi * D_bore * total_seg_length
        assert sum_Cburn_dx == pytest.approx(expected, rel=1e-3)


class TestEndfaceInjection:
    def test_bates_8_faces(self):
        """4-segment BATES with no inhibition has 8 burning faces.

        The v0.6.0 linear-distribution kernel splits each face's mass
        over 2 adjacent cells (partition of unity). Intermediate
        single-cell inter-segment gaps receive contributions from
        BOTH bordering faces, so the cell count is
        8 faces × 2 cells − 3 shared gap cells = 13 distinct cells.
        """
        geo = make_example_bates()
        ga, prop, D_port, A_port, C_burn, D_hyd, is_grain, ef, P = \
            _init_geometry_arrays(geo)
        assert int(np.sum(ef > 0)) == 13

    def test_inhibited_ends_no_faces(self):
        """Single cylinder with inhibited ends should have 0 burning faces."""
        geo = make_single_cylinder(0.040, 0.080, 0.500)
        ga, prop, D_port, A_port, C_burn, D_hyd, is_grain, ef, P = \
            _init_geometry_arrays(geo)
        assert int(np.sum(ef > 0)) == 0


class TestMassConservation:
    """
    The critical test: bore mass produced during a full burn-to-completion
    should equal the initial propellant volume × density.
    Uses constant burn rate (no pressure coupling) to isolate geometry.
    """

    def _run_burnout(self, geo, dt=1e-4, max_steps=60000):
        prop = make_hasegawa_propellant_1()
        ga = geo.compile_geometry_arrays()
        N, dx = geo.N_cells, geo.dx
        D_port = ga['D_port'].copy()
        regress = ga['regress'].copy()
        A_port, C_burn, D_hyd = np.zeros(N), np.zeros(N), np.zeros(N)
        is_grain = np.zeros(N, dtype=np.bool_)
        endface_msource = np.zeros(N)
        r_total = np.full(N, 0.005)  # constant 5 mm/s
        P = np.full(N, 3e6)

        update_cell_geometry(
            regress, D_port, ga['x_centers'], dx, N, ga['N_seg'], ga['cell_D_outer'],
            ga['seg_x_start'], ga['seg_length'],
            ga['seg_fwd_regression'], ga['seg_aft_regression'],
            ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
            ga['cell_segment_id'], P, prop.rho_propellant,
            *prop.tab_arrays(), len(prop.tabs),
            A_port, C_burn, D_hyd, is_grain, endface_msource,
            ga['cell_D_bore_init'], ga['cell_wall_web'],
            ga['cell_segment_type'], ga['cell_fmm_idx'],
            ga['fmm_offset'], ga['fmm_reg_flat'],
            ga['fmm_perim_flat'], ga['fmm_port_flat'],
        )

        initial_mass = geo.total_propellant_volume() * prop.rho_propellant
        mass_produced = 0.0

        for step in range(max_steps):
            advance_bore_regression(
                regress, r_total, dt, N,
                ga['cell_wall_web'], ga['cell_segment_id'],
            )
            update_cell_geometry(
                regress, D_port, ga['x_centers'], dx, N, ga['N_seg'], ga['cell_D_outer'],
                ga['seg_x_start'], ga['seg_length'],
                ga['seg_fwd_regression'], ga['seg_aft_regression'],
                ga['seg_inhibit_fwd'], ga['seg_inhibit_aft'],
                ga['cell_segment_id'], P, prop.rho_propellant,
                *prop.tab_arrays(), len(prop.tabs),
                A_port, C_burn, D_hyd, is_grain, endface_msource,
                ga['cell_D_bore_init'], ga['cell_wall_web'],
                ga['cell_segment_type'], ga['cell_fmm_idx'],
                ga['fmm_offset'], ga['fmm_reg_flat'],
                ga['fmm_perim_flat'], ga['fmm_port_flat'],
            )
            mass_produced += prop.rho_propellant * np.sum(
                r_total * C_burn * is_grain
            ) * dx * dt
            n_active = int(np.sum(
                (D_port < ga['D_outer']) & (ga['cell_segment_id'] >= 0)
            ))
            if n_active == 0:
                break

        return mass_produced, initial_mass

    def test_cylinder_conservation(self):
        """Single cylinder: bore mass conserved to < 0.1%."""
        geo = make_single_cylinder(0.040, 0.080, 0.500, N_cells=50)
        produced, initial = self._run_burnout(geo)
        error_pct = abs(produced - initial) / initial * 100
        assert error_pct < 0.1, f"Conservation error {error_pct:.3f}%"

    def test_bates_conservation(self):
        """BATES: bore mass conserved to < 0.1%."""
        geo = make_example_bates()
        produced, initial = self._run_burnout(geo)
        error_pct = abs(produced - initial) / initial * 100
        assert error_pct < 0.1, f"Conservation error {error_pct:.3f}%"
