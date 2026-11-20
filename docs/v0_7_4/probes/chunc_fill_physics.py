"""Fill-transient physics audit for Chunc.
(1) Geometry & area ratios -> is high bore Mach even expected (port-to-throat)?
(2) Mass balance over the fill -> generation vs nozzle-out vs accumulation; is
    the nozzle choked; is the igniter/combustion gas production too fast?
(3) Pressurization rate vs experiment (~10 ms to 8.5 MPa).
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(r"c:\Users\ejarocki\Documents\Rocketry\Code Stuff\Erosive Burning Solver\srm_1d")
sys.path.insert(0, str(ROOT))
from srm_1d.openmotor_adapter import run_from_ric

KN = dict(roughness=32e-6, kappa=0.44, pyrogen="bpnv", pyrogen_mass=None,
          T_ignition=756.0, P_cutoff=0.01e6, t_max=1.3,
          snapshot_interval=0.0005, print_interval=0.0, verbose=False)
res = run_from_ric(str(ROOT/"motors"/"machbusterNew.ric"), **KN)[0]
snaps = res["snapshots"]
dx = float(res["dx"])
t_hist = np.asarray(res["time"]); P_hist = np.asarray(res["P_head"])
mdot_ig = np.asarray(res["mdot_ig"]); mflow = np.asarray(res["massflow"])
Dthr = np.asarray(res["D_throat"])
g0 = snaps[0]["is_grain"]
x = snaps[0]["x"]

# --- (1) geometry ---
s0 = snaps[0]
Dp0 = s0["D_port"][g0]
A_throat = np.pi/4*Dthr[0]**2
print("=== Geometry (t=0, min port) ===")
print(f"  grain cells={g0.sum()}, port length={g0.sum()*dx*1e3:.0f} mm")
print(f"  D_port: head={Dp0[0]*1e3:.1f}  mid={Dp0[len(Dp0)//2]*1e3:.1f}  "
      f"aft={Dp0[-1]*1e3:.1f} mm  (fore-taper if aft>head)")
A_head=np.pi/4*Dp0[0]**2; A_aft=np.pi/4*Dp0[-1]**2
print(f"  A_port: head={A_head*1e6:.1f}  aft={A_aft*1e6:.1f} mm^2; "
      f"A_throat={A_throat*1e6:.1f} mm^2 (Dthr={Dthr[0]*1e3:.1f} mm)")
print(f"  A_port_head/A_throat = {A_head/A_throat:.2f}, A_port_aft/A_throat = {A_aft/A_throat:.2f}")
print(f"  (isentropic: A/A*=1.0 at M=1; 1.34 at M=0.5; 1.7 at M=0.36; ~2.0 at M=0.3)")

# nozzle choke threshold: P_head/P_amb > 1/crit_pr ~1.8 (gamma~1.2)
print("\n=== Mass balance over the fill (nozzle choked once P_head>~0.18 MPa) ===")
print(f"{'t_ms':>5} {'P0_MPa':>6} {'maxM':>5} {'ign_g/s':>8} {'noz_g/s':>8} "
      f"{'Mbore_g':>8} {'dMb/dt':>8} {'prop_g/s':>9}")
Mbore_prev=None; t_prev=None
for s in snaps:
    if s["t"] > 0.012: break
    A = np.pi/4*s["D_port"]**2
    Mbore = float(np.sum(s["rho"]*A*dx))     # total bore gas mass [kg]
    j = int(np.argmin(np.abs(t_hist - s["t"])))
    ig = mdot_ig[j]; noz = mflow[j]
    if Mbore_prev is not None and s["t"]>t_prev:
        dMb = (Mbore-Mbore_prev)/(s["t"]-t_prev)
    else:
        dMb = 0.0
    prop = dMb + noz - ig                     # generation closing the balance
    maxM = np.abs(s["Mach"][g0]).max()
    print(f"{s['t']*1e3:5.1f} {s['P'][g0][0]/1e6:6.2f} {maxM:5.2f} "
          f"{ig*1e3:8.2f} {noz*1e3:8.2f} {Mbore*1e3:8.3f} {dMb*1e3:8.2f} {prop*1e3:9.2f}")
    Mbore_prev=Mbore; t_prev=s["t"]

# --- (3) pressurization rate ---
print("\n=== Pressurization rate vs experiment ===")
for thr in (1,2,4,6,8):
    idx = np.argmax(P_hist/1e6 >= thr)
    tt = t_hist[idx]*1e3 if (P_hist/1e6>=thr).any() else float('nan')
    print(f"  P_head reaches {thr} MPa at {tt:.1f} ms")
print(
    "  Historical Chunc trace: starts at 0.010 s; event-zero/calibration/window "
    "review pending (not a registered no-spike score)"
)
