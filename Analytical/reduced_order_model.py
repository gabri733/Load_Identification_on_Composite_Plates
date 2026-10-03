"""
Uso:
    python inverse_fem_problem_v2.py --csv dati.csv --dt 0.002
    python inverse_fem_problem_v2.py --csv dati.csv
    python inverse_fem_problem_v2.py --csv dati.csv --n-poly 10

    # Con calibrazione:
    python inverse_fem_problem_v2.py --csv dati.csv \\
        --cal-csv dati.csv --cal-pmax 15000 --cal-phi -0.2 0 5 0
"""

import argparse, sys, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.linalg import eigh
import matplotlib.pyplot as plt

if hasattr(sys.stdout, 'reconfigure'):
    try: sys.stdout.reconfigure(encoding='utf-8')
    except: pass

# ===========================================================================
# ARGOMENTI
# ===========================================================================
parser = argparse.ArgumentParser()
parser.add_argument('--csv',        required=True)
parser.add_argument('--dt',         type=float, default=None)
parser.add_argument('--no-header',  action='store_true')
parser.add_argument('--t-start',    type=float, default=None)
parser.add_argument('--t-end',      type=float, default=None)
parser.add_argument('--k_m',     type= int, default=0)
parser.add_argument('--pmax-min',   type=float, default=1e3)
parser.add_argument('--pmax-max',   type=float, default=80e3)
parser.add_argument('--lambda-reg', type=float, default=0.5)
parser.add_argument('--T-imp',      type=float, default=0.1)
parser.add_argument('--n-poly',     type=int,   default=8,
                    help='Ordine polinomiale Ritz (default: 8, suggerito: 10 per più accuratezza)')
parser.add_argument('--grid-n',     type=int,   default=20,
                    help='Punti per asse nella grid search (default: 5)')
# --- CALIBRAZIONE ---
parser.add_argument('--cal-csv',  default=None)
parser.add_argument('--cal-pmax', type=float, default=None)
parser.add_argument('--cal-phi',  type=float, nargs=4, default=None,
                    metavar=('xs0','ys0','u','v'))
parser.add_argument('--cal-dt',   type=float, default=None)
args = parser.parse_args()

T_imp = args.T_imp
k_m = args.k_m
# ===========================================================================
# LETTURA CSV
# ===========================================================================
print(f"\nLettura: {args.csv}")
skip = 0 if args.no_header else 1
data = None
for delim in [',', ';', '\t']:
    try:
        tmp = np.loadtxt(args.csv, delimiter=delim, skiprows=skip)
        if tmp.ndim == 2 and tmp.shape[1] == 8:
            data = tmp; print(f"  Separatore: '{delim}'"); break
    except: continue
if data is None:
    print("ERRORE: CSV non leggibile con 8 colonne."); sys.exit(1)

N_rows = data.shape[0]

data_fix = data.copy()
print("  Segni CSV invariati (fisicamente corretti, nessun flip)")

# Stima dt automatica dal picco (carico sin^2 -> picco a T_imp/2)
if args.dt is None:
    rms_ch = np.sqrt(np.mean(data_fix**2, axis=0))
    ch_ref = np.argmax(rms_ch)
    idx_pk = np.argmax(np.abs(data_fix[:, ch_ref]))
    dt = (T_imp / 2) / max(idx_pk, 1)
    print(f"  dt stimato: {dt*1000:.3f} ms ({1/dt:.0f} Hz)  [picco FBG{ch_ref+1} a idx={idx_pk}]")
else:
    dt = args.dt
    print(f"  dt: {dt*1000:.3f} ms ({1/dt:.0f} Hz)")

t_full = np.arange(N_rows) * dt
i0 = 0       if args.t_start is None else max(0,      int(args.t_start / dt))
i1 = N_rows  if args.t_end   is None else min(N_rows, int(args.t_end   / dt) + 1)
eps_meas = data_fix[i0:i1, :].T   # (8, N_t)
t_meas  = t_full[i0:i1]
N_t     = eps_meas.shape[1]
print(f"  Finestra: [{t_meas[0]:.4f}, {t_meas[-1]:.4f}] s  ({N_t} campioni)")


# ===========================================================================
# PARAMETRI FISICI E ASSEMBLAGGIO MODELLO
# ===========================================================================
E_plain, E_twill, nu_f, G_f = 53e9, 59e9, 0.20, 3.3e9
rho = 960.0; L = 0.845; a = L/2; Rmin = 2

layers = [('plain',45,0.23e-3),('twill',0,0.65e-3),('twill',0,0.65e-3),
          ('twill',0,0.65e-3),('twill',0,0.65e-3),('twill',45,0.65e-3),('plain',45,0.23e-3)]

FBG_pos = np.array([
    [ 0.000, -0.3,    0],
    [ 0.3,    0.000, 90],
    [ 0.000,  0.3,    0],
    [-0.3,    0.000, 90],
    [ 0.17,-0.17, 45],
    [-0.17,-0.17, 135],
    [-0.17, 0.17, 225],
    [ 0.17, 0.17, 315],
])

def Q_loc(mat):
    E = E_plain if mat=='plain' else E_twill; d = 1-nu_f**2
    return np.array([[E/d,nu_f*E/d,0],[nu_f*E/d,E/d,0],[0,0,G_f]])

def rot_Q(Q, th):
    t=np.radians(th); c,s=np.cos(t),np.sin(t)
    T   =np.array([[c**2,s**2, 2*c*s],[s**2,c**2,-2*c*s],[-c*s, c*s,c**2-s**2]])
    Tinv=np.array([[c**2,s**2,-2*c*s],[s**2,c**2, 2*c*s],[ c*s,-c*s,c**2-s**2]])
    return Tinv@Q@T

z = [-sum(t for _,_,t in layers)/2]
for _,_,t in layers: z.append(z[-1]+t)
z = np.array(z); h_tot = z[-1]-z[0]

D = np.zeros((3,3))
for k,(mat,ang,_) in enumerate(layers):
    D += rot_Q(Q_loc(mat),ang)*(z[k+1]**3-z[k]**3)/3
D11,D22,D12,D16,D26,D66 = D[0,0],D[1,1],D[0,1],D[0,2],D[1,2],D[2,2]
print(f"\nMatrice D: D11={D11:.2f}  D22={D22:.2f}  D12={D12:.2f}  D66={D66:.2f}")

# MIGLIORAMENTO MINORE: z_fbg = coordinata z del layer sensori (interfaccia layer6/layer8)
# I sensori FBG si trovano alla quota z[6] (top del 6° strato = bottom del plain superiore)
z_fbg = z[-1]      # 1.625 mm (vs z[-1]=1.855 mm della fibra esterna)
z_surf = z[-1]    # h/2 = fibra esterna (usato come riferimento)
print(f"  z_fbg={z_fbg*1000:.3f} mm  z_surf={z_surf*1000:.3f} mm  "
      f"(differenza: {(z_surf-z_fbg)/z_surf*100:.1f}%)")

N_poly = args.n_poly
idx = [(i,j) for i in range(N_poly) for j in range(N_poly)]
n_modes = len(idx)

def phi_b(i,j,x,y): return (x**2-a**2)*(y**2-a**2)*x**i*y**j

def dx(i,j,x,y):
    gy=(y**2-a**2)*y**j; f=(x**2-a**2); df=2*x
    gx=x**i
    dgx = i*x**(i-1)       if i>=1 else np.zeros_like(x)
    return gy*(df*gx + f*dgx)

def dy(i,j,x,y):
    gx=(x**2-a**2)*x**i; f=(y**2-a**2); df=2*y
    gy=y**j
    dgy = j*y**(j-1)       if j>=1 else np.zeros_like(y)
    return gx*(df*gy+f*dgy)

def dxx(i,j,x,y):
    gy=(y**2-a**2)*y**j; f=(x**2-a**2); df=2*x; d2f=2
    gx=x**i
    dgx = i*x**(i-1)       if i>=1 else np.zeros_like(x)
    d2gx= i*(i-1)*x**(i-2) if i>=2 else np.zeros_like(x)
    return gy*(d2f*gx+2*df*dgx+f*d2gx)

def dyy(i,j,x,y):
    gx=(x**2-a**2)*x**i; f=(y**2-a**2); df=2*y; d2f=2
    gy=y**j
    dgy = j*y**(j-1)       if j>=1 else np.zeros_like(y)
    d2gy= j*(j-1)*y**(j-2) if j>=2 else np.zeros_like(y)
    return gx*(d2f*gy+2*df*dgy+f*d2gy)

def dxy(i,j,x,y):
    dfx = 2*x*x**i + (x**2-a**2)*(i*x**(i-1) if i>=1 else np.zeros_like(x))
    dfy = 2*y*y**j + (y**2-a**2)*(j*y**(j-1) if j>=1 else np.zeros_like(y))
    return dfx*dfy

print(f"\nAssemblaggio modello ({n_modes} funzioni di forma, N_poly={N_poly})...")
n_gauss = 50
xi, wi = np.polynomial.legendre.leggauss(n_gauss)
xg, wg = xi*a, wi*a
XX, YY = np.meshgrid(xg, xg, indexing='ij')
WW = np.outer(wg, wg)

PHI   = np.array([phi_b(i,j,XX,YY) for i,j in idx])
PHIx_ap = np.array([dx(i,j,a,YY)   for i,j in idx])
PHIx_am = np.array([dx(i,j,-a,YY)   for i,j in idx])
PHIy_ap = np.array([dy(i,j,XX,a)   for i,j in idx])
PHIy_am = np.array([dy(i,j,XX,-a)   for i,j in idx])
PHIxx = np.array([dxx(i,j,XX,YY)   for i,j in idx])
PHIyy = np.array([dyy(i,j,XX,YY)   for i,j in idx])
PHIxy = np.array([dxy(i,j,XX,YY)   for i,j in idx])

M = np.zeros((n_modes,n_modes))
K = np.zeros((n_modes,n_modes))
for p in range(n_modes):
    for q in range(p, n_modes):
        M[p,q] = rho*h_tot*np.sum(WW*PHI[p]*PHI[q])
        M[q,p] = M[p,q]
        k = (D11*np.sum(WW*PHIxx[p]*PHIxx[q])
           + D22*np.sum(WW*PHIyy[p]*PHIyy[q])
           + D12*np.sum(WW*(PHIxx[p]*PHIyy[q]+PHIyy[p]*PHIxx[q]))
           + 4*D66*np.sum(WW*PHIxy[p]*PHIxy[q])
           + 2*D16*np.sum(WW*(PHIxx[p]*PHIxy[q]+PHIxx[q]*PHIxy[p]))
           + 2*D26*np.sum(WW*(PHIyy[p]*PHIxy[q]+PHIyy[q]*PHIxy[p]))
           + k_m *(np.sum(WW*(PHIx_ap[p]*PHIx_ap[q]+PHIx_am[p]*PHIx_am[q]
                             +PHIy_ap[p]*PHIy_ap[q]+PHIy_am[p]*PHIy_am[q])))) 
        K[p,q] = k
        K[q,p] = k

# ===========================================================================
# MIGLIORAMENTO 2: PRE-DIAGONALIZZAZIONE M, K
# Risolve il problema agli autovalori una sola volta.
# Il forward model usa poi ODE SDOF disaccoppiate (molto più veloci).
# ===========================================================================
print("  Diagonalizzazione M, K...")
eigvals, Phi_modal = eigh(K, M)
# Filtra autovalori negativi (instabilità numerica per modi rigidi)
omega2 = np.maximum(eigvals, 0.0)
omega_r = np.sqrt(omega2)          # pulsazioni proprie (rad/s)
freqs_r = omega_r / (2*np.pi)      # frequenze proprie (Hz)

print(f"  Frequenze proprie (primi 6) [Hz]: {freqs_r[:6].round(2)}")
print(f"  f1 Ritz = {freqs_r[0]:.2f} Hz  (FEM: 50.89 Hz, errore: "
      f"{abs(freqs_r[0]-50.89)/50.89*100:.1f}%)")

# Matrice C sensori in coordinate modali: C_modal = C @ Phi_modal
# C (8, n_modes) @ Phi_modal (n_modes, n_modes) -> (8, n_modes)
C_phys = np.zeros((8, n_modes))
for fbg in range(8):
    xf,yf,th = FBG_pos[fbg]; tr=np.radians(th); ca,sa=np.cos(tr),np.sin(tr)
    pxx = np.array([dxx(i,j,float(xf),float(yf)) for i,j in idx])
    pyy = np.array([dyy(i,j,float(xf),float(yf)) for i,j in idx])
    pxy = np.array([dxy(i,j,float(xf),float(yf)) for i,j in idx])
    # Uso z_fbg (quota sensore) invece di z_surf (fibra esterna)
    C_phys[fbg] = z_fbg*(pxx*ca**2 + pyy*sa**2 + 2*pxy*ca*sa)

C_modal = C_phys @ Phi_modal   # (8, n_modes) in coordinate modali

# ===========================================================================
# FORWARD MODEL — versione modale accelerata
# ===========================================================================
def run_forward(phi_params, t_arr):

    xs0, ys0, ud, vd = phi_params

    # Pre-calcola le forze generalizzate Ritz su tutta la finestra temporale
    F_ritz = np.zeros((n_modes, len(t_arr)))
    for k, t in enumerate(t_arr):
        if 0 < t < T_imp:
            xs = xs0 + ud*t; ys = ys0 + vd*t
            R  = np.sqrt((XX-xs)**2 + (YY-ys)**2 + Rmin**2)
            pt = (Rmin/R)*np.sin(np.pi*t/T_imp)**2
            F_ritz[:, k] = np.einsum('ijk,jk->i', PHI, WW*pt)

    # Proietta in coordinate modali
    f_modal = Phi_modal.T @ F_ritz   # (n_modes, N_t)

    # Integra ogni oscillatore SDOF separatamente con RK45
    def rhs_modal(t, state):
        q  = state[:n_modes]
        qd = state[n_modes:]
        k_t = min(np.searchsorted(t_arr, t), len(t_arr)-1)
        qdd = f_modal[:, k_t] - omega2*q
        return np.concatenate([qd, qdd])

    sol = solve_ivp(rhs_modal, (t_arr[0], t_arr[-1]),
                    np.zeros(2*n_modes),
                    t_eval=t_arr, method='RK45',
                    rtol=1e-6, atol=1e-9)

    return (C_modal @ sol.y[:n_modes, :]) * 1e6   # (8, N_t)

# ===========================================================================
# STIMA p_max — Variable Projection (formula con pesi lineari)
# ===========================================================================
energia = np.array([max(np.dot(eps_meas[i], eps_meas[i]), 1e-10) for i in range(8)])
w       = 1.0 / energia   # (8,)

print(f"\nEnergia segnali per sensore [ue^2]:")
for i in range(8):
    print(f"  FBG{i+1}: {energia[i]:.2f}  w={w[i]:.2e}")

J_ref = np.sum(w * energia)
print(f"\nJ* riferimento (fit nullo) = {J_ref:.4f}  (teorico = {8.0})")

def stima_pmax(eps_tilde, eps_m):
    num = sum(w[i] * np.dot(eps_tilde[i], eps_m[i])   for i in range(8))
    den = sum(w[i] * np.dot(eps_tilde[i], eps_tilde[i]) for i in range(8))
    
    # Se il modello è quasi nullo (es. in 0,0), diamo un errore altissimo
    if den < 1e-20:
        return 0.0, 1e10  # Penalità massiva invece di J_ref
    
    p_star = num / den
    
    # Residuo classico
    J_pure = (sum(w[i] * np.dot(eps_m[i], eps_m[i]) for i in range(8)) - num**2 / den)
    
    # Normalizzazione correttiva: J pesato sull'energia spiegata
    # Questo impedisce a J di scendere solo perché il modello si spegne
    J_norm = J_pure / (num**2 / den + 1e-10) 
    
    return p_star, J_pure # Restituiamo J_pure ma possiamo usare J_norm per guidare

def J_star_reg(phi_params, t_arr):
    eps_t  = run_forward(phi_params, t_arr)
    p_est, J = stima_pmax(eps_t, eps_meas)
    # Regolarizzazione soft su p_max
    if p_est < 0:
        J += args.lambda_reg * 10.0
    elif p_est < args.pmax_min:
        J += args.lambda_reg * ((args.pmax_min - p_est) / args.pmax_min)**2
    elif p_est > args.pmax_max:
        J += args.lambda_reg * ((p_est - args.pmax_max) / args.pmax_max)**2
    return J, p_est

eps_zero = np.zeros((8, N_t))
_, J_zero = stima_pmax(eps_zero, eps_meas)
print(f"Sanity check J*(eps_tilde=0) = {J_zero:.4f}  (atteso ~{J_ref:.1f})")

# ===========================================================================
# CALIBRAZIONE (opzionale)
# ===========================================================================
k_cal = 1.0
cal_info = "nessuna (p_max Ritz non corretto)"

if args.cal_csv is not None:
    if args.cal_pmax is None or args.cal_phi is None:
        print("\nATTENZIONE: --cal-csv richiede anche --cal-pmax e --cal-phi. "
              "Calibrazione saltata.")
    else:
        print(f"\n--- CALIBRAZIONE ---")
        cal_data = None
        for delim in [',', ';', '\t']:
            try:
                tmp = np.loadtxt(args.cal_csv, delimiter=delim,
                                 skiprows=0 if args.no_header else 1)
                if tmp.ndim == 2 and tmp.shape[1] == 8:
                    cal_data = tmp; break
            except: continue
        if cal_data is None:
            print("  ERRORE: CSV calibrazione non leggibile. Saltata.")
        else:
            cal_fix  = cal_data.copy()
            cal_dt   = args.cal_dt if args.cal_dt is not None else dt
            t_cal    = np.arange(cal_fix.shape[0]) * cal_dt
            eps_cal  = (cal_fix.T - cal_fix[:1, :].T)

            eps_tilde_cal = run_forward(tuple(args.cal_phi), t_cal)
            en_cal = np.array([max(np.dot(eps_cal[i], eps_cal[i]), 1e-10) for i in range(8)])
            w_cal  = 1.0 / en_cal
            num_cal = sum(w_cal[i]*np.dot(eps_tilde_cal[i], eps_cal[i])  for i in range(8))
            den_cal = sum(w_cal[i]*np.dot(eps_tilde_cal[i], eps_tilde_cal[i]) for i in range(8))
            if abs(den_cal) < 1e-30:
                print("  ERRORE: forward cal produce risposta nulla. Saltata.")
            else:
                p_max_ritz_cal = num_cal / den_cal
                k_cal  = args.cal_pmax / p_max_ritz_cal
                cal_info = (f"k_cal = {args.cal_pmax:.0f} / {p_max_ritz_cal:.0f} "
                            f"= {k_cal:.4f}")
                print(f"  p_max Ritz (cal): {p_max_ritz_cal:.0f} Pa")
                print(f"  p_max noto (FEM): {args.cal_pmax:.0f} Pa")
                print(f"  k_cal = {k_cal:.4f}  "
                      f"({'<1: modello sovrastima' if k_cal<1 else '>1: modello sottostima'})")

# ===========================================================================
# MIGLIORAMENTO 1: GRID SEARCH AMPLIATA
# Differenze rispetto alla v1:
#   - ys0: [-1.2, 1.2] invece di [-0.7, 0.2]  -> copre Moto2 (ys0=1.0)
#   - v:   [0.0, 15.0] invece di [2.0, 15.0]  -> copre Moto1/2 (v=0)
#   - griglia più densa: 6x6x5x5 invece di 4x4x4x4
# ===========================================================================
print("\n--- GRID SEARCH (v2: range ampliato, griglia più densa) ---")
ng = args.grid_n

xs0_g  = np.linspace(-0.5,  0, ng)
# ys0_g  = np.linspace(-0.5,  0, ng-1)   # range ampliato per coprire Moto2
ys0_g  = [0]
udot_g = np.linspace( 0.0, 10.0, ng)
# vdot_g = np.linspace( 0.0, 10.0, ng-1)     # inizia da 0 (non da 2) per Moto1/2
vdot_g = [0]

combos, n_skip = [], 0
for xs0 in xs0_g:
    for ys0 in ys0_g:
        for ud in udot_g:
            for vd in vdot_g:
                xs_f = xs0 + ud*T_imp; ys_f = ys0 + vd*T_imp
                if (min(xs0,xs_f) <= a and max(xs0,xs_f) >= -a
                        and min(ys0,ys_f) <= a+0.5 and max(ys0,ys_f) >= -a-0.5):
                    combos.append((xs0, ys0, ud, vd))
                else:
                    n_skip += 1

print(f"  Combinazioni valide: {len(combos)}  (scartate: {n_skip})")
print(f"  p_max atteso: [{args.pmax_min:.0f}, {args.pmax_max:.0f}] Pa  lambda={args.lambda_reg}")

best_J, best_phi, best_p = np.inf, None, None
for n, phi_t in enumerate(combos, 1):
    J_val, p_est = J_star_reg(phi_t, t_meas)
    if J_val < best_J:
        best_J, best_phi, best_p = J_val, phi_t, p_est
        print(f"  [{n:4d}/{len(combos)}] MINIMO  "
              f"xs0={phi_t[0]:+.2f} ys0={phi_t[1]:+.2f} "
              f"u={phi_t[2]:.1f} v={phi_t[3]:+.1f}  "
              f"J*={J_val:.4f}/{J_ref:.1f}  p_max={p_est:.0f} Pa")
    elif n % max(1, len(combos)//8) == 0:
        print(f"  [{n:4d}/{len(combos)}] J*_best={best_J:.4f}  p_best={best_p:.0f} Pa")

print(f"\nGrid completata. Best: J*={best_J:.4f}  p_max={best_p:.0f} Pa")

# ===========================================================================
# NELDER-MEAD (local refinement — invariato rispetto v1)
# ===========================================================================
print("\n--- NELDER-MEAD ---")
res_nm = minimize(lambda p: J_star_reg(tuple(p), t_meas)[0],
                  list(best_phi),
                  method='Nelder-Mead',
                  options={'xatol':1e-3, 'fatol':1e-4, 'maxiter':500, 'adaptive':True})
phi_nm = tuple(res_nm.x)
print(f"  Convergenza: {res_nm.success}  iter={res_nm.nit}  nfev={res_nm.nfev}")
print(f"  phi_nm: xs0={phi_nm[0]:.4f}  ys0={phi_nm[1]:.4f}  "
      f"u={phi_nm[2]:.4f}  v={phi_nm[3]:.4f}")

# ===========================================================================
# MIGLIORAMENTO 4: RAFFINAMENTO L-BFGS-B dopo Nelder-Mead
# Bounds stretti intorno alla soluzione Nelder-Mead per convergenza finale.
# ===========================================================================
print("\n--- L-BFGS-B (raffinamento finale) ---")
delta = 1.0   # finestra di ricerca ±1 m (posizione) / ±3 m/s (velocità)
xs0_nm, ys0_nm, u_nm, v_nm = phi_nm
bounds = [
    (xs0_nm - delta,    xs0_nm + delta),
    (ys0_nm - delta,    ys0_nm + delta),
    (max(0.0, u_nm - 3), u_nm + 3),
    (max(0.0, v_nm - 3), v_nm + 3),
]
res_lbfgs = minimize(lambda p: J_star_reg(tuple(p), t_meas)[0],
                     list(phi_nm),
                     method='L-BFGS-B',
                     bounds=bounds,
                     options={'ftol':1e-7, 'gtol':1e-5, 'maxiter':200})
phi_hat = tuple(res_lbfgs.x)
print(f"  Convergenza: {res_lbfgs.success}  iter={res_lbfgs.nit}  nfev={res_lbfgs.nfev}")
print(f"  phi_hat: xs0={phi_hat[0]:.4f}  ys0={phi_hat[1]:.4f}  "
      f"u={phi_hat[2]:.4f}  v={phi_hat[3]:.4f}")

# ===========================================================================
# RECUPERO p_max FINALE
# ===========================================================================
print("\n--- RECUPERO p_max ---")
eps_tilde_hat     = run_forward(phi_hat, t_meas)
p_max_hat, J_fin  = stima_pmax(eps_tilde_hat, eps_meas)
p_max_cal         = p_max_hat * k_cal
eps_rec           = p_max_hat * eps_tilde_hat

print(f"  p_max Ritz        = {p_max_hat:.2f} Pa")
if k_cal != 1.0:
    print(f"  k_cal             = {k_cal:.4f}")
    print(f"  p_max calibrato   = {p_max_cal:.2f} Pa")
print(f"  J* finale = {J_fin:.4f} / {J_ref:.1f}  ({J_fin/J_ref*100:.1f}% del fit nullo)")

print(f"\n  Errore RMS per sensore:")
rms_vals = []
for i in range(8):
    rms  = np.sqrt(np.mean((eps_meas[i]-eps_rec[i])**2))
    ref  = np.sqrt(np.mean(eps_meas[i]**2)) + 1e-10
    rms_vals.append(rms)
    flag = " <<" if rms/ref > 0.3 else ""
    print(f"    FBG{i+1}: {rms:.3f} ue  ({rms/ref*100:.1f}%){flag}")

# ===========================================================================
# RISULTATI
# ===========================================================================
print("\n" + "="*60)
print("  RISULTATI")
print("="*60)
print(f"  xs0   = {phi_hat[0]:+.4f} m")
print(f"  ys0   = {phi_hat[1]:+.4f} m")
print(f"  u     = {phi_hat[2]:+.4f} m/s")
print(f"  v     = {phi_hat[3]:+.4f} m/s")
print(f"  p_max Ritz      = {p_max_hat:.2f} Pa")
if k_cal != 1.0:
    print(f"  k_cal           = {k_cal:.4f}  [{cal_info}]")
    print(f"  p_max calibrato = {p_max_cal:.2f} Pa")
print(f"  J*    = {J_fin:.4f} / {J_ref:.1f}  (0=perfetto, {J_ref:.0f}=nullo)")
print("="*60)

# ===========================================================================
# PLOT
# ===========================================================================
fig, axes = plt.subplots(4, 2, figsize=(14, 12))
pmax_title = (f"p_max_Ritz={p_max_hat:.0f} Pa  →  p_max_cal={p_max_cal:.0f} Pa  (k={k_cal:.3f})"
              if k_cal != 1.0 else f"p_max={p_max_hat:.0f} Pa")
fig.suptitle(
    f"Identificazione carico — v2 (K corretta + grid ampliata)\n"
    f"xs0={phi_hat[0]:.3f} m  ys0={phi_hat[1]:.3f} m  "
    f"u={phi_hat[2]:.2f} m/s  v={phi_hat[3]:.2f} m/s\n"
    f"{pmax_title}    J*={J_fin:.4f}/{J_ref:.1f}", fontsize=10)

for i in range(8):
    ax = axes[i//2, i%2]
    rms = rms_vals[i]; ref = np.sqrt(np.mean(eps_meas[i]**2)) + 1e-10
    col = 'red' if rms/ref > 0.3 else 'black'
    ax.plot(t_meas, eps_meas[i], 'k-',  lw=1.2, label='Misurato',    alpha=0.85)
    ax.plot(t_meas, eps_rec[i],  'r--', lw=1.5, label='Ricostruito', alpha=0.90)
    ax.set_title(f'FBG {i+1}   RMS={rms:.2f} ue ({rms/ref*100:.0f}%)', color=col)
    ax.set_xlabel('t [s]'); ax.set_ylabel('Strain [ue]')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
# Salva nella directory corrente (o in quella del CSV se scrivibile)
csv_dir = os.path.dirname(os.path.abspath(args.csv))
out = os.path.join(csv_dir, 'inverse_result_v2.png')
try:
    plt.savefig(out, dpi=150, bbox_inches='tight')
except OSError:
    out = os.path.join(os.getcwd(), 'Model_comparison.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"\nFigura salvata: {out}")
plt.show()