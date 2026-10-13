#!/usr/bin/env python3
"""Rigorous (arb ball arithmetic, python-flint) dual certificate for the capacity-3 energy.

Functional, rho = 1/2, lam = log q:
  B0(x)    = (x+2)log(x+2) - x log x - 2
  h_rho(x) = x log x - (x-rho)log(x-rho) - rho log rho
  V(x)     = B0(x) - 3 h_rho(x) + lam x
           = (x+2)log(x+2) - 4 x log x + 3 (x-rho)log(x-rho) + lam x - 2 + 3 rho log rho
  E(mu)    = 3/2 + int V dmu - iint log|x-y| dmu dmu,   A(lam) = inf over mu in K_3.
Dual bound, for any probability density sigma with finite energy:
  A(lam) >= 3/2 + L(sigma,sigma) + ell - 3 int_{1/2}^inf (ell - W)_+ dt,
  W = V - 2 L_sigma,  L_sigma(t) = int log|t-s| sigma(s) ds.
(uses only -L(mu-sigma,mu-sigma) >= 0 and 0 <= dmu/dx <= 3.)
sigma: piecewise constant, density n_i / 2^40 on cell [1/2 + i/H, 1/2 + (i+1)/H], i < N.
Primal upper bound A <= E(sigma) (sigma feasible: 0 <= density <= 3, mass exactly 1) is also
computed in closed form.

Upper bound for int (ell-W)_+ on [1/2, T0], T0 = 1/2 + N/H: pieces [a,b] of length 1/(H K)
(cell breakpoints are piece endpoints, never interior).  On each piece W = Cv + Kc with Cv convex,
Kc concave (classification term by term, chi''(u) = 1/u), hence
  W(t) >= Kc(a) + (Kc(b)-Kc(a))(t-a)/s + Cv(m) + Cv'(m)(t-m)  (chord + tangent at midpoint m),
a linear minorant; int over the piece of (ell - minorant)_+ is evaluated exactly (upper endpoints).
Tail t >= T0: W(t) >= g(t) := V(t) - 2 log(t - rho) (since supp sigma in [rho,T0]);
g'(t) >= lam - (3 rho + 2)/(t - rho) > 0 for t >= T0 when lam (T0 - rho) > 7/2 (checked), and
g(T0) > ell is checked in arb.  So (ell - W)_+ = 0 on [T0, inf).
"""
import sys, json, math, time
from fractions import Fraction as Fr
import flint
from flint import arb, fmpq, ctx
ctx.prec = 128
# Proven rational upper bounds; energy/verify_tau.py checks them exactly.
TAU = {67: '6.272107', 70: '6.307708'}
DEN_BITS = 40
DEN = 1 << DEN_BITS
RHO = Fr(1, 2)

def A(x):  # exact rational -> arb (dyadics here, exact)
    x = Fr(x); return arb(fmpq(x.numerator, x.denominator))

def xlogx(x):  # x exact rational >= 0
    x = Fr(x)
    if x == 0: return arb(0)
    ax = A(x); return ax*ax.log()

def chi(u):  # u log|u| - u, chi(0)=0
    u = Fr(u)
    if u == 0: return arb(0)
    au = A(u); return au*abs(au).log() - au

def rationalize(dens, H):
    N = len(dens)
    n = [min(max(int(round(v*DEN)), 0), 3*DEN) for v in dens]
    target = H*DEN
    deficit = target - sum(n)
    # put the deficit in the free (0 < n < 3 DEN) cell with the largest room in the needed direction
    free = [i for i in range(N) if 0 < n[i] < 3*DEN]
    if deficit > 0: i0 = max(free, key=lambda i: 3*DEN - n[i])
    else: i0 = max(free, key=lambda i: n[i])
    n[i0] += deficit
    assert 0 <= n[i0] <= 3*DEN and sum(n) == target
    return n

def certify(a, H, N, n, ell, K, tau=None):
    t0 = time.time()
    lam = arb(a).log()
    h = Fr(1, H); rho = RHO
    assert len(n) == N and sum(n) == H*DEN and all(0 <= v <= 3*DEN for v in n)
    sig = [Fr(v, DEN) for v in n]
    # ---------- L(sigma,sigma) exact closed form ----------
    def psi(k):  # psi(t) = t^2/2 log|t| - 3/4 t^2 at integer t
        k = abs(k)
        return (arb(k*k)/2*arb(k).log() if k > 1 else arb(0)) - arb(3*k*k)/4
    psiv = [psi(k) for k in range(N+1)]
    phi = [psiv[k+1] - 2*psiv[k] + (psiv[k-1] if k >= 1 else psiv[1]) for k in range(N)]
    Ssum = arb(0)
    for k in range(N):
        Ak = sum(n[i]*n[i+k] for i in range(N-k))  # integer, scale DEN^2
        if Ak: Ssum += (1 if k == 0 else 2)*phi[k]*arb(fmpq(Ak, DEN*DEN))
    Lss = A(h).log() + Ssum*A(h*h)
    # ---------- primal upper bound E(sigma) ----------
    def P(x):  # antiderivative of x log x, P(0)=0
        x = Fr(x)
        if x == 0: return arb(0)
        ax = A(x); return ax*ax/2*ax.log() - ax*ax/4
    const = arb(-2) + 3*A(rho)*A(rho).log()
    intV = arb(0)
    for i in range(N):
        if n[i] == 0: continue
        lo = rho + i*h; hi = lo + h
        iv = (P(hi+2) - P(lo+2)) - 4*(P(hi) - P(lo)) + 3*(P(hi-rho) - P(lo-rho)) \
             + lam*A((hi*hi - lo*lo)/2) + const*A(h)
        intV += A(sig[i])*iv
    Eup = arb(3)/2 + intV - Lss
    # ---------- W pieces ----------
    c = [rho + j*h for j in range(N+1)]
    d = [(sig[j] if j < N else 0) - (sig[j-1] if j >= 1 else 0) for j in range(N+1)]
    J = [j for j in range(N+1) if d[j] != 0]
    s = Fr(1, H*K)
    T0 = rho + N*h
    ell_a = A(ell)
    I_up = arb(0); I_mid_est = 0.0; npos = 0
    npieces = N*K
    for p in range(npieces):
        a_ = rho + p*s; b_ = a_ + s; m_ = a_ + s/2
        # V parts
        Kc_a = -4*xlogx(a_); Kc_b = -4*xlogx(b_)
        Am = A(m_)
        Cv_m = xlogx(m_+2) + 3*xlogx(m_-rho) + lam*Am + const
        Cv_d = (A(m_+2)).log() + 1 + 3*((A(m_-rho)).log() + 1) + lam
        for j in J:
            cj = c[j]; coef = -2*d[j]
            if cj <= a_: convex = (d[j] < 0)
            else:
                assert cj >= b_; convex = (d[j] > 0)
            if convex:
                u = m_ - cj
                Cv_m += A(coef)*chi(u)
                Cv_d += A(coef)*abs(A(u)).log()
            else:
                Kc_a += A(coef)*chi(a_ - cj); Kc_b += A(coef)*chi(b_ - cj)
        lo_a = Kc_a + Cv_m - Cv_d*A(s/2)
        lo_b = Kc_b + Cv_m + Cv_d*A(s/2)
        Ea = (ell_a - lo_a).upper(); Eb = (ell_a - lo_b).upper()
        ea = float(Ea.mid()); eb = float(Eb.mid())
        if Ea <= 0 and Eb <= 0:
            continue
        npos += 1
        if Ea >= 0 and Eb >= 0:
            contrib = A(s)*(Ea + Eb)/2
        else:
            Pp = Ea if Ea > 0 else Eb; Nn = Eb if Ea > 0 else Ea
            contrib = A(s)*Pp*Pp/(2*(Pp - Nn))
        I_up += contrib.upper()
    # ---------- tail ----------
    tail_slope_ok = bool(lam*A(T0 - rho) > arb(7)/2)
    gT0 = xlogx(T0+2) - 4*xlogx(T0) + 3*xlogx(T0-rho) + lam*A(T0) + const - 2*A(T0-rho).log()
    tail_ok = tail_slope_ok and bool(gT0 > ell_a)
    lower = (arb(3)/2 + Lss + ell_a - 3*I_up)
    res = dict(a=a, H=H, N=N, K=K, ell=str(ell), ell_float=float(ell),
               Lss=Lss.str(20, radius=True), Eup=Eup.str(20, radius=True),
               I_up=I_up.str(20, radius=True), pieces_positive=npos,
               gT0=gT0.str(10, radius=True), tail_ok=tail_ok,
               lower_ball=lower.str(20, radius=True),
               # These binary floats are display-only; the directed bounds are
               # carried by the exact interval strings below.
               lower_display=float(lower.lower().mid()), upper_display=float(Eup.upper().mid()),
               lower_certified_str=lower.lower().str(20), upper_certified_str=Eup.upper().str(20),
               seconds=round(time.time()-t0, 1))
    if tau is not None:
        mg = lower.lower() - A(Fr(tau))
        res['tau_cert'] = tau; res['margin_lower_bound'] = mg.lower().str(15); res['margin_positive'] = bool(mg > 0)
    return res

if __name__ == '__main__':
    fin, K, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    only = [int(x) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
    data = json.load(open(fin))
    results = []
    for r in data:
        if only and r['a'] not in only: continue
        H, N = r['H'], r['N']
        n = rationalize(r['dens'], H)
        ell = Fr(int(round(r['ell']*2**30)), 2**30)
        rr = certify(r['a'], H, N, n, ell, K, TAU.get(r['a']))
        rr['sigma_num'] = n; rr['sigma_den'] = DEN; rr['ell_num'] = ell.numerator; rr['ell_den'] = ell.denominator
        print(json.dumps({k: v for k, v in rr.items() if k != 'sigma_num'}), flush=True)
        results.append(rr)
    json.dump(results, open(out, 'w'), indent=1)
