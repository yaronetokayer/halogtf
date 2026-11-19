import numpy as np
import math
from numba import njit, float64, types, void
from halogtf.util.interpolate import interp_linear_to_interfaces
from halogtf.util.calc_linalg import solve_tridiagonal_thomas
from halogtf.util.calc_kp import relative_factors

_TINY64 = np.finfo(np.float64).tiny

@njit(void(float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, float64, types.int64, float64[:], float64[:], float64[:], float64[:], types.int64, float64, float64), cache=True)
def build_tridiag_system(r, m, rho_int, v2, dt, a_param, b_param, c_param, alph, sigma_m_0_char, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref,):
    """
    Build tridiagonal system for implicit conduction update in v2:

        a_i dv2_{i-1} + b_i dv2_i + c_i dv2_{i+1} = d_i

    Requires at least three cells. Supports positive alph and finite or
    infinite w_char. Linearizes the temperature-dependent flux once at
    the current state; timestep retries are not nonlinear Newton iterations.

    Parameters
    ----------
    r : ndarray, shape (N+1,)
        Edge radii.
    m : ndarray, shape (N+1,)
        Enclosed mass at edges.
    rho_int : ndarray, shape (N-1,)
        Interface density values.
    v2 : ndarray, shape (N,)
        Cell-centered one-dimensional velocity dispersion squared.
    dt : float
        Timestep.
    a_param, b_param, c_param : float
        Conductivity parameters.
    alph : float
        Coefficient for interpolation scheme between lmfp and smfp regimes.
        kappa = ( kappa_smfp^-alph + kappa_lmfp^-alph )^(-1/alph)
    sigma_m_0_char : float
        cross section normalization per unit mass in characteristic cross section units.
    w_char : float
        velocity scale w in characteristic velocity units.
    smfp_order : int
        Chapman-Enskog order for the normalized SMFP factor.
        See Outmezguine et al. (2023), Appendix B. Options are 1 and 2.
    a, b, c, d : ndarray, shape (N,)
        Output tridiagonal coefficients and RHS, filled in place.
    model_id : int
        Scattering model integer ID
    k_l_ref, k_s_ref : float
        Model thermal averaged factors for velocity-dependent cross section
    """
    n = v2.shape[0] # Number of cells with unknown dv2 values
    sqrt2 = math.sqrt(2.0)

    smfp_base = (a_param / b_param) * sigma_m_0_char * sigma_m_0_char
    inv_alph = 1.0 / alph

    ### First cell ###
    minusL = 0.0
    plusL = 0.0
    fluxL = 0.0

    rL, rR = r[0], r[1]
    v2C = v2[0]
    mL, mR = 0.0, m[1] 

    # Sliding window
    for j in range(n - 1):
        rRR = r[j + 2]
        v2R = v2[j + 1]
        coefR = rR * rR / (0.5 * (rRR - rL))
        
        dvR = v2R - v2C
        T_face = 0.5 * (v2C + v2R)
        k_lmfp, k_smfp, slope_lmfp, slope_smfp = relative_factors(
            T_face, w_char, model_id, smfp_order, k_l_ref, k_s_ref
        )

        # SMFP and LMFP resistances: k = sqrt(T) / (S^alph+L^alph)^(1/alph).
        S = smfp_base * k_smfp
        L = 1.0 / (c_param * rho_int[j] * T_face * k_lmfp)
        if alph == 1.0:
            D = S + L
            weight_smfp = S / D
            weight_lmfp = L / D
        else:
            # Scale before taking powers to avoid unnecessarily large powers.
            scale = max(S, L)
            scaled_smfp = (S / scale)**alph
            scaled_lmfp = (L / scale)**alph
            total = scaled_smfp + scaled_lmfp
            D = scale * total**inv_alph
            weight_smfp = scaled_smfp / total
            weight_lmfp = scaled_lmfp / total

        k_face = math.sqrt(T_face) / D
        slope_face = (0.5 - weight_smfp * slope_smfp
                      + weight_lmfp * (1.0 + slope_lmfp))

        termR = sqrt2 * k_face
        commonR = termR * dvR * slope_face / (2.0 * T_face)
        minusR = coefR * (commonR - termR)
        plusR = coefR * (commonR + termR)
        fluxR = coefR * termR * dvR

        a[j] = -minusL
        b[j] = minusR - plusL - sqrt2 * (mR - mL) / dt
        c[j] = plusR
        d[j] = fluxL - fluxR

        # Advance sliding window
        minusL, plusL, fluxL = minusR, plusR, fluxR
        rL, rR = rR, rRR
        v2C = v2R
        mL, mR = mR, m[j + 2]

    # Last cell
    a[n - 1] = -minusL
    b[n - 1] = -plusL - sqrt2 * (mR - mL) / dt
    c[n - 1] = 0.0
    d[n - 1] = fluxL

@njit(types.Tuple((float64, float64, types.int64))(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, types.int64, float64, types.int64, float64, float64), cache=True)
def conduct_implicit_nolim(v2, rho, r, m, dv2, rho_int, dt, a_param, b_param, c_param, sigma_m_0, w_char, smfp_order, alph, model_id, k_l_ref, k_s_ref,):
    """
    Implicit conduction step on v2.

    w_char is w / v_s; smfp_order selects the normalized first- or
    second-order SMFP factor. Use w_char = inf for constant scattering.

    The tridiagonal system is defined by:
        a_i dv2_i-1 + b_i dv2_i + c_i dv2_i+1 = d_i
    
    v2 is updated in-place.

    No fractional-change limiter: the caller is responsible for choosing dt.
    Returns (du_max, dt, 0). This wrapper is not used by the default driver.
    """
    N = v2.shape[0]
    du_max = 0.0

    a   = np.empty(N, dtype=np.float64)
    b   = np.empty(N, dtype=np.float64)
    c   = np.empty(N, dtype=np.float64)
    d   = np.empty(N, dtype=np.float64)
    gam = np.empty(N, dtype=np.float64)

    interp_linear_to_interfaces(r, rho, rho_int)

    build_tridiag_system(r, m, rho_int, v2, dt, a_param, b_param, c_param, alph, sigma_m_0, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref,)
    solve_tridiagonal_thomas(a, b, c, d, dv2, gam)

    tiny = _TINY64
    for i in range(N):
        dv2i    = dv2[i]
        v2i     = v2[i]

        denom = v2i if v2i > tiny else tiny
        rat = abs(dv2i) / denom
        if rat > du_max:
            du_max = rat
        v2[i] += dv2i

    return du_max, dt, 0

@njit(types.Tuple((float64, float64, types.int64))(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, types.int64, float64, types.int64, float64, float64), cache=True)
def conduct_implicit_tcool_nolim(v2, rho, r, m, dv2, t_cool, rho_int, dt, a_param, b_param, c_param, sigma_m_0, w_char, smfp_order, alph, model_id, k_l_ref, k_s_ref,):
    """
    Implicit conduction step on v2.

    w_char is w / v_s; smfp_order selects the normalized first- or
    second-order SMFP factor. Use w_char = inf for constant scattering.

    The tridiagonal system is defined by:
        a_i dv2_i-1 + b_i dv2_i + c_i dv2_i+1 = d_i
    
    v2 is updated in-place.
    t_cool is updated in-place.

    No du/u limiter: accepts the supplied dt and returns du_iter = 0.
    """
    N = v2.shape[0]
    du_max = 0.0

    a   = np.empty(N, dtype=np.float64)
    b   = np.empty(N, dtype=np.float64)
    c   = np.empty(N, dtype=np.float64)
    d   = np.empty(N, dtype=np.float64)
    gam = np.empty(N, dtype=np.float64)

    interp_linear_to_interfaces(r, rho, rho_int)

    build_tridiag_system(r, m, rho_int, v2, dt, a_param, b_param, c_param, alph, sigma_m_0, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref)
    solve_tridiagonal_thomas(a, b, c, d, dv2, gam)

    tiny = _TINY64
    for i in range(N):
        v2_old = v2[i]
        dv2i = dv2[i]

        denom = v2_old if v2_old > tiny else tiny

        rat = abs(dv2i) / denom
        if rat > du_max:
            du_max = rat

        if abs(dv2i) > tiny:
            t_cool[i] = denom * dt / abs(dv2i)
        else:
            t_cool[i] = np.inf

        v2[i] = v2_old + dv2i

    return du_max, dt, 0

@njit(types.Tuple((float64, float64, types.int64))(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, types.int64, float64, types.int64, float64, float64), cache=True)
def conduct_implicit_Theta_nolim(v2, rho, r, m, dv2, Th, rho_int, dt, a_param, b_param, c_param, sigma_m_0, w_char, smfp_order, alph, model_id, k_l_ref, k_s_ref,):
    """
    Implicit conduction step on v2.

    w_char is w / v_s; smfp_order selects the normalized first- or
    second-order SMFP factor. Use w_char = inf for constant scattering.

    The tridiagonal system is defined by:
        a_i dv2_i-1 + b_i dv2_i + c_i dv2_i+1 = d_i
    
    v2 is updated in-place. Also updates Theta in place.

    Theta = (v2 / |dv2dt|) / (a_param * sigma_m_0 * dr / sqrt(v2)).
    Both times are in t_s units; the crossing speed is the 1D dispersion.

    No fractional-change limiter: the caller is responsible for choosing dt.
    Returns (du_max, dt, 0). This wrapper is not used by the default driver.
    """
    N = v2.shape[0]
    du_max = 0.0

    a   = np.empty(N, dtype=np.float64)
    b   = np.empty(N, dtype=np.float64)
    c   = np.empty(N, dtype=np.float64)
    d   = np.empty(N, dtype=np.float64)
    gam = np.empty(N, dtype=np.float64)

    interp_linear_to_interfaces(r, rho, rho_int)

    build_tridiag_system(r, m, rho_int, v2, dt, a_param, b_param, c_param, alph, sigma_m_0, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref,)
    solve_tridiagonal_thomas(a, b, c, d, dv2, gam)

    tiny = _TINY64
    rR = r[0]
    for i in range(N):
        dv2i    = dv2[i]
        absdv2i = abs(dv2i)
        v2i     = v2[i]

        rL = rR
        rR = r[i + 1]
        dr = rR - rL

        denom = v2i if v2i > tiny else tiny
        rat = absdv2i / denom
        if rat > du_max:
            du_max = rat

        if absdv2i > tiny and dr > tiny and v2i > tiny:
            Th[i] = v2i * math.sqrt(v2i) * dt / (a_param * sigma_m_0 * dr * absdv2i)
        else:
            Th[i] = np.inf

        v2[i] += dv2i

    return du_max, dt, 0

@njit(types.Tuple((float64, float64, types.int64))(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, types.int64, float64, types.int64, float64, types.int64, float64, float64), cache=True)
def conduct_implicit_dulim(v2, rho, r, m, dv2, rho_int, dt, a_param, b_param, c_param, sigma_m_0, w_char, smfp_order, alph, model_id, eps_du, max_iter, k_l_ref, k_s_ref):
    """
    Implicit conduction step on v2.

    w_char is w / v_s; smfp_order selects the normalized first- or
    second-order SMFP factor. Use w_char = inf for constant scattering.
    Repeatedly solves the implicit system with a trial dt until the
    maximum absolute fractional change satisfies

        max_i |dv2_i| / v2_i <= eps_du

    Then updates v2 in place and returns

        (du_max, dt_used, du_iter)

    The tridiagonal system is defined by:
        a_i dv2_i-1 + b_i dv2_i + c_i dv2_i+1 = d_i
    
    v2 is updated in-place.
    """
    N = v2.shape[0]

    a   = np.empty(N, dtype=np.float64)
    b   = np.empty(N, dtype=np.float64)
    c   = np.empty(N, dtype=np.float64)
    d   = np.empty(N, dtype=np.float64)
    gam = np.empty(N, dtype=np.float64)

    interp_linear_to_interfaces(r, rho, rho_int)

    tiny = _TINY64
    safety = 0.95

    dt_trial = dt
    du_max = math.inf

    for j in range(max_iter):

        build_tridiag_system(r, m, rho_int, v2, dt_trial, a_param, b_param, c_param, alph, sigma_m_0, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref,)
        solve_tridiagonal_thomas(a, b, c, d, dv2, gam)

        du_max = 0.0
        for i in range(N):
            v2i     = v2[i]

            denom = v2i if v2i > tiny else tiny
            rat = abs(dv2[i]) / denom
            # NaN comparisons are false: explicitly reject invalid solutions.
            if not math.isfinite(rat) or v2i + dv2[i] <= 0.0:
                du_max = math.inf
                break

            if rat > du_max:
                du_max = rat
        
        if du_max <= eps_du:
            for i in range(N):
                v2[i] += dv2[i]
            return du_max, dt_trial, j

        # Saturated implicit responses must still make progress on rejection.
        dt_trial *= (min(0.5, safety * eps_du / du_max) if math.isfinite(du_max) else 0.5)

    return du_max, dt_trial, -1

@njit(types.Tuple((float64, float64, types.int64))(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, types.int64, float64, types.int64, float64, types.int64, float64, float64), cache=True)
def conduct_implicit_tcool_dulim(v2, rho, r, m, dv2, t_cool, rho_int, dt, a_param, b_param, c_param, sigma_m_0, w_char, smfp_order, alph, model_id, eps_du, max_iter, k_l_ref, k_s_ref,):
    """
    Implicit conduction step on v2.

    w_char is w / v_s; smfp_order selects the normalized first- or
    second-order SMFP factor. Use w_char = inf for constant scattering.
    Repeatedly solves the implicit system with a trial dt until the
    maximum absolute fractional change satisfies

        max_i |dv2_i| / v2_i <= eps_du

    Then updates v2 in place and returns

        (du_max, dt_used, du_iter)

    The tridiagonal system is defined by:
        a_i dv2_i-1 + b_i dv2_i + c_i dv2_i+1 = d_i
    
    v2 is updated in-place.
    """
    N = v2.shape[0]

    a   = np.empty(N, dtype=np.float64)
    b   = np.empty(N, dtype=np.float64)
    c   = np.empty(N, dtype=np.float64)
    d   = np.empty(N, dtype=np.float64)
    gam = np.empty(N, dtype=np.float64)

    interp_linear_to_interfaces(r, rho, rho_int)

    tiny = _TINY64
    safety = 0.95

    dt_trial = dt
    du_max = math.inf

    for j in range(max_iter):

        build_tridiag_system(r, m, rho_int, v2, dt_trial, a_param, b_param, c_param, alph, sigma_m_0, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref)
        solve_tridiagonal_thomas(a, b, c, d, dv2, gam)

        du_max = 0.0
        for i in range(N):
            v2i     = v2[i]

            denom = v2i if v2i > tiny else tiny
            rat = abs(dv2[i]) / denom
            # NaN comparisons are false: explicitly reject invalid solutions.
            if not math.isfinite(rat) or v2i + dv2[i] <= 0.0:
                du_max = math.inf
                break

            if rat > du_max:
                du_max = rat
        
        if du_max <= eps_du:
            for i in range(N):
                v2_old = v2[i]
                dv2i = dv2[i]

                denom = v2_old if v2_old > tiny else tiny
                if abs(dv2i) > tiny:
                    t_cool[i] = denom * dt_trial / abs(dv2i)
                else:
                    t_cool[i] = np.inf

                v2[i] = v2_old + dv2i
            return du_max, dt_trial, j

        # Saturated implicit responses must still make progress on rejection.
        dt_trial *= (min(0.5, safety * eps_du / du_max) if math.isfinite(du_max) else 0.5)

    return du_max, dt_trial, -1

@njit(types.Tuple((float64, float64, types.int64))(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, types.int64, float64, types.int64, float64, types.int64, float64, float64), cache=True)
def conduct_implicit_Theta_dulim(v2, rho, r, m, dv2, Th, rho_int, dt, a_param, b_param, c_param, sigma_m_0, w_char, smfp_order, alph, model_id, eps_du, max_iter, k_l_ref, k_s_ref,):
    """
    Implicit conduction step on v2.

    w_char is w / v_s; smfp_order selects the normalized first- or
    second-order SMFP factor. Use w_char = inf for constant scattering.
    Repeatedly solves the implicit system with a trial dt until the
    maximum absolute fractional change satisfies

        max_i |dv2_i| / v2_i <= eps_du

    Then updates v2 in place and returns

        (du_max, dt_used, du_iter).

    Also updates Theta in place.

    Theta = (v2 / |dv2dt|) / (a_param * sigma_m_0 * dr / sqrt(v2)).
    Both times are in t_s units; the crossing speed is the 1D dispersion.

    The tridiagonal system is defined by:
        a_i dv2_i-1 + b_i dv2_i + c_i dv2_i+1 = d_i
    
    v2 is updated in-place.
    """
    N = v2.shape[0]

    a   = np.empty(N, dtype=np.float64)
    b   = np.empty(N, dtype=np.float64)
    c   = np.empty(N, dtype=np.float64)
    d   = np.empty(N, dtype=np.float64)
    gam = np.empty(N, dtype=np.float64)

    interp_linear_to_interfaces(r, rho, rho_int)

    tiny = _TINY64
    safety = 0.95

    dt_trial = dt
    du_max = math.inf

    for j in range(max_iter):

        build_tridiag_system(r, m, rho_int, v2, dt_trial, a_param, b_param, c_param, alph, sigma_m_0, w_char, smfp_order, a, b, c, d, model_id, k_l_ref, k_s_ref)
        solve_tridiagonal_thomas(a, b, c, d, dv2, gam)

        du_max = 0.0
        for i in range(N):
            v2i     = v2[i]

            denom = v2i if v2i > tiny else tiny
            rat = abs(dv2[i]) / denom
            # NaN comparisons are false: explicitly reject invalid solutions.
            if not math.isfinite(rat) or v2i + dv2[i] <= 0.0:
                du_max = math.inf
                break

            if rat > du_max:
                du_max = rat
        
        if du_max <= eps_du:
            rR = r[0]
            for i in range(N):
                rL = rR
                rR = r[i + 1]
                dr = rR - rL

                v2i     = v2[i]
                dv2i    = dv2[i]
                absdv2i = abs(dv2i)

                if absdv2i > tiny and dr > tiny and v2i > tiny:
                    Th[i] = v2i * math.sqrt(v2i) * dt_trial / (a_param * sigma_m_0 * dr * absdv2i)
                else:
                    Th[i] = np.inf

                v2[i] += dv2i
            return du_max, dt_trial, j

        # Bound reduction even when the implicit response has saturated.
        fac = (min(0.5, safety * eps_du / du_max) if math.isfinite(du_max) else 0.5)
        dt_trial *= fac

    return du_max, dt_trial, -1
