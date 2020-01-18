// Evaluation of the Mittag-Leffler (ML) function with 1, 2 or 3 parameters
// by means of the OPC algorithm [1]. The routine evaluates an approximation
// Et of the ML function E such that |E-Et|/(1+|E|) approx 1.0e-15  
//    
// Adapted directly from Roberto Garrappa's MATLAB program, 
// mathworks.com/matlabcentral/fileexchange/48154-the-mittag-leffler-function
//
// TO USE THIS IN STATA, RUN THIS FILE ONCE AS IF IT WERE A DO FILE
// THIS WILL MAKE AVAILABLE THE MATA FUNCTION mlf()
//
// E = ML(z,alpha) evaluates the ML function with one parameter alpha for
// the corresponding elements of z alpha must be a real and positive
// scalar. The one parameter ML function is defined as
//
// E = sum_{k=0}^{infty} z^k/Gamma(alpha*k+1)
//
// with Gamma the Euler's gamma function.
//
//
// E = ML(z,alpha,beta) evaluates the ML function with two parameters alpha
// and beta for the corresponding elements of z alpha must be a real and
// positive scalar and beta a real scalar. The two parameters ML function is
// defined as
//
// E = sum_{k=0}^{infty} z^k/Gamma(alpha*k+beta)
//
//
// E = ML(z,alpha,beta,gamma) evaluates the ML function with three parameters
// alpha, beta and gamma for the corresponding elements of z alpha must be a
// real scalar such that 0<alpha<1, beta any real scalar and gamma a real and
// positive scalar the arguments z must satisfy |Arg(z)| > alpha*pi. The
// three parameters ML function is defined as
//
// E = sum_{k=0}^{infty} Gamma(gamma+k)*z^k/Gamma(gamma)/k!/Gamma(alpha*k+beta)
//
//
// NOTE:
// This routine implements the optimal parabolic contour (OPC) algorithm
// described in [1] and based on the inversion of the Laplace transform on a
// parabolic contour suitably choosen in one of the regions of analyticity
// of the Laplace transform.
//
//
// REFERENCES
//
//   [1] R. Garrappa, Numerical evaluation of two and three parameter
//   Mittag-Leffler functions, SIAM Journal of Numerical Analysis, 2015,
//   53(3), 1350-1369
//
//
//   Please, report any problem or comment to :
//        roberto dot garrappa at uniba dot it
//
//   Copyright (c) 2015, Roberto Garrappa, University of Bari, Italy
//   roberto dot garrappa at uniba dot it
//   Homepage: http://www.dm.uniba.it/Members/garrappa
//   Revision: 1.4 - Date: October 8 2015
	
mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

// main routine. optional args beta and gamma default to 1

numeric colvector mlf(numeric colvector z, real scalar alpha, | real scalar beta, real scalar gamma) {
	numeric colvector E; real scalar k, log_epsilon

	if (alpha <= 0 | gamma <= 0)
		return (.)

	// Check inputs
	if (args() < 4) {
		gamma = 1
		if (args() < 3)
			beta = 1
	}

	// Check parameters and arguments for the three parameter case
	if (abs(gamma-1) > epsilon(1)) {
		if (alpha > 1)
			_error(198, "With the three parameters Mittag-Leffler function the parameter ALPHA must satisfy 0 < alpha < 1.")

		if (isreal(z))
			_error(198, "With the three parameters Mittag-Leffler function this code works only when |Arg(z)| > alpha*pi.")
		if (any(abs(arg(select(z, abs(z)>epsilon(1)))) :> alpha*pi()))
			_error(198, "With the three parameters Mittag-Leffler function this code works only when |Arg(z)| > alpha*pi.")
	}

	// Target precision
	log_epsilon = ln(1e-15) 

	// Inversion of the LT for each element of z
	E = J(rows(z), 1, isreal(z)? 1/gamma(beta) : C(1/gamma(beta)))  
	for (k=rows(z); k; k--)
		if (abs(z[k]) >= 1.0e-15)
			E[k] = LTInversion(1, z[k], alpha, beta, gamma, log_epsilon)

	return(E)
}


// =========================================================================
// Version for special case: t=1, x<0, beta=0, gamma=1, 0<alpha<1, but x can be a colvector
// =========================================================================
real colvector _mlf(real colvector x, real scalar alpha) {
	real scalar N, mu, h, w, log_eps, log_epsilon; real colvector u; complex colvector z, zd; complex matrix F

	log_epsilon = ln(1e-15) 
	log_eps = ln(epsilon(1))
	mu = log_epsilon - log_eps                   // log ratio of desired to max precision
	w = sqrt(log_eps / (log_eps - log_epsilon))  // half-width of integration range, evidently needed to assure given precision
	N = ceil(-w * log_epsilon / (2*pi()))        // half the number of integration points
	h = w / N                                    // width of bars in Riemann integral
	u = rangen(-w, w, N+N+1)'                    // integration points
	z  = C(1, u); z = mu*z:*z                    // z(u) = mu * (1+ui)^2 (mu controls how close it comes to the origin)
	zd = (-mu) * C(u, -1)                        // dz/du
	F = (zd :* exp(z)) :/ (1 :- x * z:^-alpha)   // integrand: dz/du * exp(z(u)) / (1 - x * z(u))^alpha -- first line that depends on the inputs!
	return(h / pi() * Im(quadrowsum(F)))         // integral divided by 2*pi*i
}


// =========================================================================
// Evaluation of the ML function by Laplace transform inversion
// =========================================================================
numeric scalar LTInversion(real scalar t, numeric scalar lambda, real scalar alpha, real scalar beta, real scalar gamma, real scalar log_epsilon) {
	real scalar theta, kmin, kmax, tau, j1, J1, i, N, iN, mu, h, _log_epsilon, ln10
	numeric colvector s_star, z, zd, zexp, F, S
	numeric scalar Integral, Residues, ss_star
	real colvector k_vett, phi_s_star, index_s_star, index_save, p, admissible_regions, k, u
	real rowvector v
	pragma unset N

	tau = 2 * pi(); _log_epsilon = log_epsilon

	// Evaluation of the relevant poles
	theta = isreal(lambda)? (lambda<0? pi() : 0) : arg(lambda)
	kmin = ceil(-alpha/2 - theta/tau)
	kmax = floor(alpha/2 - theta/tau)
	k_vett = kmin>kmax? J(0,1,0) : kmin :: kmax
	s_star = abs(lambda)^(1/alpha) * exp(1i / alpha * (theta :+ k_vett * tau))

	// Evaluation of phi(s_star) for each pole
	phi_s_star = (Re(s_star) + abs(s_star)) * 0.5

	// Sorting of the poles according to the value of phi(s_star)
	index_s_star = order(phi_s_star, 1)
	      s_star =     s_star[index_s_star]
	  phi_s_star = phi_s_star[index_s_star]

	// Deleting possible poles with phi_s_star=0
	index_save = phi_s_star :> 1.0e-15
	    s_star = select(    s_star, index_save)
	phi_s_star = select(phi_s_star, index_save)

	// Inserting the origin in the set of the singularities
	if (rows(s_star)) {
		s_star     = 0 \     s_star
		phi_s_star = 0 \ phi_s_star
		J1 = rows(s_star)
	} else {
		s_star = phi_s_star = 0
		J1 = 1
	}

	// Strength of the singularities
	p = max((0,-2*(alpha*gamma-beta+1))) \ J(J1-1, 1, gamma)
	phi_s_star = phi_s_star \ .

	// Looking for the admissible regions with respect to round-off errors
	admissible_regions = selectindex(phi_s_star[|.\J1|] :< (log_epsilon - ln(epsilon(1))) / t :& phi_s_star[|.\J1|] :< phi_s_star[|2\.|])

	// Evaluation of parameters for inversion of LT in each admissible region
	ln10 = ln(10)
	for (; N > 200; _log_epsilon = _log_epsilon + ln10) {  // N starts at . as default value on creation
		N = .
		for (i=rows(admissible_regions);i;i--) {
			j1 = admissible_regions[i]
			v = j1 < J1? OptimalParam_RB(t, phi_s_star[j1], phi_s_star[j1+1], p[j1], gamma, _log_epsilon) :
									 OptimalParam_RU(t, phi_s_star[j1],                   p[j1],        _log_epsilon)

			// Selection of the admissible region for integration which involves the
			// minimum number of nodes
			if (v[3] < N) {
				mu = v[1]
				h  = v[2]
				N  = v[3]
				iN = j1
			}
		}
	}

	// Evaluation of the inverse Laplace transform
	k = -N :: N
	u = h*k
	z  = C(1, u); z = mu*z:*z
	zd = (-2*mu) * C(u, -1)
	zexp = exp(z*t)
	F = z :^ (alpha * gamma - beta) :/ (z:^alpha :- lambda) :^ gamma :* zd
	S = zexp :* F
	Integral = h / tau * quadcolsum(S) / 1i

	// Evaluation of residues
	if (iN < J1) {
		ss_star = s_star[|iN+1\.|]
		Residues = quadcolsum(ss_star:^(1-beta) :* exp(t * ss_star)) / alpha
	} else
		Residues = 0

	return(isreal(lambda)? Re(Integral + Residues) : Integral + Residues)
}


// =========================================================================
// Finding optimal parameters in a right-bounded region
// returns 3-rowvector with mu, h, N
// =========================================================================
real rowvector OptimalParam_RB (real scalar t, real scalar phi_s_star_j, real scalar phi_s_star_j1, real scalar pj, real scalar qj, real scalar log_epsilon) {
	real scalar log_eps,  _log_epsilon, fac, conservative_error_analysis, f_max, f_min, f_bar, fp, fq, sq_phi_star_j, threshold, sq_phi_star_j1, sq_phibar_star_j, sq_phibar_star_j1, phibar_star_j1, adm_region, w, den, muj, hj, Nj

	// Definition of some constants
	log_eps = -1.205966f2b4f12X+005  // ln(epsilon(1))
	fac = 1.01
	conservative_error_analysis = 0

	// Maximum value of fbar as the ration between tolerance and round-off unit
	f_max = exp(log_epsilon - log_eps)

	// Evaluation of the starting values for sq_phi_star_j and sq_phi_star_j1
	sq_phi_star_j = sqrt(phi_s_star_j)
	threshold = 2 * sqrt((log_epsilon - log_eps) / t)
	sq_phi_star_j1 = min((sqrt(phi_s_star_j1), threshold - sq_phi_star_j))

	// Zero or negative values of pj and qj
	if (pj < 1e-14 & qj < 1e-14) {
		sq_phibar_star_j  = sq_phi_star_j
		sq_phibar_star_j1 = sq_phi_star_j1
		adm_region = 1
	}

	// Zero or negative values of just pj
	if (pj < 1e-14 & qj >= 1e-14) {
		sq_phibar_star_j = sq_phi_star_j
		f_min = sq_phi_star_j <= 0? fac : fac * (sq_phi_star_j / (sq_phi_star_j1 - sq_phi_star_j)) ^ qj
		if (adm_region = f_min < f_max) {
			f_bar = f_min * (2 - f_min/f_max)
			fq = f_bar ^ (-1/qj)
			sq_phibar_star_j1 = (2 * sq_phi_star_j1 - fq * sq_phi_star_j) / (2 + fq)
		}
	}

	// Zero or negative values of just qj
	if (pj >= 1e-14 & qj < 1e-14) {
		sq_phibar_star_j1 = sq_phi_star_j1
		f_min = fac * (sq_phi_star_j1 / (sq_phi_star_j1 - sq_phi_star_j)) ^ pj
		if (adm_region = f_min < f_max) {
			f_bar = f_min * (2 - f_min/f_max)
			fp = f_bar ^ (-1/pj)
			sq_phibar_star_j = (2 * sq_phi_star_j + fp * sq_phi_star_j1)/(2 - fp)
		}
	}

	// Positive values of both pj and qj
	if (pj >= 1e-14 & qj >= 1e-14) {
		f_min = fac*(sq_phi_star_j + sq_phi_star_j1) / (sq_phi_star_j1 - sq_phi_star_j) ^ max((pj,qj))
		if (adm_region = f_min < f_max) {
			f_min = max((f_min,1.5))
			f_bar = f_min * (2 - f_min/f_max)
			fp = f_bar ^ (-1/pj)
			fq = f_bar ^ (-1/qj)
			w = conservative_error_analysis? -2 * phi_s_star_j1 * t / (log_epsilon - phi_s_star_j1 * t) : -phi_s_star_j1 * t / log_epsilon
			den = 2+w - (1+w)*fp + fq
			sq_phibar_star_j  = ( (2+w+fq)*sq_phi_star_j +            fp *sq_phi_star_j1)/den
			sq_phibar_star_j1 = (-(1+w)*fq*sq_phi_star_j + (2+w-(1+w)*fp)*sq_phi_star_j1)/den
		}
	}

	if (adm_region) {
		phibar_star_j1 = sq_phibar_star_j1 * sq_phibar_star_j1
		_log_epsilon = log_epsilon - ln(f_bar)  // in Mata, unlike R and MATLAB, functions can modify arguments
		w = conservative_error_analysis? -2 * phibar_star_j1 * t / (_log_epsilon - phibar_star_j1 * t) : -phibar_star_j1 * t / _log_epsilon
		muj = (((1+w)*sq_phibar_star_j + sq_phibar_star_j1)/(2+w))^2
		hj = -2*pi() / _log_epsilon * (sq_phibar_star_j1 - sq_phibar_star_j) / ((1 + w) * sq_phibar_star_j + sq_phibar_star_j1)
		Nj = ceil(sqrt(1 - _log_epsilon / t / muj) / hj)

		return((muj, hj, Nj))
	}
	return((0,0,.))
}

// =========================================================================
// Finding optimal parameters in a right-unbounded region
// returns 3-rowvector with mu, h, N
// =========================================================================
real rowvector OptimalParam_RU (real scalar t, real scalar phi_s_star_j, real scalar pj, real scalar log_epsilon) {
	real scalar phibar_star_j, sq_phi_s_star_j, sq_phibar_star_j, sq_muj, A, f_min, f_max, f_tar, fbar, phi_t, log_eps_phi_t, muj, hj, Nj, w, u, Q, log_eps, threshold

	// Evaluation of the starting values for sq_phi_star_j
	phibar_star_j = phi_s_star_j > 0? phi_s_star_j*1.01 : 0.01
	 sq_phi_s_star_j = sqrt(phi_s_star_j)
	sq_phibar_star_j = sqrt(phibar_star_j)

	// Definition of some constants
	f_min = 1; f_max = 10; f_tar = 5

	// Iterative process to look for fbar in [f_min,f_max]
	while (1) {
		phi_t = phibar_star_j*t
		log_eps_phi_t = log_epsilon/phi_t
		Nj = ceil(phi_t/pi()*(1 - 1.5 * log_eps_phi_t + sqrt(1-2*log_eps_phi_t)))
		A = pi() * Nj / phi_t
		sq_muj = sq_phibar_star_j * abs((4-A)/(7-sqrt(1+12*A)))
		fbar = ((sq_phibar_star_j - sq_phi_s_star_j) / sq_muj) ^ -pj

		if (pj < 1e-14 | (f_min < fbar & fbar < f_max))
			break

		sq_phibar_star_j = f_tar^(-1/pj)*sq_muj + sq_phi_s_star_j
		phibar_star_j = sq_phibar_star_j * sq_phibar_star_j
	}
	muj = sq_muj * sq_muj
	hj = (-3*A - 2 + 2*sqrt(1+12*A)) / (4-A) / Nj

	// Adjusting integration parameters to keep round-off errors under control
	log_eps = ln(epsilon(1))
	threshold = (log_epsilon - log_eps)/t
	if (muj > threshold) {
		Q = abs(pj) < 1e-14? 0 : f_tar^(-1/pj) * sqrt(muj)
		phibar_star_j = Q + sqrt(phi_s_star_j); phibar_star_j = phibar_star_j*phibar_star_j
		if (phibar_star_j >= threshold)
			return((muj, 0, .))
		w = sqrt(log_eps/(log_eps-log_epsilon))
		u = sqrt(-phibar_star_j*t/log_eps)
		muj = threshold
		Nj = ceil(w*log_epsilon/2/pi()/(u*w-1))
		hj = w / Nj
	}
	return((muj, hj, Nj))
}


// validate against cases generated by Tran Quoc Viet, https://github.com/tranqv/Mittag-Leffler-function-and-its-derivative/blob/master/tcases/README.md
// single argument is the path to the "tcases" directory of the Tran test files
void verify_mlf(string scalar directory) {
	real scalar i, N, alpha, beta, r; transmorphic scalar fh; complex scalar z, myE; real rowvector line

	chdir(directory)
	for (i=80;i;i--) {
		i
		fh = _fopen("tt_mlfm_c"+strofreal(i,"%02.0f")+".txt", "r")
		if (fh >= 0) {
			line = strtoreal(tokens(fget(fh)))
			N = line[2]; alpha = line[3]; beta = line[4]
			for (r=N;r;r--) {
				line = strtoreal(tokens(fget(fh)))
				z = C(line[1], line[2])
				myE = mlf(z, alpha, beta)
				if (reldif(myE, C(line[3], line[4])) > 1.75e-12)
					printf("Failure alpha = %21x, beta = %21x, z = %21x + %21x i, E = %21x + %21x i, my E = %21x + %21x i\n", alpha, beta, line[1], line[2], line[3], line[4], Re(myE), Im(myE))
			}
			fclose(fh)
		}
	}
}

mata mlib create lMittagLeffler, dir("`c(sysdir_plus)'l") replace
mata mlib add lMittagLeffler *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
