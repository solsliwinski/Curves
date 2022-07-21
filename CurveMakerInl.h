#ifndef __CURVE_MAKER_INL_H__
#define __CURVE_MAKER_INL_H__

#define BAND_SYMMETRIC_SOLVER
#if 0
#define DEBUG_BAND_SYMMETRIC_SOLVER
#endif

///////////////////////////////////////////////////////////////////////////////
//
// -------------------------------- CurveMaker --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
inline byte * CurveMaker::getCurveData(s32 i)
{
	ASSERT(i >= 0 && i < (s32)nMaxCurves);
	return &aData[nCurveOffset +  i * nCurveSize];
}
///////////////////////////////////////////////////////////////////////////////
inline s32 CurveMaker::getNextCurveIndex(s32 i) const
{
	ASSERT(i >= 0 && i < (s32)nMaxCurves);

	switch(nMaxCurves)
	{
	case 3: 
		// (i + 1) % 3 (Stolen from Christer Ericson and Chris Lentini)
		return (1 << i) & 3; 
	case 2: 
		 // (i + 1) % 2;
		return i ^ 1;   
	default:
		break;
	}

	ASSERT_FAILED;
	return -1;
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline P * CurveMaker::getSolutionMatrix()
{
	return reinterpret_cast<P *>(&aData[nSolutionMatrixOffset]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline P * CurveMaker::getSolutionVector()
{
	return reinterpret_cast<P *>(&aData[nSolutionVectorOffset]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline P * CurveMaker::getTmpMatrix()
{
	return reinterpret_cast<P *>(&aData[nTmpMatrixOffset]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline P * CurveMaker::getCholeskyDivisors()
{
	return reinterpret_cast<P *>(&aData[nCholeskyDivisorsOffset]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
void CurveSolveControls(
			const float *    pknots, 
			u32              nknots, 
			const float *    psamples, 
			u32              nsamples,
			u32              degree, 
			u32              dimensions,
			P RESTRICT_PTR   pmat,
			P RESTRICT_PTR   pvec,
			P RESTRICT_PTR   ptmp_mat,
			P RESTRICT_PTR   pdiv,
			s32 RESTRICT_PTR pknot_indices
			)
{
	PROFILE_CURVE_SOLVER_CALL(CurveSolveControls);
	ASSERT(nknots >= 2);
	ASSERT(nsamples >= 2);

	const u32 dof = degree - 1;

	// Compute dimensions of solution matrix (A).
	const u32 m = nsamples + dof;
	const u32 n = nknots + dof;
	ASSERT(n <= m);
	
	// We're going to be solving the least squares normal equations:
	// A^T * A * x = A^T * b, with x = controls, b = samples.
	// First order of business is to compute A.
	P * pm = pmat;

#if defined(DEBUG_BAND_SYMMETRIC_SOLVER)
	memset(pm,0,sizeof(P) * m * n);
#endif

	const float min_knot = pknots[0];
	const float max_knot = pknots[nknots - 1];
	ASSERT(max_knot - min_knot == float(nsamples - 1));

	// Storage for starting index of each DOF row.
	s32 dof_indices[MAX_CURVE_DEGREE - 1];
	s32 ndof_indices = 0;

	// Compute basis derivatives for first half of DOF.
	float knot_buf[MAX_CURVE_DEGREE * 2];
	u32 i, ki;
	const u32 dof_begin = (dof + 1) >> 1;
	for(i = 0; i < dof_begin; i++)
	{
		// Derivative constraints come in pairs for each sample time.
		float t = min_knot + float(i >> 1);

		ki = TCurve<float>::getKnotIndexScaled(pknots,nknots,min_knot,max_knot,t);	
		TCurve<float>::getKnotSupport(pknots,nknots,degree,ki,knot_buf);

		CurveDeriv(t,degree,1,knot_buf,pm + ki);
		pm += n;

		if(i < dof_begin - 1)
		{
			CurveDeriv(t,degree,2,knot_buf,pm + ki);
			pm += n, i++;

			ASSERT(ndof_indices < MAX_CURVE_DEGREE - 1);
			dof_indices[ndof_indices++] = ki;
		}

		ASSERT(ndof_indices < MAX_CURVE_DEGREE - 1);
		dof_indices[ndof_indices++] = ki;
	}

	// Then, compute basis coefficients for each time sample.
	const u32 ke = nknots - 1;
	ki = (u32)-1;
	for(i = 0; i < nsamples; i++, pm += n)
	{
		float t = min_knot + float(i);
		if((ki == (u32)-1) || (ki >= ke) || (t >= pknots[ki+1]))
		{
			// Previous knot index isn't valid, so search for it.
			ki = TCurve<float>::getKnotIndexScaled(pknots,nknots,min_knot,max_knot,t);
			TCurve<float>::getKnotSupport(pknots,nknots,degree,ki,knot_buf);
		}

		CurveBasis(t,degree,knot_buf,pm + ki);

		// Keep track of starting index for each sample time.
		pknot_indices[i] = ki;
		ASSERT(pknot_indices[i] >= 0);
	}

	// Finally, compute basis derivatives for second half of DOF.
	const u32 dof_end = dof - dof_begin;
	for(i = 0; i < dof_end; i++)
	{
		float t = min_knot + float(nsamples - 1 - (i >> 1));
		ki = TCurve<float>::getKnotIndexScaled(pknots,nknots,min_knot,max_knot,t);		
		TCurve<float>::getKnotSupport(pknots,nknots,degree,ki,knot_buf);

		CurveDeriv(t,degree,1,knot_buf,pm + ki);
		pm += n;

		if(i < dof_end - 1)
		{
			CurveDeriv(t,degree,2,knot_buf,pm + ki);
			pm += n; i++;

			ASSERT(ndof_indices < MAX_CURVE_DEGREE - 1);
			dof_indices[ndof_indices++] = ki;
		}

		ASSERT(ndof_indices < MAX_CURVE_DEGREE - 1);
		dof_indices[ndof_indices++] = ki;
	}

	ASSERT(pm == pmat + m * n);

	// Ok, we have A, now compute A^T * b.
	// Vector dimensionality cheat sheet:
	// b :   MxD (Never explicitly stored)
	// A^Tb: NxD
	// x:    NxD
	P * pv = pvec;
	const P * pve = pv + n * dimensions;
	while(pv != pve)
	{
		// Zero out A^T * b: Controls x Dimensions.
		*pv++ = 0;
	}

	pv = pvec;
	pm = pmat;

	// A^T * b is equivalent to a weighted sum of the bandwidth limited rows of A, 
	// with weights corresponding to the components of each constraint.
	// Since acceleration constraints are zero, we can skip every other row.
	for(i = 0; i < dof_begin; i += 2, pm += n)
	{
		u32 ic = i >> 1;
		u32 in = ic + 1;
		
		// We might not have enough samples to do this properly.
		ic = min(ic,nsamples - 1);
		in = min(in,nsamples - 1);

		const float * pcs = psamples + ic * dimensions;
		const float * pns = psamples + in * dimensions;
		
		const u32 band_start = pknot_indices[ic] * dimensions;
		const P * pband_m = pm + pknot_indices[ic];

		for(u32 j = 0; j < dimensions; j++)
		{
			// Sample time delta is always 1, so no divide is necessary.
			const P ds = P(pns[j] - pcs[j]);
			P RESTRICT_PTR pc = pv + band_start + j;

			for(u32 k = 0; k <= degree; k++, pc += dimensions)
			{
				// V(k,j) = V(k,j) + C(i,j) * A(i,k).
				*pc += (ds * pband_m[k]);
			}
		}

		if(i < dof_begin - 1)
			pm += n;
	}

	// Compute portion of A^T * b corresponding to actual samples.
	const float * ps = psamples;
	for(i = 0; i < nsamples; i++, pm += n, ps += dimensions)
	{
		const s32 band_start = pknot_indices[i] * dimensions;
		const P * pband_m = pm + pknot_indices[i];

		for(u32 j = 0; j < dimensions; j++)
		{
			const P sij = P(ps[j]);
			P RESTRICT_PTR pc = pv + j + band_start;

			for(u32 k = 0; k <= degree; k++, pc += dimensions)
			{
				// V(k,j) = V(k,j) + S(i,j) * A(i,k) 
				*pc += (sij * pband_m[k]);
			}
		}
	}
	
	// Compute final portion of A^T * b corresponding to second half of DOF constraints.
	for(i = 0; i < dof_end; i += 2, pm += n)
	{
		s32 ic = (s32)(nsamples - 1 - (i >> 1));
		s32 ip = ic - 1;
		
		// We might not have enough samples to do this properly.
		ic = max(0,min(ic,s32(nsamples) - 1));
		ip = max(0,min(ip,s32(nsamples) - 1));

		const float * pcs = psamples + ic * dimensions;
		const float * pps = psamples + ip * dimensions;

		const u32 band_start = pknot_indices[ic] * dimensions;
		const P * pband_m = pm + pknot_indices[ic];

		for(u32 j = 0; j < dimensions; j++)
		{
			const P ds = P(pcs[j] - pps[j]);
			P RESTRICT_PTR pc = pv + band_start + j;

			for(u32 k = 0; k <= degree; k++, pc += dimensions)
			{
				// V(k,j) = V(k,j) + C(i,j) * A(i,k).
				*pc += (ds * pband_m[k]);
			}
		}

		if(i < dof_end - 1)
			pm += n;
	}

	ASSERT(pm == pmat + m * n);

#if defined(DEBUG_BAND_SYMMETRIC_SOLVER)
	P * ptmp = (P *)malloc(sizeof(P) * square(max(m,n)) * 4);
	memcpy(ptmp,pmat,sizeof(P) * m * n);
	memcpy(ptmp + m * n + 2 * n * n,pvec,sizeof(P) * n * dimensions);
#endif

	// Increase bandwidth to account for DOF constraints.
	const u32 band = degree + ((degree - 1) >> 4) + 1;

	// Compute lower triangular portion of A^T * A. 
	// After this point, normal equations are ready to be solved.
	MulTransposeWMatBandDiagonal(pmat,m,n,pknot_indices,dof_indices,band,ptmp_mat);

	// Factorize M = A^T * A as L * L^T, where L is lower triangular.
	// Factorization is guaranteed to exist if M is positive definite.
	// A^T * A is positive semi-definite in general, and positive definite in our case.
	CholeskyFactorizeBandSymmetric(ptmp_mat,n,band,pdiv,pmat);
	
	// Factorization produced L, now solve L * L^T * x = A^T * b, in-place.
	CholeskySolveBandSymmetricInPlace(pmat,pdiv,n,dimensions,band,pvec);

#if defined(DEBUG_BAND_SYMMETRIC_SOLVER)
	u32 r;
	P * pata = ptmp + m * n;
	MulTransposeWMat(ptmp,m,n,pata);
	for(r = 0; r < n; r++)
	{
		for(u32 c = r - band + 1; c <= r; c++)
		{
			ASSERT(IsEqual(float(pata[r * n + c]),float(ptmp_mat[r * n + c])));
		}
	}

	P * pl = pata + n * n;
	CholeskyFactorizeBandSymmetric(pata,n,n,pdiv,pl);
	for(r = 0; r < n; r++)
	{
		for(u32 c = r - band + 1; c <= r; c++)
		{
			ASSERT(IsEqual(float(pl[r * n + c]),float(pmat[r * n + c])));
		}
	}

	P * px = pl + n * n;
	CholeskySolveBandSymmetricInPlace(pl,pdiv,n,dimensions,n,px);
	for(r = 0; r < n; r++)
	{
		for(u32 c = 0; c < dimensions; c++)
		{
			ASSERT(IsEqual(float(px[r * dimensions + c]),float(pvec[r * dimensions + c])));
		}
	}

	free(ptmp);
#endif 
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
void MulTransposeWMat(const P * pa, u32 m, u32 n, P RESTRICT_PTR pc)
{
	PROFILE_CURVE_SOLVER_CALL(MulTransposeWMat);

	const P * pai = pa;
	const P * pce = pc + n * n;

	P RESTRICT_PTR pct = pc;
	while(pct != pce)
		*pct++ = 0;
	
	for(u32 i = 0; i < m; i++, pai += n)
	{
		for(u32 j = 0; j < n; j++)
		{
			const P aij = pai[j];
			P RESTRICT_PTR pcj = pc + j * n;

			for(u32 k = 0; k <= j; k++)
			{
				// C(j,k) = C(j,k) + A(i,j) + A(i,k)
				pcj[k] += (aij * pai[k]);
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
void MulTransposeWMatBandDiagonal(
		const P *      pa, 
		u32            m, 
		u32            n, 
		const s32 *    pknot_ids,
		const s32 *    pdof_ids,
		u32            band,
		P RESTRICT_PTR pc
		)
{
	PROFILE_CURVE_SOLVER_CALL(MulTransposeWMatBandDiagonal);
	{
		PROFILE_CURVE_SOLVER_CALL(ZeroATA);

		// Since A is band-diagonal, A^T*A is band-symmetric, with bandwidth = 2 * bandwidth(A) - 1.
		// Therefore, we only have to zero out this portion of A^T*A.
		// This reduces an O(n^2) algorithm to O(n).
		P * prow = pc;
		for(s32 i = 0; i < (s32)n; i++, prow += n)
		{
			const s32 band_start = max(0,i - band + 1);
			const s32 band_end   = min(n,i + band);
			
			P RESTRICT_PTR pband = prow + band_start; 
			const P * pband_end  = prow + band_end;
			while(pband != pband_end)
				*pband++ = 0;
		}
	}

	const s32 dof_begin_count = (s32(band) - 1) >> 1;
	const s32 dof_end_count = s32(band) - 2 - dof_begin_count;

	const P * pai = pa;
	for(s32 i = 0; i < (s32)m; i++, pai += n)
	{
		// Compute start of non-zero band for ith row.
		const s32 bi = (i < dof_begin_count) ? pdof_ids[i]: 
		(
			(i + dof_end_count >= (s32)m) ? 
			pdof_ids[dof_begin_count + s32(m) - i - 1] : 
			pknot_ids[i - dof_begin_count]
		);

		ASSERT(bi >= 0 && bi < s32(n));

		for(s32 j = bi, bi_end = min(s32(n),bi + s32(band)); j < bi_end; j++)
		{
			const P aij = pai[j];
			if(aij == 0)
			{
				// There still might be zeros inside band.
				continue;
			}

			P RESTRICT_PTR pcj = pc + j * n;

			for(s32 k = bi; k <= j; k++)
			{
				// C(j,k) = C(j,k) + A(i,j) + A(i,k)
				pcj[k] += (aij * pai[k]);
			}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
void CholeskyFactorizeBandSymmetric(const P * pm, u32 n, u32 band, P RESTRICT_PTR pdiv, P RESTRICT_PTR pl)
{
	PROFILE_CURVE_SOLVER_CALL(CholeskyFactorizeBandSymmetric);

	ASSERT(pdiv != pl);
	const P * pmi = pm;
	P RESTRICT_PTR pli = pl;

	for(s32 i = 0; i < (s32)n; i++, pmi += n, pli += n)
	{
		const s32 band_start = max(0,i - s32(band) + 1);

		s32 j;
		for(j = band_start; j < i; j++)
		{
			const P * plj = pl + j * n;
			P lij = pmi[j];

			for(s32 k = band_start; k < j; k++)
			{
				// L(i,j) = (M(i,j) - sum(L(i,k) * L(j,k))) / L(j,j)
				lij -= (pli[k] * plj[k]);
			}

			pli[j] = lij * pdiv[j];
		}

		P lii2 = pmi[i];
		for(j = band_start; j < i; j++)
		{
			// L(i,i)^2 = M(i,i) - sum(L(i,j) * L(i,j))
			lii2 -= square(pli[j]);
		}

		if(lii2 > 0)
		{
			// Compute L(i,i) and place reciprocal in divisor storage area.
			pli[i]  = sqrt(lii2);
			pdiv[i] = 1.f / pli[i];
		}
		else
		{
			// This typically means M is under-determined.  Bad, but not fatal.
			// Just pick some reasonable value (zero) in this case.
			pli[i]  = 0;
			pdiv[i] = 0;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
void CholeskySolveBandSymmetricInPlace(const P * pl, const P * pdiv, u32 n, u32 d, u32 band, P RESTRICT_PTR pb)
{
	PROFILE_CURVE_SOLVER_CALL(CholeskySolveBandSymmetricInPlace);

	// M has been factorized into L * L^T, with L lower triangular.
	// Therefore we're trying to solve L * L^T * x = b. 
	// First, solve L * y = b for y.
	ASSERT(n > 0);
	const P * pli = pl;
	P RESTRICT_PTR pyi = pb;

	s32 i;
	for(i = 0; i < (s32)n; i++, pyi += d, pli += n)
	{ 
		const s32 band_start = max(0,i - s32(band) + 1);
		const P * pykj = pb + band_start * d;

		// Reverse order of loops for improved locality.
		for(s32 k = band_start; k < i; k++)
		{
			const P lik = pli[k];
			for(u32 j = 0; j < d; j++, pykj++)
			{
				// Y(i,j) = B(i,j) - sum(L(i,k) * Y(k,j))
				pyi[j] -= (lik * (*pykj));
			}
		}

		for(s32 j = 0; j < (s32)d; j++)
		{
			// Y(i,j) = Y(i,j) / L(i,i)
			pyi[j] *= pdiv[i];
		}
	}

	// Then solve L^T * x = y for x.  
	ASSERT(pyi == pb + d * n);
	P RESTRICT_PTR pxi = pyi - d;

	ASSERT(pli == pl + n * n);
	pli -= n;

	for(i = (s32)(n-1); i >= 0; i--, pxi -= d, pli -= n)
	{
		const s32 band_end = min(n,i + band);
		const P * pxkj = pxi + d;
		const P * plk  = pli + n;

		// Reverse order of loops for improved locality.
		for(s32 k = i+1; k < band_end; k++, plk += n)
		{
			const P lki = plk[i];
			for(u32 j = 0; j < d; j++, pxkj++)
			{
				// X(i,j) = Y(i,j) - sum(L(k,i) * X(k,j))
				pxi[j] -= (lki * (*pxkj));
			}
		}

		for(s32 j = 0; j < (s32)d; j++)
		{
			// X(i,j) = X(i,j) / L(i,i)
			pxi[j] *= pdiv[i];
		}	
	}
}
///////////////////////////////////////////////////////////////////////////////

#endif // __CURVE_MAKER_INL_H__