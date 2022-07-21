#ifndef __CURVE_MAKER_H__
#define __CURVE_MAKER_H__

#include "Support/SupportHeader.h"
#include "Shared/MVector.h"

enum CurveMakerFlags
{
	// All samples are constant.
	cmf_Constant = 1 << 0,

	// Allows knots to be duplicated, which will reduce continuity of the curve.
	// Useful for curves that have C0 discontinuities (ie camera cuts).
	cmf_AllowCuts = 1 << 1,
	
	// Allow knots to be removed once we've hit our error constraints.
	// This does a linear scan, so it can be rather slow for curves with many knots.
	cmf_AllowKnotReduction = 1 << 2,
	
	// Force endpoints of curve to match first and last sample.
	cmf_ConstrainEndPos = 1 << 3,
	
	// Force velocity of endpoints to match first and last sample deltas.
	cmf_ConstrainEndVel = 1 << 4,

	// Force velocity of curve end to match velocity of curve start.
	cmf_ForceMatchingEndVel = 1 << 5,

	// Use doubles instead of floats for solver.
	cmf_DoublePrecisionSolver = 1 << 6,

	// If you add flags, increment this number.
	cmf_NumMakerFlags = 7
};
///////////////////////////////////////////////////////////////////////////////
class CurveMaker
{
protected:

	enum SubdivideResult
	{
		sr_Continue = 0,
		sr_Success  = 1,
		sr_Failure  = 2 
	};

	static const float fMinAngErrTol;
	static const float fMinLinErrTol;

	const float * pSamples;
	u32           nSamples;
	u32           nDimensions;
	u32           nDegree;
	u32           nMakerFlags;
	u32           nEvalFlags;
	
	float         fErrTol;
	float         fCutTol;

	MVector<byte> aData;
	u32		      nKnotsSize;    
	u32           nKnotIndicesSize;
	u32           nKnotIndicesOffset;
	u32           nCtrlsSize;
	u32           nCtrlsOffset;
	u32           nCurveSize;
	u32           nCurveOffset; 
	u32           nSolutionMatrixSize;
	u32           nSolutionMatrixOffset;
	u32           nSolutionVectorSize;
	u32           nSolutionVectorOffset;
	u32           nTmpMatrixSize;
	u32           nTmpMatrixOffset;
	u32           nCholeskyDivisorsSize;
	u32           nCholeskyDivisorsOffset;

	float *       pCtrls;
	float *       pTmpCtrls;
	float *       pKnots;
	float *       pTmpKnots;
	s32   *       pKnotIndices;
	u32           nMaxKnots;
	u32           nMaxCurves;

	bool        areAllSamplesEqual() const;
	CurveData * fitCurve(const CurveType * ptype,u32 icur_curve, u32 * pcurve_idx,u32 * pcurve_size);
	CurveData * reduceCurveKnots(const CurveType * ptype,u32 iknot_start,u32 nknots,u32 icur_curve,u32 * pcurve_idx,u32 * pcurve_size);
	u32         subdivideKnots(const CurveType * ptype,u32 iknot_start,u32 nknots,const CurveData & cd,u32 * pcount);
	void        solveControls(const float * pknots,float * pctrls,u32 iknot_start,u32 count);
	byte *      getCurveData(s32 i);
	s32         getNextCurveIndex(s32 i) const;

	// Templated for single/double precision solver.
	template<typename P> P * getSolutionMatrix();
	template<typename P> P * getSolutionVector();
	template<typename P> P * getTmpMatrix();
	template<typename P> P * getCholeskyDivisors();

#if defined(_DEBUG)
	void ASSERT_CONSTRAINTS(const CurveData & cd, const float * pknots, u32 iknot_start, u32 nknots) const;
#else 
	void ASSERT_CONSTRAINTS(const CurveData &, const float * pknots,u32 iknot_start, u32 nknots) const { }
#endif

public:

	u32 getMakerFlags() const { return nMakerFlags; }

	CurveData * makeCompressedCurve(
					const float *       psamples,
					u32                 nsamples,
					u32                 dimensions,
					u32                 degree,
					u32                 maker_flags,
					u32                 eval_flags,
					float               err_tol,
					float               cut_tol,
					bool *              phit_tol = NULL,
					u32 *               psize = NULL,
					const CurveType **  pptypes = NULL, 
					u32                 ntypes = 0          
					);
};
///////////////////////////////////////////////////////////////////////////////
// Pre-multiplies an MxN matrix with it's transpose. 
// Since resulting matrix is symmetrical, only the lower triangle is computed.
// This approach has O(mn^2) complexity.
template<typename P>
void MulTransposeWMat(const P * pa, u32 m, u32 n, P RESTRICT_PTR pc);

// Optimized version which takes advantage of bandwidth, yielding O(m) complexity.
// Note that if A is band-diagonal, then A^t * A is band-symmetric, with band = 2 * band(A) - 1.
template<typename P>
void MulTransposeWMatBandDiagonal(const P * pa, u32 m, u32 n, const s32 * pknot_ids, const s32 * pdof_ids, u32 band, P RESTRICT_PTR pc);

// Applies Cholesky factorization to the lower triangular portion of a symmetric NxN matrix.
// The result is a lower triangular matrix L such that L * L^T = M.
// Optimized to take advantage of band symmetry, yielding O(n) complexity.
template<typename P> 
void CholeskyFactorizeBandSymmetric(const P * pm, u32 n, u32 band, P RESTRICT_PTR pdiv, P RESTRICT_PTR pl);

// Solves L * L^t * x = b, using b as the output memory area.
// Optimized to take advantage of band symmetry, yielding O(n) complexity.
template<typename P>
void CholeskySolveBandSymmetricInPlace(const P * pl, const P * pdiv, u32 n, u32 d, u32 band, P RESTRICT_PTR pb);

// Computes a least-squares-solution of control points, given a set of sample points and knot values.
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
			);

#include "CurveMakerInl.h"

#endif // __CURVE_MAKER_H__