#ifndef __CURVE_INL_H__
#define __CURVE_INL_H__

#include "Math/PMath.h"
#include "Math/GeneralizedVector.h"

#if 0
#define PROFILE_CURVE_SOLVER
#endif

#if defined(PROFILE_CURVE_SOLVER)
#include "Support/OSWrapper.h"
#include "Shared/PerfCounter.h"
#include "Shared/MString.h"
namespace profile_curve_solver
{
	struct call_timer
	{
		char       name[128];
		EventTimer timer;
		
		void start_call() { timer.startEvent(); }
		void end_call()   { timer.endEvent();   }
		void clear()      { timer.reset();      }

		void flush()
		{
			char buf[512];
			msnprintf(
				buf,sizeof(buf),
				"%s \n"
				"\ttotal secs   = %f\n"
				"\taverage secs = %f\n"
				"\tcalls        = %d\n",
				name, 
				timer.getTotalTime(),
				timer.getEventPeriod(),
				timer.getTotalEvents()
				);
			OSOutputDebugString(buf);
			clear();
		}

		call_timer(const char * pname)
		{
			mstrncpy(name,pname,sizeof(name));
			clear();
		}

		~call_timer() { flush(); }
	};

	struct call_time_tracker
	{
		call_timer & my_call_timer;
		call_time_tracker(call_timer & ct) : my_call_timer(ct) { my_call_timer.start_call(); }
		~call_time_tracker() { my_call_timer.end_call(); }
	};
}
#define PROFILE_CURVE_SOLVER_CALL(name)                 \
static profile_curve_solver::call_timer timer(#name);   \
profile_curve_solver::call_time_tracker tracker(timer);
#endif

#if !defined(PROFILE_CURVE_SOLVER_CALL)
#define PROFILE_CURVE_SOLVER_CALL(name)
#endif
///////////////////////////////////////////////////////////////////////////////
//
// -------------------------------- CurveType ---------------------------------
//
///////////////////////////////////////////////////////////////////////////////
inline const CurveType * CurveType::get(CurveTypeId id)
{
	ASSERT(id >= ct_FirstType && id < ct_MaxTypes);
	ASSERT(aTypes[id]);
	return aTypes[id];
}
///////////////////////////////////////////////////////////////////////////////
inline const CurveType * CurveType::get(const CurveData & cd)
{
	return get(cd.nType);
}
///////////////////////////////////////////////////////////////////////////////
inline void * CurveType::alloc(u32 size)
{
	return malloc(size);
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveType::dealloc(void * pv)
{
	free(pv);
}
///////////////////////////////////////////////////////////////////////////////
//
// -------------------------------- TCurve ------------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename P> const float TCurve<P>::fPrecisionMul = float((1 << (sizeof(P) * 8)) - 1);
template<typename P> const float TCurve<P>::fPrecisionInv = (1.f / float((1 << (sizeof(P) * 8)) - 1));
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
inline const P * TCurve<P>::getBaseData(const byte * pdata,u32 offset)
{
	return reinterpret_cast<const P *>(pdata + offset);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
inline P * TCurve<P>::getBaseData(byte * pdata,u32 offset)
{
	return reinterpret_cast<P *>(pdata + offset);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline const P * TCurve<P>::getKnotData(const T & cd)
{
	return getBaseData(reinterpret_cast<const byte *>(&cd),sizeof(T));
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline const P * TCurve<P>::getCtrlData(const T & cd)
{
	return getKnotData(cd) + cd.nKnots; 
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline P * TCurve<P>::getKnotData(T & cd)
{
	return getBaseData(reinterpret_cast<byte *>(&cd),sizeof(T));
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline P * TCurve<P>::getCtrlData(T & cd)
{
	return getKnotData(cd) + cd.nKnots; 
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline u32 TCurve<P>::getSize(u32 nknots, u32 nctrls)
{
	return sizeof(T) + pw_align<4>(sizeof(P) * (nknots + nctrls));
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline s32 TCurve<P>::getKnotIndex(const T & cd,float & t)
{
	return getKnotIndex(getKnotData(cd),cd.nKnots,t);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline s32 TCurve<P>::getKnotIndex(const P * pk,u32 nknots,float & t)
{
	ASSERT(nknots >= 2);
	const float min_knot = convertToF32(pk[0]);
	const float max_knot = convertToF32(pk[nknots-1]);
	t = Clamp(t,0.f,1.f);
	t = t * max_knot + (1.f - t) * min_knot;
	return getKnotIndexScaled(pk,nknots,min_knot,max_knot,t);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
s32 TCurve<P>::getKnotIndexScaled(const P * pk,u32 nknots,float min_knot,float max_knot,float t)
{
	ASSERT(t >= min_knot && t <= max_knot);

	if(t == min_knot)
	{
		// Return first knot index such that ki+1 > ki.
		return 0;
	}

	if(t == max_knot)
	{
		// Return last knot index such that ki+1 > ki.
		return nknots - 2;
	}

	if(nknots <= 2)
	{
		// Curve only contains a single segment.
		return 0;
	}

	P scaled_t = convertFromF32(t);

	// Knots are in non-decreasing order, so do a modified binary search.
	s32 hi = nknots - 2;
	s32 lo = 0;

	while(hi >= lo)
	{
		int mid = (hi + lo) >> 1;

		if(scaled_t >= pk[mid])
		{
			if(scaled_t < pk[mid + 1])
			{
				// If ki is the knot index for t, then ki <= t < ki+1.
				lo = mid;
				break;
			}
			
			// Loop invariant: scaled_t >= lo.
			ASSERT(scaled_t >= pk[mid]);
			lo = mid + 1;
		}
		else
		{
			// Loop invariant: scaled_t < hi.
			hi = mid;
		}
	}

	ASSERT(hi >= lo);
	return lo;
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline void TCurve<P>::getKnotSupport(const T & cd,s32 ki,float RESTRICT_PTR pres)
{
	extractKnots(cd,ki - cd.nDegree + 1,cd.nDegree << 1,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurve<P>::getKnotSupport(const P * pknots,u32 nknots,u32 degree,s32 ki, float RESTRICT_PTR pres)
{
	ASSERT(ki >= 0 && ki < (s32)(nknots - 1));
	extractKnots(pknots,nknots,degree,ki - degree + 1,degree << 1,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline void TCurve<P>::extractKnots(const T & cd,const s32 base,const u32 count, float RESTRICT_PTR pres)
{
	extractKnots(getKnotData(cd),cd.nKnots,cd.nDegree,base,count,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
void TCurve<P>::extractKnots(const P * pk,u32 nknots,u32 degree,const s32 ibase,const u32 count, float RESTRICT_PTR pres)
{
	if(count == 0)
		return;

	s32 end = ibase + count - 1;
	s32 cur = ibase;

	ASSERT(cur > -s32(degree));
	ASSERT(end < s32(nknots + degree));

	if(cur < 0)
	{
		float k = convertToF32(pk[0]);
		do
		{
			// Duplicate first knot value.
			s32 i = cur - ibase;
			ASSERT(i >= 0 && i < s32(degree << 1));
			pres[i] = k;
			cur++;
		}
		while((cur < 0) && (cur <= end));
	}

	if((end >= (s32)nknots) && (cur <= end))
	{
		float k = convertToF32(pk[nknots - 1]);
		do
		{
			// Duplicate last knot value.
			s32 i = end - ibase;
			ASSERT(i >= 0 && i < s32(degree << 1));
			pres[i] = k;
			end--;
		}
		while((end >= (s32)nknots) && (cur <= end));
	}

	for(; cur <= end; cur++)
	{
		// Copy interior knot values.
		s32 i = cur - ibase;
		ASSERT(i >= 0 && i < s32(degree << 1));
		pres[i] = convertToF32(pk[cur]);
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> template<typename T>
inline void TCurve<P>::scaleBasisDerivatives(const T & cd, u32 k,float RESTRICT_PTR pbasis)
{
	const P * pk = getKnotData(cd);
	const float min_knot = convertToF32(pk[0]);
	const float max_knot = convertToF32(pk[cd.nKnots-1]);
	const float knot_scale = max_knot - min_knot;
	float scale_factor = knot_scale;

	u32 i;
	for(i = 1; i < k; i++)
		scale_factor *= knot_scale;

	const u32 degree = cd.nDegree;
	for(i = 0; i <= degree; i++)
		pbasis[i] *= scale_factor;
}
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
inline float TCurve<P>::convertToF32(P p)
{
	return (float)p;
}
///////////////////////////////////////////////////////////////////////////////	
template<typename P> 
inline P TCurve<P>::convertFromF32(float f)
{
	return (P)f;
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveDataDv --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
inline u32 CurveDataDv::getSize() const
{
	return TCurve<float>::getSize<CurveDataDv>(nKnots,nDimAndCtrl & nCtrlMask);
}
///////////////////////////////////////////////////////////////////////////////
inline s32 CurveDataDv::getKnotIndex(float & t) const
{
	return TCurve<float>::getKnotIndex(*this,t);
}
///////////////////////////////////////////////////////////////////////////////
inline float CurveDataDv::getKnotScale() const
{
	return TCurve<float>::getKnotData(*this)[nKnots - 1];
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataDv::getKnotSupport(s32 ki, float RESTRICT_PTR pres) const
{
	TCurve<float>::getKnotSupport(*this,ki,pres);
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataDv::getCtrlSupport(s32 ki, u32 flags, float RESTRICT_PTR pres) const
{
	// Degree has already been removed from knot index.
	ASSERT(ki >= 0 && ki < nKnots - 1);
	ASSERT((nDimAndCtrl >> nDimShift) * (ki + nDegree) < (nDimAndCtrl & nCtrlMask));

	const float * pc = TCurve<float>::getCtrlData(*this);
	const u32 d = nDimAndCtrl >> nDimShift, k = nDegree, si = ki * d;
	u32 di = 0;

	for(u32 i = 0; i <= k; i++)
	{
		for(u32 n = di + d; di < n; di++)
		{
			ASSERT(di < d * (nDegree + 1));
			ASSERT(si + di < (nDimAndCtrl & nCtrlMask));
			pres[di] = pc[si + di];
		}
	}

	ASSERT((flags & cef_QuatContinuity) == 0 || d == 4);
	if((flags & cef_QuatContinuity) && (d == 4))
		CurveEnsureQuatContinuity(pres,k + 1);
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveDataD3 --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataD3::fixEndian(CurveDataD3 * pcd)
{
	stwbrx((u32 *)pcd->ptBias.v + 0);
	stwbrx((u32 *)pcd->ptBias.v + 1);
	stwbrx((u32 *)pcd->ptBias.v + 2);
	stwbrx((u32 *)pcd->ptScale.v + 0);
	stwbrx((u32 *)pcd->ptScale.v + 1);
	stwbrx((u32 *)pcd->ptScale.v + 2);
	sthbrx(&pcd->nKnots);
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveDataD4 --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataD4::fixEndian(CurveDataD4 * pcd)
{
	for(int i = 0; i < 4; i++)
	{
		stwbrx((u32 *)pcd->aBias + i);
		stwbrx((u32 *)pcd->aScale + i);
	}

	sthbrx(&pcd->nKnots);
}
///////////////////////////////////////////////////////////////////////////////
//
// -------------------------------- TCurveDc ----------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D> template <typename T>
inline u32 TCurveDc<P,D>::getSize(const T & cd)
{
	// This bit of bizarre syntax is required by the C++ standard.
	// Some compilers enforce it (GCC), and some do not (MSVC).
	return TCurve<P>::template getSize<T>(cd.nKnots,getControlCount(cd));
}
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D> template <typename T>
inline u32 TCurveDc<P,D>::getControlCount(const T & cd)
{
	return D * (cd.nKnots + cd.nDegree - 1);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D> template<typename T>
inline void TCurveDc<P,D>::getKnotCtrlSupport(const T & cd,float & t,u32 flags,float RESTRICT_PTR pknots,float RESTRICT_PTR pctrls)
{
	ASSERT(t >= 0.f && t <= 1.f);
	const s32 ki = cd.getKnotIndex(t);
	cd.getKnotSupport(ki,pknots);
	cd.getCtrlSupport(ki,flags,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D>
inline void TCurveDc<P,D>::getKnotCtrlSupport(const TCurveDataD4n<P> & cd,float & t,u32 flags,float RESTRICT_PTR pknots,float RESTRICT_PTR pctrls)
{
	COMPILE_TIME_ASSERT(D == 4);
	ASSERT(t >= 0.f && t <= 1.f);
	const s32 ki = cd.getKnotIndex(t);
	cd.getKnotSupport(ki,pknots);

	float bias[4], scale[4];
	CurveDataD4n::uncompressBiasScale(cd.nQuantizedMinMax,TCurveDataD4n<P>::fNormInv,bias,scale);
	cd.getCtrlSupport(ki,flags,bias,scale,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D> template<typename T>
inline void TCurveDc<P,D>::evaluate(const T & cd,float t,u32 flags,CurveEvalBuf * pres)
{
	const u32 deg = cd.nDegree;
	getKnotCtrlSupport(cd,t,flags,pres->aKnots,pres->aCtrls);
	CurveBasis(t,deg,pres->aKnots,pres->aBasis);
	CurveApplyBasis(pres->aCtrls,pres->aBasis,D,deg,pres->aCurve);

	if(flags & cef_Normalize)
		NormalizeVectorSafeN(pres->aCurve,D);	
}
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D> template<typename T>
inline void TCurveDc<P,D>::differentiate(const T & cd,float t, u32 k, CurveEvalBuf * pres)
{
	const u32 deg = cd.nDegree;
	getKnotCtrlSupport(cd,t,0,pres->aKnots,pres->aCtrls);
	CurveDeriv(t,deg,k,pres->aKnots,pres->aBasis);
	TCurve<P>::scaleBasisDerivatives(cd,k,pres->aBasis);
	CurveApplyBasis(pres->aCtrls,pres->aBasis,D,deg,pres->aCurve);	
}
///////////////////////////////////////////////////////////////////////////////
template<typename P, int D> template<typename T>
void TCurveDc<P,D>::calcError(
						const T &     cd,
						const float * psamples,
						u32           nsamples,
						float         fstart_time,
						float         ftime_scale,
						u32           eval_flags,
						CurveError *  pce,
						const float * ptimes
						)
{
	const float oo_time_scale = (ftime_scale > 0.f) ? 1.f / ftime_scale : 0.f;
	CurveEvalBuf buf;
	
	pce->fTotErrSq = 0.f;
	pce->fMaxErrSq = 0.f;
	pce->iMaxErrSample = -1;
	
	for(u32 i = 0; i < nsamples; i++, psamples += D)
	{
		float t = ptimes ? ptimes[i] : float(i);
		t = (t + fstart_time) * oo_time_scale;
		evaluate(cd,t,eval_flags,&buf);

		float e2 = CurveCalcErrorSq(buf.aCurve,psamples,D,eval_flags);
		pce->fTotErrSq += e2;
		if(e2 > pce->fMaxErrSq)
		{
			pce->fMaxErrSq = e2;
			pce->iMaxErrSample = (s32)i;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------ TCurveDataD3 --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
inline u32 TCurveDataD3<P>::getSize() const
{
	return TCurveDc<P,3>::getSize(*this);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline s32 TCurveDataD3<P>::getKnotIndex(float & t) const
{
	return TCurve<P>::getKnotIndex(*this,t);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline float TCurveDataD3<P>::getKnotScale() const
{
	return TCurve<P>::convertToF32(TCurve<P>::getKnotData(*this)[nKnots - 1]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD3<P>::getKnotSupport(s32 ki, float RESTRICT_PTR pres) const
{
	TCurve<P>::getKnotSupport(*this,ki,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD3<P>::getCtrlSupport(s32 ki,u32 flags, float RESTRICT_PTR pres) const
{
	ASSERT(ki >= 0 && ki < nKnots - 1);
	ASSERT( (3 * (u32)(ki + nDegree) < TCurveDc<P,3>::getControlCount(*this)) );

	const P * pc = TCurve<P>::getCtrlData(*this);
	const s32 k = nDegree * 3;

	for(s32 di = 0, si = ki * 3; di <= k; di += 3, si += 3)
	{
		uncompress(pc + si,ptBias,ptScale,pres + di);
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD3<P>::compress(
								const Point3P & c,
								const Point3P & bias_inv, 
								const Point3P & scale_inv,
								P RESTRICT_PTR  pres
								)
{
	
	ASSERT(c.compareAllGE(bias_inv * -1.f));

	Point3P t;
	t.x = (c.x + bias_inv.x) * scale_inv.x;
	t.y = (c.y + bias_inv.y) * scale_inv.y;
	t.z = (c.z + bias_inv.z) * scale_inv.z;

#if defined(_DEBUG)
	const float m = TCurve<P>::fPrecisionMul + 0.5f;
	ASSERT(t.compareAllLE(makePoint3P(m,m,m)));
	ASSERT(t.compareAllGE(makePoint3P(0.f,0.f,0.f)));
#endif

	// Round to nearest integer.
	const float r = 0.5f;
	pres[0] = TCurve<P>::convertFromF32(t.x + r);
	pres[1] = TCurve<P>::convertFromF32(t.y + r);
	pres[2] = TCurve<P>::convertFromF32(t.z + r);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD3<P>::uncompress(
								const P RESTRICT_PTR pc, 
								const Point3P & bias, 
								const Point3P & scale,
								float RESTRICT_PTR pres
								)
{
	pres[0] = TCurve<P>::convertToF32(pc[0]) * scale.x + bias.x; 
	pres[1] = TCurve<P>::convertToF32(pc[1]) * scale.y + bias.y;
	pres[2] = TCurve<P>::convertToF32(pc[2]) * scale.z + bias.z;

#if defined(_DEBUG)
	float s2 = scale.length2();
	const float m = (s2 != 0.f) ? TCurve<P>::fPrecisionMul + rsqrt(s2) * ZERO_F : 0.f;
	ASSERT(makePoint3P(pres[0],pres[1],pres[2]).compareAllGE(bias));
	ASSERT(makePoint3P(pres[0],pres[1],pres[2]).compareAllLE(bias + scale * m));
#endif
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
CurveDataD3 * TCurveDataD3<P>::create(
								const float * pknots,
								const float * pctrls,
								u32           nknots,
								u32           degree,
								u32           type,
								u32 *         psize,
								byte *        pmem_loc,
								u32           mem_loc_size
								)
{
	const u32 nctrls = (nknots + degree - 1) * 3;
	const u32 nbsize = TCurve<P>::template getSize<CurveDataD3>(nknots,nctrls);

	if(!pmem_loc)
	{
		pmem_loc = (byte *)CurveType::alloc(nbsize);
		mem_loc_size = nbsize;
	}

	if(nbsize > mem_loc_size)
		return NULL;

	CurveDataD3 * pcd = reinterpret_cast<CurveDataD3 *>(pmem_loc);
	pcd->nType   = (u8)type;
	pcd->nDegree = (u8)degree;
	pcd->nKnots  = (u16)nknots;
	
	// Compute AABB for control points.
	Point3P min_pt = makePoint3P(+FLT_MAX,+FLT_MAX,+FLT_MAX);
	Point3P max_pt = makePoint3P(-FLT_MAX,-FLT_MAX,-FLT_MAX);

	for(u32 i = 0; i < nctrls; i += 3)
	{
		min_pt.v[0] = min(min_pt.v[0],pctrls[i+0]);
		min_pt.v[1] = min(min_pt.v[1],pctrls[i+1]);
		min_pt.v[2] = min(min_pt.v[2],pctrls[i+2]);
		
		max_pt.v[0] = max(max_pt.v[0],pctrls[i+0]);
		max_pt.v[1] = max(max_pt.v[1],pctrls[i+1]);
		max_pt.v[2] = max(max_pt.v[2],pctrls[i+2]);
	}

	pcd->ptBias  = min_pt;
	pcd->ptScale = (max_pt - min_pt) * TCurve<P>::fPrecisionInv;

	// Compute inverse bias and scale.
	const Point3P bias_inv  = min_pt * -1.f;
	const Point3P scale_inv = makePoint3P(
		(max_pt.v[0] != min_pt.v[0]) ? TCurve<P>::fPrecisionMul / (max_pt.v[0] - min_pt.v[0]) : 1.f, 
		(max_pt.v[1] != min_pt.v[1]) ? TCurve<P>::fPrecisionMul / (max_pt.v[1] - min_pt.v[1]) : 1.f, 
		(max_pt.v[2] != min_pt.v[2]) ? TCurve<P>::fPrecisionMul / (max_pt.v[2] - min_pt.v[2]) : 1.f);
	
	P RESTRICT_PTR pk = TCurve<P>::getKnotData(*pcd);
	const P * pe = pk + nknots;
	while(pk != pe)
	{
		ASSERT((*pknots >= 0.f) && (floorf(*pknots) == *pknots));
		*pk++ = TCurve<P>::convertFromF32(*pknots++);
	}

	P RESTRICT_PTR pc = TCurve<P>::getCtrlData(*pcd);
	pe = pc + nctrls;
	while(pc != pe)
	{
		TCurveDataD3<P>::compress(makePoint3P(pctrls[0],pctrls[1],pctrls[2]),bias_inv,scale_inv,pc);
		pctrls += 3; pc += 3;
	}

	// Zero out any padding we may have.
	while(u32(pc) < u32(pmem_loc + nbsize))
		*pc++ = 0;

	ASSERT((u32(pc) - u32(pcd)) == nbsize);
	if(psize)
		*psize = nbsize;

	return pcd;
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveDataD4n -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
inline u32 CurveDataD4n::compressMinMax(const float * pmin, const float * pmax)
{
#if defined(_DEBUG)
	s32 i;
	for(i = 0; i < 4; i++)
	{
		ASSERT(pmin[i] <= pmax[i]);
		u32 cmin = (u32)(pmin[0] * 7.5f + 7.5f);
		ASSERT(cmin >= 0 && cmin <= 15);
		u32 cmax = (u32)(pmin[0] * 7.5f + 7.5f);
		ASSERT(cmax >= 0 && cmax <= 15);
	}
#endif

	u32 q;

	// We want quantized min <= min.
	// Default rounding (truncate) should be sufficient.
	q  = (u32)(pmin[0] * 7.5f + 7.5f);
	q |= (u32)(pmin[1] * 7.5f + 7.5f) << 4;
	q |= (u32)(pmin[2] * 7.5f + 7.5f) << 8;
	q |= (u32)(pmin[3] * 7.5f + 7.5f) << 12;

	// Use to ceilf to ensure quantized max >= max.
	q |= (u32)(ceilf(pmax[0] * 7.5f + 7.5f)) << 16;
	q |= (u32)(ceilf(pmax[1] * 7.5f + 7.5f)) << 20;
	q |= (u32)(ceilf(pmax[2] * 7.5f + 7.5f)) << 24;
	q |= (u32)(ceilf(pmax[3] * 7.5f + 7.5f)) << 28;

#if defined(_DEBUG)
	float tmin[4];
	float tmax[4];
	uncompressMinMax(q, tmin, tmax);
	for(i = 0; i < 4; i++)
	{
		ASSERT(tmin[i] <= tmax[i]);
		ASSERT(tmin[i] <= pmin[i]);
		ASSERT(tmax[i] >= pmax[i]);
	}
#endif

	return q;
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataD4n::uncompressMinMax(u32 q, float RESTRICT_PTR pmin, float RESTRICT_PTR pmax)
{
	ASSERT(pmin != pmax);
	pmin[0] = (U4ToF32[(q >>  0) & 0xF] - 7.5f) * (1.f / 7.5f); 
	pmin[1] = (U4ToF32[(q >>  4) & 0xF] - 7.5f) * (1.f / 7.5f); 
	pmin[2] = (U4ToF32[(q >>  8) & 0xF] - 7.5f) * (1.f / 7.5f); 
	pmin[3] = (U4ToF32[(q >> 12) & 0xF] - 7.5f) * (1.f / 7.5f); 
	pmax[0] = (U4ToF32[(q >> 16) & 0xF] - 7.5f) * (1.f / 7.5f); 
	pmax[1] = (U4ToF32[(q >> 20) & 0xF] - 7.5f) * (1.f / 7.5f); 
	pmax[2] = (U4ToF32[(q >> 24) & 0xF] - 7.5f) * (1.f / 7.5f);
	pmax[3] = (U4ToF32[(q >> 28) & 0xF] - 7.5f) * (1.f / 7.5f);
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataD4n::uncompressBiasScale(u32 q, float scale,float RESTRICT_PTR pbias, float RESTRICT_PTR pscale)
{
	ASSERT(pbias != pscale);
	uncompressMinMax(q,pbias,pscale);
	pscale[0] = (pscale[0] - pbias[0]) * scale; 
	pscale[1] = (pscale[1] - pbias[1]) * scale; 
	pscale[2] = (pscale[2] - pbias[2]) * scale; 
	pscale[3] = (pscale[3] - pbias[3]) * scale; 
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveDataD4n::fixEndian(CurveDataD4n * pcd)
{
	stwbrx(&pcd->nQuantizedMinMax);
	sthbrx(&pcd->nKnots);
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------ TCurveDataD4n -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename P> const float TCurveDataD4n<P>::fNormMul   = float((1 << (sizeof(P) * 8 - 1)) - 1);
template<typename P> const float TCurveDataD4n<P>::fNormInv   = 1.f / float((1 << (sizeof(P) * 8 - 1)) - 1);
template<typename P> const u32   TCurveDataD4n<P>::nNormShift = sizeof(P) * 8 - 1;
template<typename P> const u32   TCurveDataD4n<P>::nNormMask  = (1 << (sizeof(P) * 8 - 1)) - 1;
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
inline u32 TCurveDataD4n<P>::getSize() const
{
	return TCurveDc<P,3>::getSize(*this);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline s32 TCurveDataD4n<P>::getKnotIndex(float & t) const
{
	return TCurve<P>::getKnotIndex(*this,t);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline float TCurveDataD4n<P>::getKnotScale() const
{
	return TCurve<P>::convertToF32(TCurve<P>::getKnotData(*this)[nKnots - 1]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4n<P>::getKnotSupport(s32 ki, float RESTRICT_PTR pres) const
{
	TCurve<P>::getKnotSupport(*this,ki,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4n<P>::getCtrlSupport(
								s32                ki, 
								u32                flags,
								const float *      pbias,
								const float *      pscale,
								float RESTRICT_PTR pres
								) const
{
	ASSERT(pbias != pscale);
	ASSERT(pscale != pres);
	ASSERT(ki >= 0 && ki < nKnots - 1);
	ASSERT((3 * (u32)(ki + nDegree) < TCurveDc<P,3>::getControlCount(*this)));

	const P * pc = TCurve<P>::getCtrlData(*this);
	const s32 k = nDegree * 4;

	for(s32 di = 0, si = ki * 3; di <= k; di += 4, si += 3)
	{
		uncompress(pc + si,pbias,pscale,pres + di);
	}

	if(flags & cef_QuatContinuity)
		CurveEnsureQuatContinuity(pres,nDegree + 1);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4n<P>::compressWMinErr(
								u32            flags,
								const float *  pc,
								const float *  pbias,
								const float *  pscale,
								const float *  pbias_inv,
								const float *  pscale_inv,
								P RESTRICT_PTR pres
								)
{
	ASSERT(pc[0] >= pbias[0]);
	ASSERT(pc[1] >= pbias[1]);
	ASSERT(pc[2] >= pbias[2]);
	ASSERT(pc[3] >= pbias[3]);
	ASSERT(pc[0] <= pbias[0] + pscale[0] * fNormMul + ZERO_F);
	ASSERT(pc[1] <= pbias[1] + pscale[1] * fNormMul + ZERO_F);
	ASSERT(pc[2] <= pbias[2] + pscale[2] * fNormMul + ZERO_F);
	ASSERT(pc[3] <= pbias[3] + pscale[3] * fNormMul + ZERO_F);

	// Select implicit coordinate which minimizes error.
	float uncomp_ctrl[4];
	compress(pc,pbias_inv,pscale_inv,0,pres);
	uncompress(pres,pbias,pscale,uncomp_ctrl);

	float min_err_sq = CurveCalcErrorSq(pc,uncomp_ctrl,4,flags);

	for(s32 i = 1; i < 4; i++)
	{
		P comp_ctrl[3];
		compress(pc,pbias_inv,pscale_inv,i,comp_ctrl);
		uncompress(comp_ctrl,pbias,pscale,uncomp_ctrl);

		float err_sq = CurveCalcErrorSq(pc,uncomp_ctrl,4,flags);

		if(err_sq < min_err_sq)
		{
			pres[0] = comp_ctrl[0];
			pres[1] = comp_ctrl[1];
			pres[2] = comp_ctrl[2];
			min_err_sq = err_sq;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4n<P>::compress(
								const float *  pc,
								const float *  pbias_inv,
								const float *  pscale_inv,
								s32            imissing,
								P RESTRICT_PTR pres
								)
{
	ASSERT(IsEqual(Len2VectorN(pc,4),1.f));
	ASSERT(imissing >= 0 && imissing <= 3);
	
	s32 cur = (imissing + 1) & 3;
	const float c0 = (pc[cur] + pbias_inv[cur]) * pscale_inv[cur];
	
	cur = (cur + 1) & 3;
	const float c1 = (pc[cur] + pbias_inv[cur]) * pscale_inv[cur];
	
	cur = (cur + 1) & 3;
	const float c2 = (pc[cur] + pbias_inv[cur]) * pscale_inv[cur];

	// Compress explicit coordinates, rounding to nearest integer.
	const float r = 0.5f;
	pres[0] = TCurve<P>::convertFromF32(c0 + r);
	pres[1] = TCurve<P>::convertFromF32(c1 + r);
	pres[2] = TCurve<P>::convertFromF32(c2 + r);

	// High bit should not be set.
	ASSERT((pres[0] & nNormMask) == pres[0]);
	ASSERT((pres[1] & nNormMask) == pres[1]);
	ASSERT((pres[2] & nNormMask) == pres[2]);

	// Encode sign and index of missing coordinate.
	pres[0] |= (pc[imissing] < 0.f) << nNormShift;
	pres[1] |= (imissing &  1) << nNormShift;
	pres[2] |= (imissing >> 1) << nNormShift;
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4n<P>::uncompress(
								const P *                pc, 
								const float *            pbias,
								const float *            pscale,
								float RESTRICT_PTR       pres
								)
{
	ASSERT(pbias != pscale);
	ASSERT(pscale != pres);
	const s32 imissing = (pc[1] >> nNormShift) | ((pc[2] >> nNormShift) << 1);
	
	s32 cur = (imissing + 1) & 3;
	pres[cur] = TCurve<P>::convertToF32(pc[0] & nNormMask) * pscale[cur] + pbias[cur];
	float sum_sq = square(pres[cur]);
	
	cur = (cur + 1) & 3;
	pres[cur] = TCurve<P>::convertToF32(pc[1] & nNormMask) * pscale[cur] + pbias[cur];
	sum_sq += square(pres[cur]);
	
	cur = (cur + 1) & 3;
	pres[cur] = TCurve<P>::convertToF32(pc[2] & nNormMask) * pscale[cur] + pbias[cur];
	sum_sq += square(pres[cur]);

	float missing_sq = 1.f - Clamp(sum_sq,0.f,1.f);
	if(missing_sq > 0.f)
	{
		// Use rsqrt to save a sqrt call.
		// rsqrt might have a platform-specific implementation, but assume accuracy is similar.
		pres[imissing] = rsqrt(missing_sq) * missing_sq;
		pres[imissing] = (pc[0] >> nNormShift) ? -pres[imissing] : pres[imissing];
	}
	else
	{
		// Sign bit doesn't matter.
		ASSERT(missing_sq == 0.f);
		pres[imissing] = 0.f;
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
CurveDataD4n * TCurveDataD4n<P>::create(
									const float * pknots,
									const float * pctrls,
									u32           nknots,
									u32           degree,
									u32           flags,
									u32           type,
									u32 *         psize,
									byte *        pmem_loc,
									u32           mem_loc_size
									)
{
	const u32 nctrls = (nknots + degree - 1) * 3;
	const u32 nbsize = TCurve<P>::template getSize<CurveDataD4n>(nknots,nctrls);

	if(!pmem_loc)
	{
		pmem_loc = (byte *)CurveType::alloc(nbsize);
		mem_loc_size = nbsize;
	}

	if(nbsize > mem_loc_size)
		return NULL;

	CurveDataD4n * pcd = reinterpret_cast<CurveDataD4n *>(pmem_loc);
	pcd->nType   = (u8)type;
	pcd->nDegree = (u8)degree;
	pcd->nKnots  = (u16)nknots;
	
	float mina[4] = {+FLT_MAX,+FLT_MAX,+FLT_MAX,+FLT_MAX};
	float maxa[4] = {-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

	float norm_ctrl[4];
	for(u32 i = 0, n = nctrls + nknots + degree - 1; i < n; i += 4)
	{
		CopyVectorN(pctrls + i,norm_ctrl,4);
		NormalizeVectorSafeN(norm_ctrl,4);

		mina[0] = min(mina[0],norm_ctrl[0]);
		mina[1] = min(mina[1],norm_ctrl[1]);
		mina[2] = min(mina[2],norm_ctrl[2]);
		mina[3] = min(mina[3],norm_ctrl[3]);
		
		maxa[0] = max(maxa[0],norm_ctrl[0]);
		maxa[1] = max(maxa[1],norm_ctrl[1]);
		maxa[2] = max(maxa[2],norm_ctrl[2]);
		maxa[3] = max(maxa[3],norm_ctrl[3]);
	}

	pcd->nQuantizedMinMax = CurveDataD4n::compressMinMax(mina,maxa);

	// Make sure we use quantized min and max for compression.
	float bias[4], scale[4];
	CurveDataD4n::uncompressBiasScale(pcd->nQuantizedMinMax,fNormInv,bias,scale);

	float bias_inv[4]  = { -bias[0], -bias[1], -bias[2], -bias[3] };
	float scale_inv[4] = 
	{ 
		(scale[0] != 0.f) ? 1.f / scale[0] : 1.f,
		(scale[1] != 0.f) ? 1.f / scale[1] : 1.f,
		(scale[2] != 0.f) ? 1.f / scale[2] : 1.f,
		(scale[3] != 0.f) ? 1.f / scale[3] : 1.f
	};

	P RESTRICT_PTR pk = TCurve<P>::getKnotData(*pcd);
	const P * pe = pk + nknots;
	while(pk != pe)
	{
		ASSERT((*pknots >= 0.f) && (floorf(*pknots) == *pknots));
		*pk++ = TCurve<P>::convertFromF32(*pknots++);
	}

	P RESTRICT_PTR pc = TCurve<P>::getCtrlData(*pcd);
	pe = pc + nctrls;
	while(pc != pe)
	{
		CopyVectorN(pctrls,norm_ctrl,4);
		NormalizeVectorSafeN(norm_ctrl,4);

		TCurveDataD4n<P>::compressWMinErr(flags,norm_ctrl,bias,scale,bias_inv,scale_inv,pc);
		pctrls += 4; // Source controls have 4 coordinates.
		pc += 3;     // Destination controls only have 3 explicit coordinates. 
	}

	// Zero out any padding we may have.
	while(u32(pc) < u32(pmem_loc + nbsize))
		*pc++ = 0;

	ASSERT((u32(pc) - u32(pcd)) == nbsize);
	if(psize)
		*psize = nbsize;

	return pcd;
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------ TCurveDataD4 --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename P> 
inline u32 TCurveDataD4<P>::getSize() const
{
	return TCurveDc<P,4>::getSize(*this);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline s32 TCurveDataD4<P>::getKnotIndex(float & t) const
{
	return TCurve<P>::getKnotIndex(*this,t);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline float TCurveDataD4<P>::getKnotScale() const
{
	return TCurve<P>::convertToF32(TCurve<P>::getKnotData(*this)[nKnots - 1]);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4<P>::getKnotSupport(s32 ki, float RESTRICT_PTR pres) const
{
	TCurve<P>::getKnotSupport(*this,ki,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4<P>::getCtrlSupport(s32 ki,u32 flags, float RESTRICT_PTR pres) const
{
	ASSERT(ki >= 0 && ki < nKnots - 1);
	ASSERT( (4 * (u32)(ki + nDegree) < TCurveDc<P,4>::getControlCount(*this)) );

	const P * pc = TCurve<P>::getCtrlData(*this);
	const s32 k = nDegree * 4;

	for(s32 di = 0, si = ki * 4; di <= k; di += 4, si += 4)
	{
		uncompress(pc + si,aBias,aScale,pres + di);
	}

	if(flags & cef_QuatContinuity)
		CurveEnsureQuatContinuity(pres,nDegree + 1);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4<P>::compress(
								const float * pc,
								const float * pbias_inv,
								const float * pscale_inv,
								P RESTRICT_PTR  pres
								)
{
	float t[4];
	t[0] = (pc[0] + pbias_inv[0]) * pscale_inv[0];
	t[1] = (pc[1] + pbias_inv[1]) * pscale_inv[1];
	t[2] = (pc[2] + pbias_inv[2]) * pscale_inv[2];
	t[3] = (pc[3] + pbias_inv[3]) * pscale_inv[3];

	const float r = 0.5f;
	pres[0] = TCurve<P>::convertFromF32(t[0] + r);
	pres[1] = TCurve<P>::convertFromF32(t[1] + r);
	pres[2] = TCurve<P>::convertFromF32(t[2] + r);
	pres[3] = TCurve<P>::convertFromF32(t[3] + r);
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
inline void TCurveDataD4<P>::uncompress(
								const P * pc, 
								const float * pbias,
								const float * pscale,
								float RESTRICT_PTR pres
								)
{
	pres[0] = TCurve<P>::convertToF32(pc[0]) * pscale[0] + pbias[0]; 
	pres[1] = TCurve<P>::convertToF32(pc[1]) * pscale[1] + pbias[1];
	pres[2] = TCurve<P>::convertToF32(pc[2]) * pscale[2] + pbias[2];
	pres[3] = TCurve<P>::convertToF32(pc[3]) * pscale[3] + pbias[3];
}
///////////////////////////////////////////////////////////////////////////////
template<typename P>
CurveDataD4 * TCurveDataD4<P>::create(
								const float * pknots,
								const float * pctrls,
								u32           nknots,
								u32           degree,
								u32           type,
								u32 *         psize,
								byte *        pmem_loc,
								u32           mem_loc_size
								)
{
	const u32 nctrls = (nknots + degree - 1) * 4;
	const u32 nbsize = TCurve<P>::template getSize<CurveDataD4>(nknots,nctrls);

	if(!pmem_loc)
	{
		pmem_loc = (byte *)CurveType::alloc(nbsize);
		mem_loc_size = nbsize;
	}

	if(nbsize > mem_loc_size)
		return NULL;

	CurveDataD4 * pcd = reinterpret_cast<CurveDataD4 *>(pmem_loc);
	pcd->nType   = (u8)type;
	pcd->nDegree = (u8)degree;
	pcd->nKnots  = (u16)nknots;
	
	// Compute AABB for control points.
	float min_axes[4] = {+FLT_MAX,+FLT_MAX,+FLT_MAX,+FLT_MAX};
	float max_axes[4] = {-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX};

	for(u32 i = 0; i < nctrls; i += 4)
	{
		min_axes[0] = min(min_axes[0],pctrls[i+0]);
		min_axes[1] = min(min_axes[1],pctrls[i+1]);
		min_axes[2] = min(min_axes[2],pctrls[i+2]);
		min_axes[3] = min(min_axes[3],pctrls[i+3]);
		
		max_axes[0] = max(max_axes[0],pctrls[i+0]);
		max_axes[1] = max(max_axes[1],pctrls[i+1]);
		max_axes[2] = max(max_axes[2],pctrls[i+2]);
		max_axes[3] = max(max_axes[3],pctrls[i+3]);
	}

	pcd->aBias[0] = min_axes[0];
	pcd->aBias[1] = min_axes[1];
	pcd->aBias[2] = min_axes[2];
	pcd->aBias[3] = min_axes[3];

	pcd->aScale[0] = (max_axes[0] - min_axes[0]) * TCurve<P>::fPrecisionInv;
	pcd->aScale[1] = (max_axes[1] - min_axes[1]) * TCurve<P>::fPrecisionInv;
	pcd->aScale[2] = (max_axes[2] - min_axes[2]) * TCurve<P>::fPrecisionInv;
	pcd->aScale[3] = (max_axes[3] - min_axes[3]) * TCurve<P>::fPrecisionInv;

	// Compute inverse bias and scale.
	const float bias_inv[4]  = {-min_axes[0],-min_axes[1],-min_axes[2],-min_axes[3]};
	const float scale_inv[4] = {
		(max_axes[0] != min_axes[0]) ? TCurve<P>::fPrecisionMul / (max_axes[0] - min_axes[0]) : 1.f, 
		(max_axes[1] != min_axes[1]) ? TCurve<P>::fPrecisionMul / (max_axes[1] - min_axes[1]) : 1.f, 
		(max_axes[2] != min_axes[2]) ? TCurve<P>::fPrecisionMul / (max_axes[2] - min_axes[2]) : 1.f,
		(max_axes[3] != min_axes[3]) ? TCurve<P>::fPrecisionMul / (max_axes[3] - min_axes[3]) : 1.f};
	
	P RESTRICT_PTR pk = TCurve<P>::getKnotData(*pcd);
	const P * pe = pk + nknots;
	while(pk != pe)
	{
		ASSERT((*pknots >= 0.f) && (floorf(*pknots) == *pknots));
		*pk++ = TCurve<P>::convertFromF32(*pknots++);
	}

	P RESTRICT_PTR pc = TCurve<P>::getCtrlData(*pcd);
	pe = pc + nctrls;
	while(pc != pe)
	{
		TCurveDataD4<P>::compress(pctrls,bias_inv,scale_inv,pc);
		pctrls += 4; pc += 4;
	}

	// Zero out any padding we may have.
	while(u32(pc) < u32(pmem_loc + nbsize))
		*pc++ = 0;

	ASSERT((u32(pc) - u32(pcd)) == nbsize);
	if(psize)
		*psize = nbsize;

	return pcd;
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------ TCurveDataCF32 -----------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<int D>
void TCurveDataCF32<D>::calcError(const float * psamples,u32 nsamples,CurveError * pce) const
{
	pce->fTotErrSq = 0.f;
	pce->fMaxErrSq = 0.f;
	pce->iMaxErrSample = -1;

	for(u32 i = 0; i < nsamples; i++, psamples += D)
	{
		float d2 = CurveCalcErrorSq(psamples,aConstVals,D,0);
		pce->fTotErrSq += d2;
		if(d2 > pce->fMaxErrSq)
		{
			pce->fMaxErrSq = d2;
			pce->iMaxErrSample = (s32)i;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
template<int D>
TCurveDataCF32<D> * TCurveDataCF32<D>::create(                             
											const float * pconst_vals,                                                                                                  
											u32           degree,  
											u32           type,
											u32 *         psize,                                         
											byte *        pmem_loc,
											u32           mem_loc_size
											)
{
	const u32 nbsize = sizeof(TCurveDataCF32<D>);

	if(!pmem_loc)
	{
		pmem_loc = (byte *)CurveType::alloc(nbsize);
		mem_loc_size = nbsize;
	}

	if(nbsize > mem_loc_size)
		return NULL;

	TCurveDataCF32<D> * pcd = reinterpret_cast<TCurveDataCF32<D> *>(pmem_loc);
	pcd->nType   = (u8)type;
	pcd->nDegree = 0;
	pcd->nKnots  = 0;

	for(int i = 0; i < D; i++)
		pcd->aConstVals[i] = pconst_vals[i];

	if(psize)
		*psize = nbsize;

	return pcd;
}
///////////////////////////////////////////////////////////////////////////////
template<int D>
inline void TCurveDataCF32<D>::fixEndian(TCurveDataCF32<D> * pcd)
{
	for(int i = 0; i < D; i++)
		stwbrx((u32 *)pcd->aConstVals + i);
}
///////////////////////////////////////////////////////////////////////////////
//
// -------------------------------- CurveBasis --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P, typename OUT_P>
inline void CurveBasis(IN_P t,u32 d,const IN_P * k,OUT_P RESTRICT_PTR b)
{
	ASSERT(k != (IN_P *)b);

	switch(d)
	{
	case 1:
		CurveBasisLinear(t,k[0],k[1],b,b+1);
		break;
	case 2:
		CurveBasisQuadratic(t,k[0],k[1],k[2],k[3],b,b+1,b+2);
		break;
	case 3:
		CurveBasisCubic(t,k[0],k[1],k[2],k[3],k[4],k[5],b,b+1,b+2,b+3);
		break;
	default:
		CurveBasisGeneral(t,d,k,b);
		break;
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P, typename OUT_P>
inline void CurveBasisLinear(
				IN_P  t, 
				IN_P  ki, 
				IN_P  ki_p1, 
				OUT_P RESTRICT_PTR const bi_m1,
				OUT_P RESTRICT_PTR const bi
				)
{
	ASSERT(ki_p1 > ki);
	ASSERT(bi_m1 != bi);
	const OUT_P r = OUT_P(ki_p1 - t) / OUT_P((ki_p1 - ki));

	*bi_m1 = r;
	*bi    = OUT_P(1.0) - r;
	
	ASSERT(IsEqual(float(*bi_m1 + *bi),1.f));
}
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P, typename OUT_P>
inline void CurveBasisQuadratic(
				IN_P  t,
				IN_P  ki_m1,
				IN_P  ki,
				IN_P  ki_p1,
				IN_P  ki_p2,
				OUT_P RESTRICT_PTR const bi_m2,
				OUT_P RESTRICT_PTR const bi_m1,
				OUT_P RESTRICT_PTR const bi
				)
{
	// + = 1
	// - = 7
	// * = 6
	// / = 2

	ASSERT(bi_m2 != bi_m1);
	ASSERT(bi_m2 != bi);
	ASSERT(bi_m1 != bi);

	const OUT_P d10 = OUT_P(ki_p1 - ki);
	const OUT_P d11 = OUT_P(ki_p1 - ki_m1);
	const OUT_P d20 = OUT_P(ki_p2 - ki);
	const OUT_P d1t = OUT_P(ki_p1 - t);
	const OUT_P dt0 = OUT_P(t - ki);
	const OUT_P dt1 = OUT_P(t - ki_m1);
	const OUT_P d2t = OUT_P(ki_p2 - t);

	ASSERT(d11 * d10 > 0);
	ASSERT(d20 * d10 > 0);

	const OUT_P l = d1t / (d11 * d10);
	const OUT_P r = dt0 / (d20 * d10);

	*bi_m2 = d1t * l;
	*bi_m1 = dt1 * l + d2t * r;
	*bi    = dt0 * r;

	ASSERT(IsEqual(float(*bi_m2 + *bi_m1 + *bi),1.f));
}
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P, typename OUT_P>
inline void CurveBasisCubic(
				IN_P  t,
				IN_P  ki_m2,
				IN_P  ki_m1,
				IN_P  ki,
				IN_P  ki_p1,
				IN_P  ki_p2,
				IN_P  ki_p3,
				OUT_P RESTRICT_PTR const bi_m3,
				OUT_P RESTRICT_PTR const bi_m2,
				OUT_P RESTRICT_PTR const bi_m1,
				OUT_P RESTRICT_PTR const bi
				)
{
	// + = 3
	// - = 12
	// * = 17
	// / = 6

	ASSERT(bi_m3 != bi_m2);
	ASSERT(bi_m3 != bi_m1);
	ASSERT(bi_m3 != bi);
	ASSERT(bi_m2 != bi_m1);
	ASSERT(bi_m2 != bi);
	ASSERT(bi_m1 != bi);

	ASSERT(ki_p1 > ki);
	ASSERT(ki_p1 > ki_m1);
	ASSERT(ki_p1 > ki_m2);
	ASSERT(ki_p2 > ki);
	ASSERT(ki_p2 > ki_m1);
	ASSERT(ki_p3 > ki);

	const OUT_P d10 = OUT_P(1.0) / OUT_P(ki_p1 - ki);
	const OUT_P d11 = OUT_P(1.0) / OUT_P(ki_p1 - ki_m1);
	const OUT_P d12 = OUT_P(1.0) / OUT_P(ki_p1 - ki_m2);
	const OUT_P d20 = OUT_P(1.0) / OUT_P(ki_p2 - ki);
	const OUT_P d21 = OUT_P(1.0) / OUT_P(ki_p2 - ki_m1);
	const OUT_P d30 = OUT_P(1.0) / OUT_P(ki_p3 - ki);

	const OUT_P nt0 = OUT_P(t - ki);
	const OUT_P n1t = OUT_P(ki_p1 - t);
	const OUT_P nt1 = OUT_P(t - ki_m1);
	const OUT_P n2t = OUT_P(ki_p2 - t);
	const OUT_P nt2 = OUT_P(t - ki_m2);
	const OUT_P n3t = OUT_P(ki_p3 - t);

	const OUT_P l0 = n1t * d11 * d10;
	const OUT_P r0 = nt0 * d20 * d10;

	const OUT_P l1 = d12 * n1t * l0;
	const OUT_P l2 = d21 * (nt1 * l0 + n2t * r0);
	const OUT_P r1 = d30 * nt0 * r0;

	*bi_m3 = n1t * l1;
	*bi_m2 = nt2 * l1 + n2t * l2;
	*bi_m1 = nt1 * l2 + n3t * r1;
	*bi    = nt0 * r1;

	ASSERT(IsEqual(float(*bi_m3 + *bi_m2 + *bi_m1 + *bi),1.f));
}
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P, typename OUT_P>
void CurveBasisGeneral(IN_P t,u32 d,const IN_P * k,OUT_P RESTRICT_PTR b)
{
	ASSERT(k != (IN_P *)b);
	ASSERT(d > 0 && d <= MAX_CURVE_DEGREE);
	const s32 c = d - 1;
	b[d] = OUT_P(1.0);

	// Compute basis functions in increasing order, and store in-place.
	for(u32 i = 1; i <= d; i++)
	{
		// Compute Bd-i,i (left-most term)
		ASSERT(k[c+1] > k[c-i+1]);
		b[d-i] = (OUT_P(k[c+1] - t) / OUT_P(k[c+1] - k[c-i+1])) * b[d-i+1];
		
		// Compute [Bd-i+1,i,Bd-1,i] (interior terms)
		for(u32 j = 1; j < i; j++)
		{
			ASSERT(k[c+j] > k[c-i+j]);
			ASSERT(k[c+j+1] > k[c-i+j+1]);
			const OUT_P l = OUT_P(t - k[c-i+j]) / OUT_P(k[c+j] - k[c-i+j]);
			const OUT_P r = OUT_P(k[c+j+1] - t) / OUT_P(k[c+j+1] - k[c-i+j+1]);
			b[d-i+j] = l * b[d-i+j] + r * b[d-i+j+1];
		}

		// Compute Bd,i (right-most term)
		ASSERT(k[c+i] > k[c]);
		b[d] = (OUT_P(t - k[c]) / OUT_P(k[c+i] - k[c])) * b[d];
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P,typename OUT_P>
void CurveBasisGeneralAll(IN_P t,u32 d,const IN_P * k,OUT_P RESTRICT_PTR b)
{
	ASSERT(k != (IN_P *)b);
	ASSERT(d > 0 && d <= MAX_CURVE_DEGREE);
	const s32 c = d - 1;
	b[d] = OUT_P(1.0);

	OUT_P * cb = b + d + 1;
	OUT_P * pb = b;

	// Compute basis functions in increasing order, storing each set in corresponding row.
	for(u32 i = 1; i <= d; i++)
	{
		// Compute Bc-i,i (left-most term)
		ASSERT(k[c+1] > k[c-i+1]);
		cb[d-i] = (OUT_P(k[c+1] - t) / OUT_P(k[c+1] - k[c-i+1])) * pb[d-i+1];
		
		// Compute [Bc-i+1,i,Bc-1,i] (interior terms)
		for(u32 j = 1; j < i; j++)
		{
			ASSERT(k[c+j] > k[c-i+j]);
			ASSERT(k[c+j+1] > k[c-i+j+1]);
			const OUT_P l = OUT_P(t - k[c-i+j]) / OUT_P(k[c+j] - k[c-i+j]);
			const OUT_P r = OUT_P(k[c+j+1] - t) / OUT_P(k[c+j+1] - k[c-i+j+1]);
			cb[d-i+j] = l * pb[d-i+j] + r * pb[d-i+j+1];
		}

		// Compute Bc,i (right-most term)
		ASSERT(k[c+i] > k[c]);
		cb[d] = (OUT_P(t - k[c]) / OUT_P(k[c+i] - k[c])) * pb[d];

		// Advance previous and current basis pointers.
		pb += (d + 1);
		cb += (d + 1); 
	}
}
///////////////////////////////////////////////////////////////////////////////
//
// -------------------------------- CurveDeriv --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename IN_P, typename OUT_P>
void CurveDeriv(IN_P t,u32 d,u32 n,const IN_P * k,OUT_P RESTRICT_PTR b)
{
	ASSERT(d > 0 && d <= MAX_CURVE_DEGREE);
	ASSERT(d >= n);
	const s32 c = d - 1;

	OUT_P buf[(MAX_CURVE_DEGREE+1) * (MAX_CURVE_DEGREE+1) * 2];
	CurveBasisGeneralAll(t,d,k,buf);
	ASSERT(buf[d] == OUT_P(1.0));
	
	OUT_P * prev = buf;
	OUT_P * next = buf + (MAX_CURVE_DEGREE+1) * (MAX_CURVE_DEGREE+1);
	next[d] = OUT_P(1.0);

	for(u32 p = 0; p < n; p++)
	{
		OUT_P * cb = next + d + 1;
		OUT_P * pb = prev;

		for(u32 i = 1; i <= d; i++)
		{
			// Compute dBd-i,i
			ASSERT(k[c+1] > k[c-i+1]);
			cb[d-i] = (OUT_P(-s32(i)) / OUT_P(k[c+1] - k[c-i+1])) * pb[d-i+1];

			// Compute [dBd-i+1,i,dBd-1,i]
			for(u32 j = 1; j < i; j++)
			{
				ASSERT(k[c+j] > k[c-i+j]);
				ASSERT(k[c+j+1] > k[c-i+j+1]);
				const OUT_P l = OUT_P(+i) / OUT_P(k[c+j] - k[c-i+j]);
				const OUT_P r = OUT_P(-s32(i)) / OUT_P(k[c+j+1] - k[c-i+j+1]);
				cb[d-i+j] = l * pb[d-i+j] + r * pb[d-i+j+1];
			}

			// Compute dBd,i
			ASSERT(k[c+i] > k[c]);
			cb[d] = (OUT_P(i) / OUT_P(k[c+i] - k[c])) * pb[d];

			// Advance row pointers for prev and next buffers.
			cb += (d + 1);
			pb += (d + 1);
		}

		mswap(prev,next);
	}

	// Copy terms from last row of previous buffer.
	prev += d * (d + 1);
	for(s32 i = (s32)d; i >= 0; i--)
		b[i] = prev[i];
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveApplyBasis(const float * pctrls,const float * pbasis,u32 dim, u32 degree, float RESTRICT_PTR pres)
{
	u32 i;
	for(i = 0; i < dim; i++, pctrls++)
		pres[i] = (*pctrls) * pbasis[0];

	for(i = 1; i <= degree; i++)
	{
		const float bi = pbasis[i];
		for(u32 j = 0; j < dim; j++, pctrls++)
			pres[j] += ((*pctrls) * bi);
	}
}
///////////////////////////////////////////////////////////////////////////////
inline float CurveCalcErrorSq(const float * pcurve,const float * psample,u32 dim, u32 flags)
{
	if(flags & cef_QuatErrorMetric)
	{
		ASSERT(dim == 4);
		if(DotVectorsN(pcurve,psample,4) < 0.f)
		{
			// Just return sin(Pi/2) ^ 2 if quats are in opposite hemispheres.
			return 1.f;
		}

		Quat qc,qs,qd;
		qc.pt.set(pcurve[0],pcurve[1],pcurve[2],pcurve[3]);
		qc.pt.normalize();
		qs.pt.set(psample[0],psample[1],psample[2],psample[3]);
		qs.pt.normalize();

		// Assume curve = error * sample, so error = curve * sample ^ -1.
		qd = MulQuatQuatInv(qc,qs);

		// Return sin(theta/2) ^ 2.
		return qd.pt.asPoint3().length2();
	}
	else
	{
		// Return distance squared.
		return Dist2VectorsN(pcurve,psample,dim);
	}
}
///////////////////////////////////////////////////////////////////////////////
inline float CurveRadToQuatErrorMetric(float radians)
{
	// See CurveCalcErrorSq for details.
	return sinf(radians * 0.5f);
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveEnsureQuatContinuity(float RESTRICT_PTR pres,u32 count)
{
	float RESTRICT_PTR pc = pres;
	const float * pe = pres + count * 4;	

	float xp = *pc++;
	float yp = *pc++;
	float zp = *pc++;
	float wp = *pc++;

	while(pc < pe)
	{
		const float x = pc[0];
		const float y = pc[1];
		const float z = pc[2];
		const float w = pc[3];

		if((x*xp + y*yp + z*zp + w*wp) < 0.f)
		{
			pc[0] = -x;
			pc[1] = -y;
			pc[2] = -z;
			pc[3] = -w;
		}

		xp = pc[0];
		yp = pc[1];
		zp = pc[2];
		wp = pc[3];

		pc += 4;
	}
}
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveSample ----------------------------------
//
///////////////////////////////////////////////////////////////////////////////
template<typename EMIT_TYPE,int DEG,int DIM,int FLAGS>
void CurveSample(float t,const CurveData ** ppcd,u32 count,byte RESTRICT_PTR pemit,byte RESTRICT_PTR pcache)
{
	COMPILE_TIME_ASSERT(DEG >= 1 && DEG <= MAX_CURVE_DEGREE);
	float curve_buf[DIM]; 
	float basis_buf[DEG + 1];

	for(u32 i = 0; i < count; i++, pcache += sizeof(TCurveSupportCache<DEG,DIM>))
	{
		TCurveSupportCache<DEG,DIM> & cache = *reinterpret_cast<TCurveSupportCache<DEG,DIM> *>(pcache);
		const float scaled_t = t * cache.fKnotScale;
		
		// The required support for k[j] <= scaled_t < k[j+1] is (k[j-DEG+1] ... k[j+DEG]).
		// Cached knots are stored with the assumption that j = DEG-1.
		if((scaled_t < cache.aKnots[DEG-1]) || (scaled_t >= cache.aKnots[DEG]))
		{
			// Cache isn't correct, so update it.
			const CurveType * ptype = CurveType::get(*ppcd[cache.iCurve]);
			ASSERT(ptype);
			ptype->getKnotCtrlSupport(*ppcd[cache.iCurve],t,FLAGS,cache.aKnots,cache.aCtrls);
		}

		CurveBasis(scaled_t,DEG,cache.aKnots,basis_buf);
		CurveApplyBasisT<DEG,DIM>(cache.aCtrls,basis_buf,curve_buf);
		CurveEmit<EMIT_TYPE,DIM,FLAGS & cef_Normalize>(pemit + cache.nOffset,curve_buf);
	}
}
///////////////////////////////////////////////////////////////////////////////
template<typename EMIT_TYPE,int DIM>
void CurveSampleConst(float t,const CurveData ** ppcd,u32 count,byte RESTRICT_PTR pemit,byte RESTRICT_PTR pcache)
{
	for(u32 i = 0; i < count; i++, pcache += sizeof(CurveCache))
	{
		const CurveCache & cache = *reinterpret_cast<const CurveCache *>(pcache);
		const TCurveDataCF32<DIM> & cd = *static_cast<const TCurveDataCF32<DIM> *>(ppcd[cache.iCurve]);
		CurveEmit<EMIT_TYPE,DIM,0>(pemit + cache.nOffset,cd.aConstVals);
	}
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveEmit<QuatXForm,3,0>(byte * pdst,const float * psrc)
{
	QuatXForm * pxf = (QuatXForm *)pdst;
	pxf->transPart.set(psrc[0],psrc[1],psrc[2]);
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveEmit<QuatXForm,4,cef_Normalize>(byte * pdst,const float * psrc)
{
	QuatXForm * pxf = (QuatXForm *)pdst;
	pxf->quatPart.pt.set(psrc[0],psrc[1],psrc[2],psrc[3]);
	pxf->quatPart.normalize();
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveEmit<QuatXForm,4,0>(byte * pdst,const float * psrc)
{
	QuatXForm * pxf = (QuatXForm *)pdst;
	pxf->quatPart.pt.set(psrc[0],psrc[1],psrc[2],psrc[3]);
	ASSERT(IsEqual(pxf->quatPart.length2(),1.f));
}
///////////////////////////////////////////////////////////////////////////////
template<int DEG,int DIM>
inline void CurveApplyBasisT(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	CurveApplyBasis(pctrls,pbasis,DIM,DEG,pres);
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveApplyBasisT<1,3>(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	pres[0]  = (*pctrls++) * pbasis[0];
	pres[1]  = (*pctrls++) * pbasis[0];
	pres[2]  = (*pctrls++) * pbasis[0];

	pres[0] += (*pctrls++) * pbasis[1];
	pres[1] += (*pctrls++) * pbasis[1];
	pres[2] += (*pctrls++) * pbasis[1];
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveApplyBasisT<2,3>(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	pres[0]  = (*pctrls++) * pbasis[0];
	pres[1]  = (*pctrls++) * pbasis[0];
	pres[2]  = (*pctrls++) * pbasis[0];

	pres[0] += (*pctrls++) * pbasis[1];
	pres[1] += (*pctrls++) * pbasis[1];
	pres[2] += (*pctrls++) * pbasis[1];

	pres[0] += (*pctrls++) * pbasis[2];
	pres[1] += (*pctrls++) * pbasis[2];
	pres[2] += (*pctrls++) * pbasis[2];
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveApplyBasisT<3,3>(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	pres[0]  = (*pctrls++) * pbasis[0];
	pres[1]  = (*pctrls++) * pbasis[0];
	pres[2]  = (*pctrls++) * pbasis[0];

	pres[0] += (*pctrls++) * pbasis[1];
	pres[1] += (*pctrls++) * pbasis[1];
	pres[2] += (*pctrls++) * pbasis[1];

	pres[0] += (*pctrls++) * pbasis[2];
	pres[1] += (*pctrls++) * pbasis[2];
	pres[2] += (*pctrls++) * pbasis[2];

	pres[0] += (*pctrls++) * pbasis[3];
	pres[1] += (*pctrls++) * pbasis[3];
	pres[2] += (*pctrls++) * pbasis[3];
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveApplyBasisT<1,4>(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	pres[0]  = (*pctrls++) * pbasis[0];
	pres[1]  = (*pctrls++) * pbasis[0];
	pres[2]  = (*pctrls++) * pbasis[0];
	pres[3]  = (*pctrls++) * pbasis[0];

	pres[0] += (*pctrls++) * pbasis[1];
	pres[1] += (*pctrls++) * pbasis[1];
	pres[2] += (*pctrls++) * pbasis[1];
	pres[3] += (*pctrls++) * pbasis[1];
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveApplyBasisT<2,4>(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	pres[0]  = (*pctrls++) * pbasis[0];
	pres[1]  = (*pctrls++) * pbasis[0];
	pres[2]  = (*pctrls++) * pbasis[0];
	pres[3]  = (*pctrls++) * pbasis[0];

	pres[0] += (*pctrls++) * pbasis[1];
	pres[1] += (*pctrls++) * pbasis[1];
	pres[2] += (*pctrls++) * pbasis[1];
	pres[3] += (*pctrls++) * pbasis[1];

	pres[0] += (*pctrls++) * pbasis[2];
	pres[1] += (*pctrls++) * pbasis[2];
	pres[2] += (*pctrls++) * pbasis[2];
	pres[3] += (*pctrls++) * pbasis[2];
}
///////////////////////////////////////////////////////////////////////////////
template<>
inline void CurveApplyBasisT<3,4>(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres)
{
	pres[0]  = (*pctrls++) * pbasis[0];
	pres[1]  = (*pctrls++) * pbasis[0];
	pres[2]  = (*pctrls++) * pbasis[0];
	pres[3]  = (*pctrls++) * pbasis[0];

	pres[0] += (*pctrls++) * pbasis[1];
	pres[1] += (*pctrls++) * pbasis[1];
	pres[2] += (*pctrls++) * pbasis[1];
	pres[3] += (*pctrls++) * pbasis[1];

	pres[0] += (*pctrls++) * pbasis[2];
	pres[1] += (*pctrls++) * pbasis[2];
	pres[2] += (*pctrls++) * pbasis[2];
	pres[3] += (*pctrls++) * pbasis[2];

	pres[0] += (*pctrls++) * pbasis[3];
	pres[1] += (*pctrls++) * pbasis[3];
	pres[2] += (*pctrls++) * pbasis[3];
	pres[3] += (*pctrls++) * pbasis[3];
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveCache::init(s32 icurve,u32 offset, byte * pb)
{
	CurveCache & cache = *reinterpret_cast<CurveCache *>(pb);
	cache.iCurve  = (u16)icurve;
	cache.nOffset = (u16)offset;
}
///////////////////////////////////////////////////////////////////////////////
inline void CurveSupportCache::init(u32 deg,s32 icurve,u32 emit_offset,float knot_scale, byte * pb)
{
	CurveCache::init(icurve,emit_offset,pb);
	reinterpret_cast<CurveSupportCache *>(pb)->fKnotScale = knot_scale;
	float * pknots = reinterpret_cast<float *>(pb + sizeof(CurveSupportCache));
	for(u32 i = 0, n = deg * 2; i < n; i++)
		pknots[i] = -1.f;
}
///////////////////////////////////////////////////////////////////////////////
inline u32 CurveSupportCache::getSize(u32 deg,u32 dim)
{
	return sizeof(CurveSupportCache) + sizeof(float) * (deg * 2 + deg * dim + dim);
}
///////////////////////////////////////////////////////////////////////////////

#endif // __CURVE_INL_H__
