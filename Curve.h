#ifndef __CURVE_H__
#define __CURVE_H__

#include "Math/PMath.h"
#include "Math/FastMath.h"

struct CurveData;
struct CurveEvalBuf;
template<typename P> struct TCurveDataD4n;

///////////////////////////////////////////////////////////////////////////////
enum CurveTypes
{
	// This format supports controls of any dimension,
	// with 32-bit float precision for knots and controls.
	ct_Dv_F32,
	
	// Each control is three-dimensional, 
	// with knots and control components stored as 16-bit unsigned integers.
	ct_D3_U16,

	// Same as above, with less precision.
	ct_D3_U8,
	
	// This is essentially a compressed quaternion format.
	// For each control, we store 3 explicit 16-bit values.
	// The lower 15 bits are the compressed values.
	// The first high bit indicates the sign of the remaining component.
	// The next two high bits form the index, which is selected to minimize error.
	ct_D4n_U16,

	// Same as above, with less precision.
	ct_D4n_U8,

	// Similar to D3, except controls are 4 dimensional.
	ct_D4_U16,
	ct_D4_U8,

	// Optimized format for a single 3d sample with 32-bit float components.
	ct_D3_ConstF32,

	// Optimized format for a single 4d sample with 32-bit float components.
	ct_D4_ConstF32,

	ct_FirstType = ct_Dv_F32,
	ct_LastType  = ct_D4_ConstF32,
	ct_MaxTypes  = 32
};
///////////////////////////////////////////////////////////////////////////////
enum CurveAttributeFlags
{
	// Hint that curve type will normalize all control points.
	caf_NormalizedCtrls = 1 << 0,

	// Hint that curve type is constant.
	caf_Constant = 1 << 1,

	// Hint that curve type uses aggressive compression.
	caf_HiCompression = 1 << 2,
};
///////////////////////////////////////////////////////////////////////////////
enum CurveEvalFlags
{
	// Indicates that evaluation result should be normalized.
	cef_Normalize = 1 << 0,

	// Indicates that error metric should be non-Euclidean (ie quaternions)
	cef_QuatErrorMetric = 1 << 1,

	// Ensures adjacent control points maintain hemispherical continuity.
	cef_QuatContinuity = 1 << 2,

	// Convenient mask for normal quaternion evaluation.
	cef_Quaternion = cef_QuatErrorMetric | cef_QuatContinuity | cef_Normalize,
};

// Used for evaluating error of a curve over a given set of samples.
struct CurveError
{
	float fMaxErrSq;
	float fTotErrSq;
	s32   iMaxErrSample;
};

///////////////////////////////////////////////////////////////////////////////
typedef u32 CurveTypeId;
typedef void (* CurveSampler)(float t,const CurveData ** ppcd,u32 count,byte RESTRICT_PTR pemit,byte RESTRICT_PTR pcache);

class CurveType : public VirtualBase
{
private:

	// Disable default constructor, copy constructor, and assignment operator.
	CurveType();								     
	CurveType(const CurveType & src);                
	CurveType & operator = (const CurveType & src); 

protected:

	static const CurveType * aTypes[ct_MaxTypes];

public:

	static const CurveType * get(CurveTypeId id);
	static const CurveType * get(const CurveData & cd);
	
	// Get all curve types that can support requested curve specs.
	static u32 getSupportTypes(u32 dim,u32 nknots,CurveType const * * pptypes);

	// Use these methods when allocating memory for curve data.
	static void * alloc(u32 size);
	static void dealloc(void * pv);

	CurveType(CurveTypeId id);
	
	// Returns string name of curve type.
	virtual const char * getName() const = 0;
	
	// Get total size of the curve, including header.
	virtual u32 getSize(const CurveData & cd) const = 0;

	// Get CurveTypeFlags associated with this type.
	virtual u32 getFlags() const = 0;

	// Evaluate curve at time t with 0 <= t <= 1.  
	// Evaluation can be kind of slow, because it needs to fetch knot and curve data.
	virtual void evaluate(const CurveData & cd, float t, u32 flags, CurveEvalBuf * pbuf) const = 0;

	// Similar to evaluate, but computes kth order derivative of the curve.
	// Note that k must be <= degree of the curve.
	virtual void differentiate(const CurveData & cd, float t, u32 k, CurveEvalBuf * pbuf) const = 0;

	// Returns true is this type can support requested specs.
	virtual bool checkSupport(u32 dim,u32 nknots) const = 0;

	// Return appropriate sampling function for given curve data.
	virtual CurveSampler getSampler(const CurveData & cd,u32 flags) const = 0;

	// Get cache size required for sampling this curve.
	virtual u32 getCacheSize(const CurveData & cd) const = 0;

	// Initialize cache for sampling.
	virtual void initCache(const CurveData & cd, s32 icurve,u32 emit_offset,byte * pcache) const = 0;

	// Get knot and control support for normalized time t.
	virtual void getKnotCtrlSupport(
					const CurveData &  cd, 
					float              t,
					u32                flags,
					float RESTRICT_PTR pknots,   
					float RESTRICT_PTR pctrls 
					) const = 0;

	// Create a new curve, possibly applying quantization to knots and controls.  
	// If a memory location of sufficient size is supplied, creation will occur in-place,
	// otherwise memory is allocated.
	virtual CurveData * createCurve(
							const float * pknots,
							const float * pctrls,
							u32           nknots,
							u32           dimensions,
							u32           degree,
							u32           flags,
							u32 *         psize = NULL,
							byte *        pmem_loc = NULL,
							u32           mem_loc_size = 0
							) const = 0;

	// Compute error of the curve, for a given set of samples.
	virtual void calcError(
					const CurveData & cd,
					const float *     psamples,
					u32               nsamples,
					float             start_time,
					float             time_scale,
					u32               flags,
					CurveError *      pce,
					const float *     ptimes = NULL
					) const = 0;

	// Byte swap curve data.
	virtual void fixEndian(CurveData * pcd) const = 0;
};
///////////////////////////////////////////////////////////////////////////////
#define DECLARE_CURVE_TYPE(type,id)                                                                     \
class CurveType##type : public CurveType                                                                \
{                                                                                                       \
public:                                                                                                 \
	static CurveType##type instance;                                                                    \
	CurveType##type(CurveTypeId tid) : CurveType(tid) {}                                                \
	virtual const char * getName() const { return #type; }                                              \
	virtual u32  getSize(const CurveData & cd) const;                                                   \
	virtual u32  getFlags() const;                                                                      \
	virtual void evaluate(const CurveData & cd,float t,u32 flags,CurveEvalBuf * pbuf) const;            \
	virtual void differentiate(const CurveData & cd,float t,u32 k,CurveEvalBuf * pbuf) const;           \
	virtual CurveSampler getSampler(const CurveData & cd, u32 flags) const;                             \
	virtual u32  getCacheSize(const CurveData & cd) const;                                              \
	virtual void initCache(const CurveData & cd, s32 icurve,u32 emit_offset,byte * pcache) const;       \
	virtual void getKnotCtrlSupport(                                                                    \
					const CurveData & cd,                                                               \
					float t,                                                                            \
					u32   flags,		                                                                \
					float RESTRICT_PTR pknots,                                                          \
					float RESTRICT_PTR pctrls                                                           \
					) const;                                                                            \
	virtual bool checkSupport(u32 dim,u32 nknots) const;                                                \
	virtual CurveData * createCurve(                                                                    \
							const float * pknots,                                                       \
							const float * pctrls,                                                       \
							u32           nknots,                                                       \
							u32           dimensions,                                                   \
							u32           degree,                                                       \
							u32           flags,                                                        \
							u32 *         psize = NULL,                                                 \
							byte *        pmem_loc = NULL,                                              \
							u32           mem_loc_size = 0                                              \
							) const;                                                                    \
	virtual void calcError(                                                                             \
					const CurveData & cd,                                                               \
					const float *     psamples,                                                         \
					u32               nsamples,                                                         \
					float             start_time,                                                       \
					float             time_scale,                                                       \
					u32               flags,                                                            \
					CurveError *      pce,                                                              \
					const float *     ptimes = NULL                                                     \
					) const;                                                                            \
 	virtual void fixEndian(CurveData *) const;                                                          \
};                                                                                                      \
CurveType##type CurveType##type::instance(id);

#define DECLARE_BUILTIN_CURVE_TYPE(type) DECLARE_CURVE_TYPE(type,ct_##type)
//////////////////////////////////////////////////////////////////////////////
#define MAX_CURVE_DEGREE 7
#define MAX_CURVE_DIMENSIONS 4
#define MAX_CURVE_HEADER_SIZE 512

struct CurveEvalBuf
{
	float aKnots[MAX_CURVE_DEGREE * 2];
	float aCtrls[MAX_CURVE_DIMENSIONS * (MAX_CURVE_DEGREE + 1)];
	float aBasis[MAX_CURVE_DEGREE + 1];
	float aCurve[MAX_CURVE_DIMENSIONS];
};
///////////////////////////////////////////////////////////////////////////////
struct CurveData
{
	// Should be one of the CurveType enums.
	u8 nType;

	// Degree of spline equation (typically 3 or less).
	u8 nDegree;

	// A closed, non-uniform knot vector looks like [t0,t1 .. tk,tk+1 ... tn, tn+1,tn+2 .. tn+k+1],
	// with t0 = t1 = .. tk, and tn+1 = tn+2 = ... tn+k+1.
	// As a space optimization, only the first instance of each knot is stored.
	// Therefore nKnots = n - k + 2.
	u16 nKnots;
};
///////////////////////////////////////////////////////////////////////////////
template<typename P>
class TCurve
{
	// Disable default contructor.
	TCurve() { NYI; }

public:

	// Useful if we don't have an explicit type.
	static const P * getBaseData(const byte * pdata,u32 offset);
	static P * getBaseData(byte * pdata,u32 offset);

	// Get starting address of knot data.
	template<typename T> static const P * getKnotData(const T & cd);
	template<typename T> static P * getKnotData(T & cd);

	// Get starting address of control data.
	template<typename T> static const P * getCtrlData(const T & cd);
	template<typename T> static P * getCtrlData(T & cd);

	// Get size of curve, including header.
	template<typename T> static u32 getSize(u32 nknots,u32 nctrls);

	// Find ki such that ki <= t < ki+1. t is assumed to be in [0,1].
	template<typename T> static s32 getKnotIndex(const T & cd,float & t);
	static s32 getKnotIndex(const P * pk,u32 nknots,float & t);

	// Similar to above, but assumes t has already been scaled and clamped.
	static s32 getKnotIndexScaled(const P * pk,u32 nknots,float min_knot,float max_knot,float t);
	
	// Get buffer of knot values needed to evaluate basis functions at ki.
	// The support ranges from ki - degree + 1 to ki + degree.
	template<typename T> static void getKnotSupport(const T & cd, s32 ki, float RESTRICT_PTR pres);
	static void getKnotSupport(const P * pknots,u32 nknots,u32 degree,s32 ki,float RESTRICT_PTR pres);
	
	// Extract knot values from base to base + count - 1.
	template<typename T> static void extractKnots(const T & cd,const s32 base,const u32 count, float RESTRICT_PTR pres);
	static void extractKnots(const P * pknots,u32 nknots,u32 degree,const s32 base,const u32 count, float RESTRICT_PTR pres);

	// Internally, knot values can, and typically do, range from M to M + N-1.
	// Externally, curves must be sampled with normalized time.
	// Therefore, the differential has a scale factor that needs to be removed,
	// if the resulting basis function is going to make sense with normalized time.
	// In general, the kth basis derivative needs to remove (knot_scale ^ k).
	template<typename T> static void scaleBasisDerivatives(const T & cd, u32 k,float RESTRICT_PTR pbasis);

	// Allow customization of type conversion.
	static float convertToF32(P k);
	static P convertFromF32(float k);

	// Scalars for compressing/decompressing values to a certain precision.
	static const float fPrecisionMul;
	static const float fPrecisionInv;
};
///////////////////////////////////////////////////////////////////////////////
struct CurveDataDv : public CurveData
{ 
	// Control count is stored in lower 24 bits, dimensions in upper 8 bits.
	// Number of dimensions also implied by (nKnots - 1 + nDegree) / nControls.
	u32 nDimAndCtrl;

	enum
	{
		nDimShift = 24,
		nCtrlMask = 0xFFFFFF
	};

	// Knot and control values are always 32-bit floats for this curve type.
	u32   getSize() const;
	s32   getKnotIndex(float & t) const;
	float getKnotScale() const;
	void  getKnotSupport(s32 ki, float RESTRICT_PTR pres) const;
	void  getCtrlSupport(s32 ki, u32 flags, float RESTRICT_PTR pres) const;
	void  getKnotCtrlSupport(float & t,u32 flags,float RESTRICT_PTR pknots,float RESTRICT_PTR pctrls) const;
	void  evaluate(float t,u32 flags,CurveEvalBuf * pres) const;
	void  differentiate(float t, u32 k, CurveEvalBuf * pres) const;
	void  calcError(
			const float * psamples,
			u32           nsamples,
			float         fstart_time,
			float         ftime_scale,
			u32           eval_flags,
			CurveError *  perr,
			const float * ptimes = NULL
			) const;

	static CurveDataDv * create(
					const float * pknots,
					const float * pctrls,
					u32           nknots,
					u32           dimensions,
					u32           degree,
					u32 *         psize,
					byte *        pmem_loc,
					u32           mem_loc_size
					);

	static void fixEndian(CurveDataDv * pcd);
};
///////////////////////////////////////////////////////////////////////////////
template<int D>
struct TCurveDataCF32 : public CurveData
{
	float aConstVals[D];

	void calcError(const float * psamples,u32 nsamples,CurveError * pce) const;

	static TCurveDataCF32 * create(                             
								const float * pconst_vals,                                                                                                 
								u32           degree, 
								u32           type,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								);

	static void fixEndian(TCurveDataCF32 * pcd);
};

typedef TCurveDataCF32<3> CurveDataC3F32;
typedef TCurveDataCF32<4> CurveDataC4F32;
///////////////////////////////////////////////////////////////////////////////
template<typename P,int D>
class TCurveDc 
{
public:

	template<typename T> static u32  getSize(const T & cd);
	template<typename T> static u32  getControlCount(const T & cd);
	template<typename T> static void getKnotCtrlSupport(const T & cd,float & t,u32 flags,float RESTRICT_PTR pknots,float RESTRICT_PTR pctrls);
	template<typename T> static void evaluate(const T & cd,float t,u32 flags,CurveEvalBuf * pres);
	template<typename T> static void differentiate(const T & cd,float t, u32 k, CurveEvalBuf * pres);

	template<typename T> static void calcError(
			const T &     cd,
			const float * psamples,
			u32           nsamples,
			float         fstart_time,
			float         ftime_scale,
			u32           eval_flags,
			CurveError *  perr,
			const float * ptimes = NULL
			);

	// Overload D4n types since they extract controls differently.
	static void getKnotCtrlSupport(const TCurveDataD4n<P> & cd,float & t,u32 flags,float RESTRICT_PTR pknots,float RESTRICT_PTR pctrls);
};
///////////////////////////////////////////////////////////////////////////////
struct CurveDataD3 : public CurveData
{
	// Scale and bias for each coordinate axis.
	Point3P ptScale;
	Point3P ptBias;

	static void fixEndian(CurveDataD3 * pc);
};

template<typename P>
struct TCurveDataD3 : public CurveDataD3
{
	u32   getSize() const;
	s32   getKnotIndex(float & t) const;
	float getKnotScale() const;
	void  getKnotSupport(s32 ki, float RESTRICT_PTR pres) const;
	void  getCtrlSupport(s32 ki, u32 flags, float RESTRICT_PTR pres) const;

	// static so an instance isn't required to compress/uncompress control values.
	static void compress(
			const Point3P & c,
			const Point3P & bias_inv, 
			const Point3P & scale_inv,
			P RESTRICT_PTR pres
			);

	static void uncompress(
			const P * pc, 
			const Point3P & bias, 
			const Point3P & scale, 
			float RESTRICT_PTR pres
			);

	static CurveDataD3 * create(
			const float * pknots,
			const float * pctrls,
			u32           nknots,
			u32           degree,
			u32           type,
			u32 *         psize,
			byte *        pmem_loc,
			u32           mem_loc_size
			);
};

typedef TCurveDataD3<u16> CurveDataD3U16;
typedef TCurveDataD3<u8>  CurveDataD3U8;
///////////////////////////////////////////////////////////////////////////////
struct CurveDataD4n : public CurveData
{
	// Encodes 4-bit min and max values for each axis.
	u32 nQuantizedMinMax;

	static u32  compressMinMax(const float * pmin, const float * pmax);
	static void uncompressMinMax(u32 qmin_max,float RESTRICT_PTR pmin,float RESTRICT_PTR pmax);
	static void uncompressBiasScale(u32 qmin_max,float scale,float RESTRICT_PTR pbias, float RESTRICT_PTR pscale);
	static void fixEndian(CurveDataD4n * pc);

	static const float U4ToF32[16];
};

template<typename P>
struct TCurveDataD4n : public CurveDataD4n
{
	u32   getSize() const;
	s32   getKnotIndex(float & t) const;
	float getKnotScale() const;
	void  getKnotSupport(s32 ki, float RESTRICT_PTR pres) const;
	void  getCtrlSupport(
			s32                ki, 
			u32                flags,
			const float *      pbias,
			const float *      pscale,
			float RESTRICT_PTR pres
			) const;

	// Compress once for each coordinate, and select version which minimizes error.
	static void compressWMinErr(
			u32            flags,
			const float *  pc,
			const float *  pbias,
			const float *  pscale,
			const float *  pbias_inv,
			const float *  pscale_inv,
			P RESTRICT_PTR pres
			);

	// Compress, given an implicit coordinate. 
	static void compress(
			const float *  pc,
			const float *  pbias_inv,
			const float *  pscale_inv,
			s32            implicit,
			P RESTRICT_PTR pres
			);

	static void uncompress(
			const P *           pc, 
			const float *       pbias,
			const float *       pscale,
			float RESTRICT_PTR  pres
			);

	static CurveDataD4n * create(
			const float * pknots,
			const float * pctrls,
			u32           nknots,
			u32           degree,
			u32           flags,
			u32           type,
			u32 *         psize,
			byte *        pmem_loc,
			u32           mem_loc_size
			);

	static const float fNormMul;
	static const float fNormInv;
	static const u32   nNormShift;
	static const u32   nNormMask;
};

typedef TCurveDataD4n<u16> CurveDataD4nU16;
typedef TCurveDataD4n<u8>  CurveDataD4nU8;
///////////////////////////////////////////////////////////////////////////////
struct CurveDataD4 : public CurveData
{
	// Scale and bias for each coordinate axis.
	float aBias[4];
	float aScale[4];

	static void fixEndian(CurveDataD4 * pc);
};

template<typename P>
struct TCurveDataD4 : public CurveDataD4
{
	u32   getSize() const;
	s32   getKnotIndex(float & t) const;
	float getKnotScale() const;
	void  getKnotSupport(s32 ki, float RESTRICT_PTR pres) const;
	void  getCtrlSupport(s32 ki, u32 flags, float RESTRICT_PTR pres) const;

	// static so an instance isn't required to compress/uncompress control values.
	static void compress(
			const float * pc,
			const float * pbias_inv,
			const float * pscale_inv,
			P RESTRICT_PTR pres
			);

	static void uncompress(
			const P * pc, 
			const float * pbias, 
			const float * pscale, 
			float RESTRICT_PTR pres
			);

	static CurveDataD4 * create(
			const float * pknots,
			const float * pctrls,
			u32           nknots,
			u32           degree,
			u32           type,
			u32 *         psize,
			byte *        pmem_loc,
			u32           mem_loc_size
			);
};

typedef TCurveDataD4<u16> CurveDataD4U16;
typedef TCurveDataD4<u8>  CurveDataD4U8;
///////////////////////////////////////////////////////////////////////////////
// Compute supporting BSpline basis functions for a given sample time, with:
// b[0] ... b[d] -> [Bi-d,d ... Bi,d] and 
// Bi,d = (t - k[i]) / (k[i+d] - k[i]) * Bi,d-1 + (k[i+d+1] - t) / (k[i+d+1] - ki[i+1]) * Bi+1,d-1. 
template<typename IN_P, typename OUT_P>
void CurveBasis(IN_P t,u32 degree,const IN_P * knots,OUT_P RESTRICT_PTR basis);

template<typename IN_P, typename OUT_P>
void CurveBasisLinear(
		IN_P  t, 
		IN_P  ki, 
		IN_P  ki_p1, 
		OUT_P RESTRICT_PTR const bi_m1,
		OUT_P RESTRICT_PTR const bi
		);

template<typename IN_P,typename OUT_P>
void CurveBasisQuadratic(
		IN_P  t,
		IN_P  ki_m1,
		IN_P  ki,
		IN_P  ki_p1,
		IN_P  ki_p2,
		OUT_P RESTRICT_PTR const bi_m2,
		OUT_P RESTRICT_PTR const bi_m1,
		OUT_P RESTRICT_PTR const bi
		);

template<typename IN_P,typename OUT_P>
void CurveBasisCubic(
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
		);

template<typename IN_P,typename OUT_P>
void CurveBasisGeneral(IN_P t,u32 d,const IN_P * knots,OUT_P RESTRICT_PTR basis);

// Similar to preceding function, except we return the basis for 0 <= i <= d.
template<typename IN_P,typename OUT_P>
void CurveBasisGeneralAll(IN_P t,u32 d,const IN_P * knots,OUT_P RESTRICT_PTR basis);

// Derivative of basis function is given as:
// (dBi,d)/dt = d / (k[i+d] - ki[i]) * Bi,d-1 - d / (k[i+d+1] - k[i+1]) * Bi+1,d-1.
// Nth order derivative follows immediately:
// Bi,d,n = d / (k[i+d] - k[i]) * Bi,d-1,n-1 - d / (k[i+d+1] - k[i+1]) * Bi+1,d-1,n-1
template<typename IN_P, typename OUT_P>
void CurveDeriv(IN_P t,u32 d,u32 n,const IN_P * pk,OUT_P RESTRICT_PTR pb);

///////////////////////////////////////////////////////////////////////////////
// Some helper functions:
float CurveCalcErrorSq(const float * pcurve,const float * psample,u32 dim, u32 flags);
float CurveRadToQuatErrorMetric(float radians);
void  CurveEnsureQuatContinuity(float RESTRICT_PTR pres,u32 count);
void  CurveApplyBasis(const float * pctrls,const float * pbasis,u32 dim, u32 degree, float RESTRICT_PTR pres);

// CurveSample operates on an array of curves and cache data.  This is much faster than calling evaluate.
template<typename EMIT_TYPE,int DEG,int DIM,int FLAGS>
void CurveSample(float t,const CurveData ** ppcd,byte RESTRICT_PTR pemit, byte RESTRICT_PTR pcache);

template<typename EMIT_TYPE,int DIM>
void CurveSampleConst(float t,const CurveData ** ppcd,byte RESTRICT_PTR pemit, byte RESTRICT_PTR pcache);

template<typename EMIT_TYPE,int DIM,int FLAGS>
void CurveEmit(byte * pdst,const float * psrc);

template<int DEG,int DIM>
void CurveApplyBasisT(const float * pctrls,const float * pbasis,float RESTRICT_PTR pres);
///////////////////////////////////////////////////////////////////////////////
struct CurveSamplerCache
{
	CurveSampler pSampler;
	u16          nEntries;
	u16          nTotSize;
};
///////////////////////////////////////////////////////////////////////////////
struct CurveCache
{
	u16 iCurve;
	u16 nOffset;

	static void init(s32 icurve,u32 offset, byte * pb);
};
///////////////////////////////////////////////////////////////////////////////
struct CurveSupportCache : public CurveCache
{
	float fKnotScale;

	static void init(u32 deg,s32 icurve, u32 offset,float knot_scale, byte * pb);
	static u32  getSize(u32 deg,u32 dim);
};
///////////////////////////////////////////////////////////////////////////////
template<int DEG,int DIM>
struct TCurveSupportCache : public CurveSupportCache
{
	float aKnots[DEG * 2];
	float aCtrls[DEG * DIM + DIM];
};
///////////////////////////////////////////////////////////////////////////////

#include "CurveInl.h"

#endif // __CURVE_H__