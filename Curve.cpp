#include "Support/SupportHeader.h"
#include "Curve.h"

const CurveType * CurveType::aTypes[ct_MaxTypes];

///////////////////////////////////////////////////////////////////////////////
CurveType::CurveType(CurveTypeId id)
{
	ASSERT_PRINTF(!aTypes[id], "Multiple curve types declared with id = %d.\n",id);
	aTypes[id] = this;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveType::getSupportTypes(u32 dim,u32 nknots,CurveType const * * pptypes)
{
	u32 n = 0;
	for(u32 i = ct_FirstType; i <= ct_LastType; i++)
	{
		ASSERT(aTypes[i]);
		if(aTypes[i]->checkSupport(dim,nknots))
			pptypes[n++] = aTypes[i];
	}
	
	return n;
}
///////////////////////////////////////////////////////////////////////////////
enum CurveSamplerTypes
{
	cst_Point3   = 0,
	cst_Quat     = 1,
	cst_NumTypes = 2, 
};
///////////////////////////////////////////////////////////////////////////////
static CurveSampler CurveSamplerTable[MAX_CURVE_DEGREE + 1][cst_NumTypes] = 
{
	&CurveSampleConst<QuatXForm,3>,&CurveSampleConst<QuatXForm,4>,
	&CurveSample<QuatXForm,1,3,0>,&CurveSample<QuatXForm,1,4,cef_Normalize>,
	&CurveSample<QuatXForm,2,3,0>,&CurveSample<QuatXForm,2,4,cef_Normalize>,
	&CurveSample<QuatXForm,3,3,0>,&CurveSample<QuatXForm,3,4,cef_Normalize>,
	&CurveSample<QuatXForm,4,3,0>,&CurveSample<QuatXForm,4,4,cef_Normalize>,
	&CurveSample<QuatXForm,5,3,0>,&CurveSample<QuatXForm,5,4,cef_Normalize>,
	&CurveSample<QuatXForm,6,3,0>,&CurveSample<QuatXForm,6,4,cef_Normalize>,
	&CurveSample<QuatXForm,7,3,0>,&CurveSample<QuatXForm,7,4,cef_Normalize>
};
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveTypeDv_F32 ------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(Dv_F32);

u32 CurveTypeDv_F32::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	COMPILE_TIME_ASSERT(sizeof(CurveDataDv) <= MAX_CURVE_HEADER_SIZE);
	return static_cast<const CurveDataDv &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeDv_F32::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	static_cast<const CurveDataDv &>(cd).evaluate(t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeDv_F32::differentiate(const CurveData & cd, float t,u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	static_cast<const CurveDataDv &>(cd).differentiate(t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeDv_F32::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	const u32 dim = (static_cast<const CurveDataDv &>(cd).nDimAndCtrl >> CurveDataDv::nDimShift);
	switch(dim)
	{
	case 3: return (flags == 0) ? CurveSamplerTable[cd.nDegree][cst_Point3] : NULL;
	case 4: return (flags == cef_Quaternion) ? CurveSamplerTable[cd.nDegree][cst_Quat] : NULL;
	}

	ASSERT_FAILED;
	return NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeDv_F32::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	const u32 dim = (static_cast<const CurveDataDv &>(cd).nDimAndCtrl >> CurveDataDv::nDimShift);
	return CurveSupportCache::getSize(cd.nDegree,dim);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeDv_F32::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataDv &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeDv_F32::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_Dv_F32);
	static_cast<const CurveDataDv &>(cd).getKnotCtrlSupport(t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeDv_F32::checkSupport(u32 dim,u32 nknots) const
{
	return (nknots <= 0xFFFF);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeDv_F32::getFlags() const
{
	return 0;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeDv_F32::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree, 
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT(degree <= MAX_CURVE_DEGREE && dimensions <= MAX_CURVE_DIMENSIONS);
	return CurveDataDv::create(pknots,pctrls,nknots,dimensions,degree,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeDv_F32::calcError(
						const CurveData & cd,
						const float *     psamples,
						u32               nsamples,
						float             start_time,
						float             time_scale,
						u32               flags,
						CurveError *      pce,
						const float *     ptimes
						) const 
{
	ASSERT(cd.nType == ct_Dv_F32);
	static_cast<const CurveDataDv &>(cd).calcError(psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeDv_F32::fixEndian(CurveData * pc) const
{
	ASSERT(pc->nType == ct_Dv_F32);
	CurveDataDv::fixEndian(static_cast<CurveDataDv *>(pc));
}
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveTypeD3_U16 ------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D3_U16);

u32 CurveTypeD3_U16::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D3_U16);
	COMPILE_TIME_ASSERT(sizeof(CurveDataD3U16) <= MAX_CURVE_HEADER_SIZE);
	return static_cast<const CurveDataD3U16 &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U16::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D3_U16);
	TCurveDc<u16,3>::evaluate(static_cast<const CurveDataD3U16 &>(cd),t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U16::differentiate(const CurveData & cd, float t,u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D3_U16);
	TCurveDc<u16,3>::differentiate(static_cast<const CurveDataD3U16 &>(cd),t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD3_U16::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D3_U16);
	return (flags == 0) ? CurveSamplerTable[cd.nDegree][cst_Point3] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD3_U16::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D3_U16);
	return CurveSupportCache::getSize(cd.nDegree,3);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U16::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D3_U16);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataD3U16 &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U16::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_D3_U16);
	TCurveDc<u16,3>::getKnotCtrlSupport(static_cast<const CurveDataD3U16 &>(cd),t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD3_U16::checkSupport(u32 dim,u32 nknots) const
{
	return (nknots <= 0xFFFF) && (dim == 3) ;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD3_U16::getFlags() const
{
	return 0;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD3_U16::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree, 
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 3));
	return CurveDataD3U16::create(pknots,pctrls,nknots,degree,ct_D3_U16,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U16::calcError(
						const CurveData & cd, 
						const float *     psamples, 
						u32               nsamples, 
						float             start_time, 
						float             time_scale,
						u32               flags,
						CurveError *      pce,
						const float *     ptimes
						) const 
{
	ASSERT(cd.nType == ct_D3_U16);
	TCurveDc<u16,3>::calcError(static_cast<const CurveDataD3U16 &>(cd),psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U16::fixEndian(CurveData * pc) const
{
	ASSERT(pc->nType == ct_D3_U16);
	CurveDataD3U16 * pcd = static_cast<CurveDataD3U16 *>(pc);
	
	u16 * pknots = TCurve<u16>::getKnotData(*pcd);
	const u16 * pe = pknots + pcd->nKnots;
	while(pknots != pe)
		sthbrx(pknots++);

	const u32 nctrls = TCurveDc<u16,3>::getControlCount(*pcd);
	u16 * pctrls = TCurve<u16>::getCtrlData(*pcd);
	pe = pctrls + nctrls;
	while(pctrls != pe)
		sthbrx(pctrls++);

	CurveDataD3::fixEndian(pcd);
}
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveTypeD3_U8 -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D3_U8);

u32 CurveTypeD3_U8::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D3_U8);
	return static_cast<const CurveDataD3U8 &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U8::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D3_U8);
	TCurveDc<u8,3>::evaluate(static_cast<const CurveDataD3U8 &>(cd),t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U8::differentiate(const CurveData & cd, float t, u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D3_U8);
	TCurveDc<u8,3>::differentiate(static_cast<const CurveDataD3U8 &>(cd),t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD3_U8::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D3_U8);
	return (flags == 0) ? CurveSamplerTable[cd.nDegree][cst_Point3] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD3_U8::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D3_U8);
	return CurveSupportCache::getSize(cd.nDegree,3);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U8::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D3_U8);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataD3U8 &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U8::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_D3_U8);
	TCurveDc<u8,3>::getKnotCtrlSupport(static_cast<const CurveDataD3U8 &>(cd),t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD3_U8::checkSupport(u32 dim,u32 nknots) const
{
	return (nknots <= 0xFF) && (dim == 3);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD3_U8::getFlags() const
{
	return caf_HiCompression;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD3_U8::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree,
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 3));
	return CurveDataD3U8::create(pknots,pctrls,nknots,degree,ct_D3_U8,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U8::calcError(
						const CurveData & cd, 
						const float *     psamples, 
						u32               nsamples, 
						float             start_time,
						float             time_scale,
						u32               flags,
						CurveError *      pce,
						const float *     ptimes
						) const 
{
	ASSERT(cd.nType == ct_D3_U8);
	TCurveDc<u8,3>::calcError(static_cast<const CurveDataD3U8 &>(cd),psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_U8::fixEndian(CurveData * pcd) const
{
	ASSERT(pcd->nType == ct_D3_U8);
	CurveDataD3::fixEndian(static_cast<CurveDataD3 *>(pcd));
}
///////////////////////////////////////////////////////////////////////////////
//
// ---------------------------- CurveTypeD4n_U16 ------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D4n_U16);

u32 CurveTypeD4n_U16::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	COMPILE_TIME_ASSERT(sizeof(CurveDataD4nU16) <= MAX_CURVE_HEADER_SIZE);
	return static_cast<const CurveDataD4nU16 &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U16::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	TCurveDc<u16,4>::evaluate(static_cast<const CurveDataD4nU16 &>(cd),t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U16::differentiate(const CurveData & cd, float t,u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	TCurveDc<u16,4>::differentiate(static_cast<const CurveDataD4nU16 &>(cd),t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD4n_U16::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	return (flags == cef_Quaternion) ? CurveSamplerTable[cd.nDegree][cst_Quat] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4n_U16::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	return CurveSupportCache::getSize(cd.nDegree,4);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U16::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataD4nU16 &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U16::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	TCurveDc<u16,4>::getKnotCtrlSupport(static_cast<const CurveDataD4nU16 &>(cd),t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD4n_U16::checkSupport(u32 dim,u32 nknots) const
{
	return (nknots <= 0xFFFF) && (dim == 4);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4n_U16::getFlags() const
{
	return caf_NormalizedCtrls;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD4n_U16::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree, 
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 4));
	return CurveDataD4nU16::create(pknots,pctrls,nknots,degree,flags,ct_D4n_U16,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U16::calcError(
							const CurveData & cd, 
							const float *     psamples, 
							u32               nsamples, 
							float             start_time,
							float             time_scale,
							u32               flags,
							CurveError *      pce,
							const float *     ptimes
							) const
{
	ASSERT(cd.nType == ct_D4n_U16);
	TCurveDc<u16,4>::calcError(static_cast<const CurveDataD4nU16 &>(cd),psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U16::fixEndian(CurveData * pc) const
{
	ASSERT(pc->nType == ct_D4n_U16);
	CurveDataD4nU16 * pcd = static_cast<CurveDataD4nU16 *>(pc);

	u16 * pknots = TCurve<u16>::getKnotData(*pcd);
	const u16 * pe = pknots + pcd->nKnots;
	while(pknots != pe)
		sthbrx(pknots++);

	const u32 nctrls = TCurveDc<u16,3>::getControlCount(*pcd);
	u16 * pctrls = TCurve<u16>::getCtrlData(*pcd);
	pe = pctrls + nctrls;
	while(pctrls != pe)
		sthbrx(pctrls++);

	CurveDataD4nU16::fixEndian(pcd);
}
///////////////////////////////////////////////////////////////////////////////
//
// ---------------------------- CurveTypeD4n_U8 -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D4n_U8);

u32 CurveTypeD4n_U8::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	COMPILE_TIME_ASSERT(sizeof(CurveDataD4nU8) <= MAX_CURVE_HEADER_SIZE);
	return static_cast<const CurveDataD4nU8 &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U8::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	TCurveDc<u8,4>::evaluate(static_cast<const CurveDataD4nU8 &>(cd),t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U8::differentiate(const CurveData & cd, float t,u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	TCurveDc<u8,4>::differentiate(static_cast<const CurveDataD4nU8 &>(cd),t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD4n_U8::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	return (flags == cef_Quaternion) ? CurveSamplerTable[cd.nDegree][cst_Quat] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4n_U8::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	return CurveSupportCache::getSize(cd.nDegree,4);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U8::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataD4nU8 &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U8::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	TCurveDc<u8,4>::getKnotCtrlSupport(static_cast<const CurveDataD4nU8 &>(cd),t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD4n_U8::checkSupport(u32 dim, u32 nknots) const
{
	return (nknots <= 0xFF) && (dim == 4);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4n_U8::getFlags() const
{
	return (caf_NormalizedCtrls | caf_HiCompression);
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD4n_U8::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree, 
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 4));
	return CurveDataD4nU8::create(pknots,pctrls,nknots,degree,flags,ct_D4n_U8,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U8::calcError(
						const CurveData & cd, 
						const float *     psamples, 
						u32               nsamples, 
						float             start_time,
						float             time_scale,
						u32               flags,
						CurveError *      pce,
						const float *     ptimes
						) const
{
	ASSERT(cd.nType == ct_D4n_U8);
	TCurveDc<u8,4>::calcError(static_cast<const CurveDataD4nU8 &>(cd),psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4n_U8::fixEndian(CurveData * pcd) const
{
	ASSERT(pcd->nType == ct_D4n_U8);
	CurveDataD4n::fixEndian(static_cast<CurveDataD4n *>(pcd));
}
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveTypeD4_U16 ------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D4_U16);

u32 CurveTypeD4_U16::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4_U16);
	COMPILE_TIME_ASSERT(sizeof(CurveDataD4U16) <= MAX_CURVE_HEADER_SIZE);
	return static_cast<const CurveDataD4U16 &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U16::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4_U16);
	TCurveDc<u16,4>::evaluate(static_cast<const CurveDataD4U16 &>(cd),t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U16::differentiate(const CurveData & cd, float t,u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4_U16);
	TCurveDc<u16,4>::differentiate(static_cast<const CurveDataD4U16 &>(cd),t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD4_U16::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D4_U16);
	return (flags == cef_Quaternion) ? CurveSamplerTable[cd.nDegree][cst_Quat] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4_U16::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4_U16);
	return CurveSupportCache::getSize(cd.nDegree,4);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U16::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D4_U16);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataD4U16 &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U16::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_D4_U16);
	TCurveDc<u16,4>::getKnotCtrlSupport(static_cast<const CurveDataD4U16 &>(cd),t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD4_U16::checkSupport(u32 dim,u32 nknots) const
{
	return (nknots <= 0xFFFF) && (dim == 4);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4_U16::getFlags() const
{
	return 0;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD4_U16::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree, 
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 4));
	return CurveDataD4U16::create(pknots,pctrls,nknots,degree,ct_D4_U16,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U16::calcError(
						const CurveData & cd, 
						const float *     psamples, 
						u32               nsamples, 
						float             start_time, 
						float             time_scale,
						u32               flags,
						CurveError *      pce,
						const float *     ptimes
						) const 
{
	ASSERT(cd.nType == ct_D4_U16);
	TCurveDc<u16,4>::calcError(static_cast<const CurveDataD4U16 &>(cd),psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U16::fixEndian(CurveData * pc) const
{
	ASSERT(pc->nType == ct_D4_U16);
	CurveDataD4U16 * pcd = static_cast<CurveDataD4U16 *>(pc);
	
	u16 * pknots = TCurve<u16>::getKnotData(*pcd);
	const u16 * pe = pknots + pcd->nKnots;
	while(pknots != pe)
		sthbrx(pknots++);

	const u32 nctrls = TCurveDc<u16,4>::getControlCount(*pcd);
	u16 * pctrls = TCurve<u16>::getCtrlData(*pcd);
	pe = pctrls + nctrls;
	while(pctrls != pe)
		sthbrx(pctrls++);

	CurveDataD4::fixEndian(pcd);
}
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveTypeD4_U8 -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D4_U8);

u32 CurveTypeD4_U8::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4_U8);
	return static_cast<const CurveDataD4U8 &>(cd).getSize();
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U8::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4_U8);
	TCurveDc<u8,4>::evaluate(static_cast<const CurveDataD4U8 &>(cd),t,flags,pres);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U8::differentiate(const CurveData & cd, float t, u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4_U8);
	TCurveDc<u8,4>::differentiate(static_cast<const CurveDataD4U8 &>(cd),t,k,pres);
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD4_U8::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D4_U8);
	return (flags == cef_Quaternion) ? CurveSamplerTable[cd.nDegree][cst_Quat] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4_U8::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4_U8);
	return CurveSupportCache::getSize(cd.nDegree,4);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U8::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D4_U8);
	CurveSupportCache::init(cd.nDegree,icurve,emit_offset,static_cast<const CurveDataD4U8 &>(cd).getKnotScale(),pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U8::getKnotCtrlSupport(
						const CurveData &  cd, 
						float              t, 
						u32                flags,
						float RESTRICT_PTR pknots, 
						float RESTRICT_PTR pctrls
						) const
{
	ASSERT(cd.nType == ct_D4_U8);
	TCurveDc<u8,4>::getKnotCtrlSupport(static_cast<const CurveDataD4U8 &>(cd),t,flags,pknots,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD4_U8::checkSupport(u32 dim,u32 nknots) const
{
	return (nknots <= 0xFF) && (dim == 4);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4_U8::getFlags() const
{
	return caf_HiCompression;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD4_U8::createCurve(                                                                                      
								const float * pknots,                                        
								const float * pctrls,                                        
								u32           nknots,                                    
								u32           dimensions,                                    
								u32           degree,
								u32           flags,
								u32 *         psize,                                         
								byte *        pmem_loc,
								u32           mem_loc_size
								) const 
{
	ASSERT(pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 4));
	return CurveDataD4U8::create(pknots,pctrls,nknots,degree,ct_D4_U8,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U8::calcError(
						const CurveData & cd, 
						const float *     psamples, 
						u32               nsamples, 
						float             start_time,
						float             time_scale,
						u32               flags,
						CurveError *      pce,
						const float *     ptimes
						) const 
{
	ASSERT(cd.nType == ct_D4_U8);
	TCurveDc<u8,4>::calcError(static_cast<const CurveDataD4U8 &>(cd),psamples,nsamples,start_time,time_scale,flags,pce,ptimes);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_U8::fixEndian(CurveData * pcd) const
{
	ASSERT(pcd->nType == ct_D4_U8);
	CurveDataD4::fixEndian(static_cast<CurveDataD4 *>(pcd));
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------- CurveTypeD3_ConstF32 -----------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D3_ConstF32);

u32 CurveTypeD3_ConstF32::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	COMPILE_TIME_ASSERT(sizeof(const CurveDataC3F32) <= MAX_CURVE_HEADER_SIZE);
	return sizeof(static_cast<const CurveDataC3F32 &>(cd));
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_ConstF32::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	const CurveDataC3F32 & tcd = static_cast<const CurveDataC3F32 &>(cd);
	pres->aCurve[0] = tcd.aConstVals[0];
	pres->aCurve[1] = tcd.aConstVals[1];
	pres->aCurve[2] = tcd.aConstVals[2];
	ASSERT((flags & cef_Normalize) == 0 || IsEqual(Len2VectorN(pres->aCurve,3),1.f));
	ASSERT((flags & cef_Quaternion) == (flags & cef_Normalize));
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_ConstF32::differentiate(const CurveData & cd,float t,u32 k, CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	pres->aCurve[0] = 0.f;
	pres->aCurve[1] = 0.f;
	pres->aCurve[2] = 0.f;
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD3_ConstF32::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	return (flags == 0) ? CurveSamplerTable[0][cst_Point3] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD3_ConstF32::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	return sizeof(CurveCache);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_ConstF32::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	CurveCache::init(icurve,emit_offset,pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_ConstF32::getKnotCtrlSupport(
							const CurveData &  cd, 
							float              t, 
							u32                flags,
							float RESTRICT_PTR pknots, 
							float RESTRICT_PTR pctrls
							) const
{
	ASSERT_FAILED_MSG("CurveTypeD3_ConstF32::getKnotCtrlSupport should never be called.\n");
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD3_ConstF32::checkSupport(u32 dim,u32 nknots) const
{
	return (dim == 3);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD3_ConstF32::getFlags() const
{
	return caf_Constant;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD3_ConstF32::createCurve(                                                                                      
									const float * pknots,                                        
									const float * pctrls,                                        
									u32           knot_count,                                    
									u32           dimensions,                                    
									u32           degree, 
									u32           flags,
									u32 *         psize,                                         
									byte *        pmem_loc,
									u32           mem_loc_size
									) const 
{
	ASSERT(pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 3));
	return CurveDataC3F32::create(pctrls,degree,ct_D3_ConstF32,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_ConstF32::calcError(
								const CurveData & cd,
								const float *     psamples,
								u32               nsamples, 
								float             start_time,
								float             time_scale,
								u32               flags,
								CurveError *      pce,
								const float *     ptimes
								) const 
{
	ASSERT(cd.nType == ct_D3_ConstF32);
	static_cast<const CurveDataC3F32 &>(cd).calcError(psamples,nsamples,pce);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD3_ConstF32::fixEndian(CurveData * pcd) const
{
	ASSERT(pcd->nType == ct_D3_ConstF32);
	CurveDataC3F32::fixEndian(static_cast<CurveDataC3F32 *>(pcd));
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------- CurveTypeD4_ConstF32 -----------------------------
//
///////////////////////////////////////////////////////////////////////////////
DECLARE_BUILTIN_CURVE_TYPE(D4_ConstF32);

u32 CurveTypeD4_ConstF32::getSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	COMPILE_TIME_ASSERT(sizeof(const CurveDataC4F32) <= MAX_CURVE_HEADER_SIZE);
	return sizeof(static_cast<const CurveDataC4F32 &>(cd));
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_ConstF32::evaluate(const CurveData & cd, float t,u32 flags,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	const CurveDataC4F32 & tcd = static_cast<const CurveDataC4F32 &>(cd);
	pres->aCurve[0] = tcd.aConstVals[0];
	pres->aCurve[1] = tcd.aConstVals[1];
	pres->aCurve[2] = tcd.aConstVals[2];
	pres->aCurve[3] = tcd.aConstVals[3];
	ASSERT((flags & cef_Normalize) == 0 || IsEqual(Len2VectorN(pres->aCurve,4),1.f));
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_ConstF32::differentiate(const CurveData & cd, float t,u32 k,CurveEvalBuf * pres) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	pres->aCurve[0] = 0.f;
	pres->aCurve[1] = 0.f;
	pres->aCurve[2] = 0.f;
	pres->aCurve[3] = 0.f;
}
///////////////////////////////////////////////////////////////////////////////
CurveSampler CurveTypeD4_ConstF32::getSampler(const CurveData & cd,u32 flags) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	return (flags == cef_Quaternion) ? CurveSamplerTable[0][cst_Quat] : NULL;
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4_ConstF32::getCacheSize(const CurveData & cd) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	return sizeof(CurveCache);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_ConstF32::initCache(const CurveData & cd,s32 icurve,u32 emit_offset,byte * pb) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	CurveCache::init(icurve,emit_offset,pb);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_ConstF32::getKnotCtrlSupport(
							const CurveData &  cd, 
							float              t, 
							u32                flags,
							float RESTRICT_PTR pknots, 
							float RESTRICT_PTR pctrls
							) const
{
	ASSERT_FAILED_MSG("CurveTypeD4_ConstF32::getKnotCtrlSupport should never be called.\n");
}
///////////////////////////////////////////////////////////////////////////////
bool CurveTypeD4_ConstF32::checkSupport(u32 dim,u32 nknots) const
{
	return (dim == 4);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveTypeD4_ConstF32::getFlags() const
{
	return caf_Constant;
}
///////////////////////////////////////////////////////////////////////////////
CurveData * CurveTypeD4_ConstF32::createCurve(                                                                                      
									const float * pknots,                                        
									const float * pctrls,                                        
									u32           knot_count,                                    
									u32           dimensions,                                    
									u32           degree,
									u32           flags,
									u32 *         psize,                                         
									byte *        pmem_loc,
									u32           mem_loc_size
									) const 
{
	ASSERT(!pknots && pctrls);
	ASSERT((degree <= MAX_CURVE_DEGREE) && (dimensions == 4));
	return CurveDataC4F32::create(pctrls,degree,ct_D4_ConstF32,psize,pmem_loc,mem_loc_size);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_ConstF32::calcError(
								const CurveData & cd, 
								const float *     psamples, 
								u32               nsamples, 
								float             start_time,
								float             time_scale,
								u32               flags,
								CurveError *      pce,
								const float *     ptimes
								) const
{
	ASSERT(cd.nType == ct_D4_ConstF32);
	static_cast<const CurveDataC4F32 &>(cd).calcError(psamples,nsamples,pce);
}
///////////////////////////////////////////////////////////////////////////////
void CurveTypeD4_ConstF32::fixEndian(CurveData * pcd) const
{
	ASSERT(pcd->nType == ct_D4_ConstF32);
	CurveDataC4F32::fixEndian(static_cast<CurveDataC4F32 *>(pcd));
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveDataDv --------------------------------
//
///////////////////////////////////////////////////////////////////////////////
void CurveDataDv::getKnotCtrlSupport(float & t,u32 flags,float RESTRICT_PTR pknots,float RESTRICT_PTR pctrls) const
{
	ASSERT(t >= 0.f && t <= 1.f);
	const s32 ki = getKnotIndex(t);
	getKnotSupport(ki,pknots);
	getCtrlSupport(ki,flags,pctrls);
}
///////////////////////////////////////////////////////////////////////////////
void CurveDataDv::evaluate(float t,u32 flags,CurveEvalBuf * pres) const
{
	getKnotCtrlSupport(t,flags,pres->aKnots,pres->aCtrls);
	CurveBasis(t,nDegree,pres->aKnots,pres->aBasis);
	
	const u32 dim = (nDimAndCtrl >> nDimShift);
	CurveApplyBasis(pres->aCtrls,pres->aBasis,dim,nDegree,pres->aCurve);

	if(flags & cef_Normalize)
		NormalizeVectorSafeN(pres->aCurve,dim);
}
///////////////////////////////////////////////////////////////////////////////
void CurveDataDv::differentiate(float t, u32 k, CurveEvalBuf * pres) const
{
	getKnotCtrlSupport(t,0,pres->aKnots,pres->aCtrls);
	CurveDeriv(t,nDegree,k,pres->aKnots,pres->aBasis);
	TCurve<float>::scaleBasisDerivatives(*this,k,pres->aBasis);

	const u32 dim = (nDimAndCtrl >> nDimShift);
	CurveApplyBasis(pres->aCtrls,pres->aBasis,dim,nDegree,pres->aCurve);
}
///////////////////////////////////////////////////////////////////////////////
void CurveDataDv::calcError(
					const float * psamples,
					u32           nsamples,
					float         start_time,
					float         time_scale,
					u32           flags,
					CurveError *  pce,
					const float * ptimes
					) const
{
	const float oo_time_scale = (time_scale > 0.f) ? 1.f / time_scale : 0.f;
	const u32  dim = nDimAndCtrl >> nDimShift;

	CurveEvalBuf buf;
	pce->fTotErrSq = 0.f;
	pce->fMaxErrSq = 0.f;
	pce->iMaxErrSample = -1;
	
	for(u32 i = 0; i < nsamples; i++, psamples += dim)
	{
		float t = ptimes ? ptimes[i] : float(i);
		t = (t + start_time) * oo_time_scale;
		evaluate(t,flags,&buf);

		float e2 = CurveCalcErrorSq(buf.aCurve,psamples,dim,flags);
		pce->fTotErrSq += e2;
		if(e2 > pce->fMaxErrSq)
		{
			pce->fMaxErrSq = e2;
			pce->iMaxErrSample = (s32)i;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
CurveDataDv * CurveDataDv::create(	
							const float * pknots,
							const float * pctrls,
							u32           nknots,
							u32           dimensions,
							u32           degree,
							u32 *         psize,
							byte *        pmem_loc,
							u32           mem_loc_size
							)
{
	const u32 nctrls = (nknots + degree - 1) * dimensions;
	const u32 nbsize = TCurve<float>::getSize<CurveDataDv>(nknots,nctrls);

	if(!pmem_loc)
	{
		pmem_loc = (byte *)CurveType::alloc(nbsize);
		mem_loc_size = nbsize;
	}
	
	if(nbsize > mem_loc_size)
		return NULL;
	
	CurveDataDv * pcd = reinterpret_cast<CurveDataDv *>(pmem_loc);
	pcd->nType       = ct_Dv_F32;
	pcd->nDegree     = (u8)degree;
	pcd->nKnots      = (u16)nknots;
	pcd->nDimAndCtrl = (dimensions << CurveDataDv::nDimShift) | nctrls;

	float RESTRICT_PTR pk = TCurve<float>::getKnotData(*pcd);
	const float * pe = pk + nknots;
	while(pk != pe)
		*pk++ = *pknots++;
	
	float RESTRICT_PTR pc = TCurve<float>::getCtrlData(*pcd);
	pe = pc + nctrls;
	while(pc != pe)
		*pc++ = *pctrls++;

	ASSERT((u32(pc) - u32(pcd)) == nbsize);
	if(psize)
		*psize = nbsize;

	return pcd;
}
///////////////////////////////////////////////////////////////////////////////
void CurveDataDv::fixEndian(CurveDataDv * pcd)
{
	float * pk = TCurve<float>::getKnotData(*pcd);
	const float * pe = pk + pcd->nKnots;
	while(pk != pe)
		stwbrx((u32 *)pk++);

	const u32 nctrls = pcd->nDimAndCtrl & nCtrlMask;
	float * pc = TCurve<float>::getCtrlData(*pcd);
	pe = pc + nctrls;
	while(pc != pe)
		stwbrx((u32 *)pc++);

	sthbrx(&pcd->nKnots);
	stwbrx(&pcd->nDimAndCtrl);
}
///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveDataD4n -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
const float CurveDataD4n::U4ToF32[16] = 
{
	0.f,1.f,2.f,3.f,4.f,5.f,6.f,7.f,8.f,9.f,10.f,11.f,12.f,13.f,14.f,15.f
};
///////////////////////////////////////////////////////////////////////////////




