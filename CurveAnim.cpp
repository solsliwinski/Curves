#include "Support/SupportHeader.h"
#include "Shared/BitUtils.h"
#include "Curve.h"
#include "CurveAnim.h"

#if defined(PWTOOL)
#include "Shared/QuickSort.h"
#include "Shared/BufferTokenizer.h"
#endif

IMPLEMENT_LOG_CATEGORY("Anim", lg_Anim, true);

COMPILE_TIME_ASSERT(nu_HasRotEntry == 1);		
COMPILE_TIME_ASSERT(nu_HasTransEntry == 2);
COMPILE_TIME_ASSERT(nu_RefRotEntry == 4);
COMPILE_TIME_ASSERT(nu_RefTransEntry == 8);

///////////////////////////////////////////////////////////////////////////////
//
// ------------------------------- CurveAnim ----------------------------------
//
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::calcNTP(s32 inode, QuatXForm * pxf, float fpos) const
{
	ASSERT((getNodeUsage(inode) & nu_HasRotAndTrans) == nu_HasRotAndTrans);
	calcNTPQuat(inode,&pxf->quatPart, fpos);
	calcNTPTrans(inode, &pxf->transPart, fpos);
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::calcNTPQuat(s32 inode, Quat * pquat, float fpos) const
{
	const u32 nu = getNodeUsage(inode);
	ASSERT(nu & nu_HasRotEntry);

	if(nu & nu_RefRotEntry)
	{
		*pquat = pNodeTree->getRefNTP(inode).quatPart;
		return;
	}

	s32 icurve = getCurveIndex(inode | cac_QuatMask);
	ASSERT(icurve != -1);

	const CurveData & cd = *ppCurves[icurve];
	DEBUG_ASSERT(CurveType::get(cd)->checkSupport(4,cd.nKnots));

	CurveEvalBuf buf;
	CurveType::get(cd)->evaluate(cd,fpos,cef_Quaternion,&buf);
	pquat->pt.set(buf.aCurve[0],buf.aCurve[1],buf.aCurve[2],buf.aCurve[3]);
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::calcNTPTrans(s32 inode, Point3 * ppt, float fpos) const
{
	const u32 nu = getNodeUsage(inode);
	ASSERT(nu & nu_HasTransEntry);

	if(nu & nu_RefTransEntry)
	{
		*ppt = pNodeTree->getRefNTP(inode).transPart;
		return;
	}

	s32 icurve = getCurveIndex(inode);
	ASSERT(icurve != -1);

	const CurveData & cd = *ppCurves[icurve];
	DEBUG_ASSERT(CurveType::get(cd)->checkSupport(3,cd.nKnots));

	CurveEvalBuf buf;
	CurveType::get(cd)->evaluate(cd,fpos,0,&buf);
	ppt->set(buf.aCurve[0],buf.aCurve[1],buf.aCurve[2]);
}
///////////////////////////////////////////////////////////////////////////////
u32 CurveAnim::nodeUsage(int inode) const
{
	return getNodeUsage(inode);
}
///////////////////////////////////////////////////////////////////////////////
s32 CurveAnim::getCurveIndex(u32 inode) const
{
	for(u32 i = 0; i < nCurves; i++)
	{
		if(pNodeIds[i] == inode)
			return i;
	}

	return -1;
}
///////////////////////////////////////////////////////////////////////////////
int CurveAnim::getMaxInstSize() const
{
	return nMaxInstSize;
}
///////////////////////////////////////////////////////////////////////////////
int CurveAnim::getInstSize(int nnodes, const u32 * pnode_usage) const
{
	u32 nb_inst = sizeof(InstHeader);

	int i, n;
	for(i = 0; i < nnodes; i++)
	{
		if(pnode_usage && !TestBitArrayBit(pnode_usage,i))
			continue;

		u32 nu = getNodeUsage(i);
		if(nu & nu_RefTransEntry)
			nb_inst += sizeof(s16);

		if(nu & nu_RefRotEntry)
			nb_inst += sizeof(s16);
	}

	nb_inst = pw_align<4>(nb_inst);

	CurveSampler pls = NULL;
	u32 nb_cache = 0;
	for(i = 0, n = nCurves; i < n; i++)
	{
		u32 inode = pNodeIds[i] & cac_NodeMask;
		if(pnode_usage && !TestBitArrayBit(pnode_usage,inode))
			continue;

		const CurveData & curve = *ppCurves[i];
		const CurveType * ptype = CurveType::get(curve);

		CurveSampler ps = ptype->getSampler(curve,(pNodeIds[i] & cac_QuatMask) ? cef_Quaternion : 0);
		ASSERT(ps);

		if(ps != pls)
		{
			pls = ps;
			nb_cache = ptype->getCacheSize(curve);
			nb_inst += sizeof(CurveSamplerCache);
		}

		nb_inst += nb_cache;
	}

	ASSERT(nb_inst == pw_align<4>(nb_inst));
	ASSERT(nb_inst <= (u32)getMaxInstSize());
	return nb_inst;
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::initInst(byte * pb, int nnodes, const u32 * pnode_usage) const
{
	InstHeader * phdr = (InstHeader *)pb;
	phdr->nPosRef   = 0;
	phdr->nRotRef   = 0;
	phdr->nSamplers = 0;
	phdr->nTotSize  = sizeof(InstHeader);

	int i,n;
	for(i = 0; i < nnodes; i++)
	{
		if(pnode_usage && !TestBitArrayBit(pnode_usage,i))
			continue;

		u32 ub = getNodeUsage(i);
		if(ub & nu_RefTransEntry)
		{
			ASSERT(ub & nu_HasTransEntry);
			ASSERT(phdr->nPosRef < 255);
			phdr->nPosRef++;
		}

		if(ub & nu_RefRotEntry)
		{
			ASSERT(ub & nu_HasRotEntry);
			ASSERT(phdr->nRotRef < 255);
			phdr->nRotRef++;
		}
	}

	pb += sizeof(InstHeader);
	s16 * ppos_ref_ids = (s16 *)pb;
	s16 * prot_ref_ids = ppos_ref_ids + phdr->nPosRef;

	for(i = 0; i < nnodes; i++)
	{
		if(pnode_usage && !TestBitArrayBit(pnode_usage,i))
			continue;

		u32 ub = getNodeUsage(i);
		if(ub & nu_RefTransEntry)
			*ppos_ref_ids++ = (s16)i;

		if(ub & nu_RefRotEntry)
			*prot_ref_ids++ = (s16)i;
	}

	pb += pw_align<4>(sizeof(s16) * (phdr->nRotRef + phdr->nPosRef));
	phdr->nTotSize += pw_align<4>(sizeof(s16) * (phdr->nRotRef + phdr->nPosRef));

	CurveSamplerCache * pcache = (CurveSamplerCache *)pb;
	pcache->pSampler = NULL;
	u32 cache_size = 0;

	for(i = 0, n = nCurves; i < n; i++)
	{
		s32 inode = pNodeIds[i] & cac_NodeMask;
		if(pnode_usage && !TestBitArrayBit(pnode_usage,inode))
			continue;

		const CurveData & curve = *ppCurves[i];
		const CurveType * ptype = CurveType::get(curve);

		CurveSampler pcs = ptype->getSampler(curve,(pNodeIds[i] & cac_QuatMask) ? cef_Quaternion : 0);
		ASSERT(pcs);

		if(pcs != pcache->pSampler)
		{
			pcache = (CurveSamplerCache *)pb;
			pcache->pSampler = pcs;
			pcache->nEntries = 0;
			pcache->nTotSize = sizeof(CurveSamplerCache);
			
			phdr->nTotSize += sizeof(CurveSamplerCache);
			phdr->nSamplers++;

			cache_size = ptype->getCacheSize(curve);
			pb += sizeof(CurveSamplerCache);
		}

		ptype->initCache(curve,i,inode * sizeof(QuatXForm),pb);
		pcache->nEntries++;
		pcache->nTotSize += (u16)cache_size;
		phdr->nTotSize += cache_size;
		pb += cache_size;
	}

	DEBUG_ASSERT(phdr->nTotSize <= (u32)getMaxInstSize());
	DEBUG_ASSERT(phdr->nTotSize == (u32)getInstSize(nnodes,pnode_usage));
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::shutdownInst(byte * pinst_data) const
{
#if defined(_DEBUG)
	InstHeader * phdr = (InstHeader *)pinst_data;
	memset(pinst_data,0xbe,phdr->nTotSize);
#endif
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::calcNTPInst(QuatXForm * pxfs, float fpos, byte * pinst_data) const
{
	const InstHeader * phdr = (InstHeader *)pinst_data;
	s16 RESTRICT_PTR pref_ids = (s16 *)(pinst_data + sizeof(InstHeader));
	
	u32 i,n;
	for(i = 0, n = phdr->nPosRef; i < n; i++)
	{
		s16 id = *pref_ids++;
		pxfs[id].transPart = pNodeTree->getRefNTP(id).transPart;
	}

	for(i = 0, n = phdr->nRotRef; i < n; i++)
	{
		s16 id = *pref_ids++;
		pxfs[id].quatPart = pNodeTree->getRefNTP(id).quatPart;
	}

	pinst_data += sizeof(InstHeader) + pw_align<4>(sizeof(s16) * (phdr->nPosRef + phdr->nRotRef));

	for(i = 0, n = phdr->nSamplers; i < n; i++)
	{
		CurveSamplerCache * pcache = (CurveSamplerCache *)pinst_data;
		pcache->pSampler(fpos,ppCurves,pcache->nEntries,(byte *)pxfs,pinst_data + sizeof(CurveSamplerCache));
		pinst_data += pcache->nTotSize;
	}
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::fixup(NodeTree * pnt)
{
	NodeTreeAnim::fixup(pnt);

	byte * pthis = (byte *)this;

	ppCurves   = (const CurveData **)(pthis + nCurvesOffset);
	pNodeUsage = (const u32 *)(pthis + nNodeUsageOffset);
	pNodeIds   = (const u16 *)(pthis + nNodeIdsOffset);

	for(u32 i = 0, n = nCurves; i < n; i++)
	{
		ppCurves[i] = (const CurveData *)(pthis + ((u32)ppCurves[i]));
		ASSERT(ppCurves[i]->nType <= ct_LastType);
	}
}
///////////////////////////////////////////////////////////////////////////////
CurveAnim * CurveAnim::initCurveAnimFromImage(byte * pdata, NodeTree * pnt)
{
	CurveAnim * pca = new(pdata) CurveAnim;
	pca->fixup(pnt);
	return pca;
}
///////////////////////////////////////////////////////////////////////////////
#if defined(PWTOOL)
CurveAnim * CurveAnim::makeCurveAnimFromAssetFile(AssetFileResourcePtr passet,const CurveAnimMaker & cam, NodeTree * pnt)
{
	ASSERT(passet->getNodeTree());

	u32 nsize = sizeof(CurveAnim) + 
		cam.getCurveCount() * sizeof(CurveData *) + 
		pw_align<4>(cam.getCurveCount() * sizeof(s16)) + 
		cam.getNodeUsageSize() + 
		cam.getCurveSize();

	byte * pdata = (byte *)pw_aligned_malloc(nsize,16);
	memset(pdata,0,nsize);

	CurveAnim * panim = new (pdata) CurveAnim;
	panim->loadAssetFile(passet,pnt,sizeof(CurveAnim));

	panim->nCurves      = (u16)cam.getCurveCount();
	panim->nNodes       = (u16)pnt->noNodes();
	panim->nMaxInstSize = sizeof(CurveAnim::InstHeader) + pw_align<4>(sizeof(s16) * cam.getRefCount());
	panim->nTrueSize    = nsize;
	panim->animType     = na_Curve;
	panim->ppCurves     = (const CurveData **)(pdata + sizeof(CurveAnim));

	const CurveData * const * ppcurves = cam.getCurves();
	byte * pb = pdata + sizeof(CurveAnim) + sizeof(CurveData *) * cam.getCurveCount();

	panim->pNodeIds = (u16 *)pb;

	int i;
	for(i = 0; i < (int)cam.getCurveCount(); i++)
	{
		(u16 &)(panim->pNodeIds[i]) = (u16)cam.getNodeId(i);
	}

	pb += pw_align<4>(sizeof(u16) * cam.getCurveCount());

	panim->pNodeUsage = (u32 *)pb;
	memcpy(pb,cam.getNodeUsage(),cam.getNodeUsageSize());

	pb += cam.getNodeUsageSize();

	CurveSampler pls = NULL;
	for(i = 0; i < panim->nCurves; i++)
	{
		panim->ppCurves[i] = (CurveData *)pb;

		const CurveType * ptype = CurveType::get(*ppcurves[i]);
		ASSERT(ptype);

		u32 curve_size = ptype->getSize(*ppcurves[i]);
		memcpy(pb,ppcurves[i],curve_size);
		pb += curve_size;

		u32 flags = (cam.getNodeId(i) & cac_QuatMask) ? cef_Quaternion : 0;
		CurveSampler ps = ptype->getSampler(*ppcurves[i],flags);
		if(ps != pls)
		{
			panim->nMaxInstSize += sizeof(CurveSamplerCache);
			pls = ps;
		}

		ASSERT(pls);
		panim->nMaxInstSize += ptype->getCacheSize(*ppcurves[i]);
	}

	ASSERT((u32)pb - (u32)pdata == panim->nTrueSize);
	return panim;
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::fixEndian(byte * pdata)
{
	CurveAnim * pca = (CurveAnim *)pdata;
	u32 * pcurve_offsets = (u32 *)(pdata + pca->nCurvesOffset);

	u32 i,n;
	for(i = 0, n = pca->nCurves; i < n; i++)
	{
		CurveData & cd = *(CurveData *)(pdata + pcurve_offsets[i]);
		const CurveType * ptype = CurveType::get(cd);
		ASSERT(ptype);
		ptype->fixEndian(&cd);
		stwbrx(pcurve_offsets + i);
	}

	stwbrx(&pca->nCurvesOffset);

	// pNodeIds is curve-to-node index translation table, so iterate over number of curves.
	s16 * pnode_ids = (s16 *)(pdata + pca->nNodeIdsOffset);
	for(i = 0; i < (u32)pca->nCurves; i++)
	{
		sthbrx(pnode_ids + i);
	}

	u32 * pnode_usage = (u32 *)(pdata + pca->nNodeUsageOffset);
	u32 ndw = (pca->nNodes * 4 + 31) >> 5;
	for(i = 0; i < ndw; i++)
	{
		stwbrx(pnode_usage + i);
	}

	stwbrx(&pca->nNodeIdsOffset);
	stwbrx(&pca->nNodeUsageOffset);

	sthbrx(&pca->nCurves);
	sthbrx(&pca->nNodes);
	stwbrx(&pca->nMaxInstSize);
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnim::finalWritePrep(u32 id_node_tree)
{
	NodeTreeAnim::finalWritePrep(id_node_tree);

	byte * pthis = (byte *)this;

	for(u32 i = 0, n = nCurves; i < n; i++)
	{
		u32 offset = u32(ppCurves[i]) - u32(pthis);
		ppCurves[i] = (CurveData *)offset;
	}

	nCurvesOffset = u32(ppCurves) - u32(pthis);
	nNodeIdsOffset = u32(pNodeIds) - u32(pthis);
	nNodeUsageOffset = u32(pNodeUsage) - u32(pthis);
}
///////////////////////////////////////////////////////////////////////////////
//
// ----------------------------- CurveAnimMaker -------------------------------
//
///////////////////////////////////////////////////////////////////////////////
#define CURVE_ERR_TOL_TOKEN "err_tol"
#define CURVE_CUT_TOL_TOKEN "cut_tol"
#define CURVE_DEGREE_TOKEN  "degree"
#define CURVE_FLAGS_TOKEN   "flags"
#define CURVE_TYPES_TOKEN   "types"
#define CURVE_POS_SUFFIX    "_pos"
#define CURVE_ROT_SUFFIX    "_rot"
#define CURVE_REC_SUFFIX    "_rec"
#define CURVE_DEFAULT_TOKEN "default"
#define ADD_NAMED_ENUM(enum_val) { #enum_val, enum_val }
///////////////////////////////////////////////////////////////////////////////
struct CurveSortRec
{
	CurveData * pCurve;
	u32         nKeyVal;

	void init(CurveData * pcurve,u32 inode)
	{
		const u32 dim = (inode & cac_QuatMask) ? 4 : 3;
		nKeyVal = (dim << 24) | (pcurve->nDegree << 16) | inode;
		pCurve  = pcurve;
		ASSERT((inode & cac_QuatMask) == (getNodeId() & cac_QuatMask));
	}

	u32 getNodeId() const
	{
		return nKeyVal & 0xFFFF;
	}
};
///////////////////////////////////////////////////////////////////////////////
class CurveSortRecComparator : public Comparator<const CurveSortRec,const CurveSortRec>
{
public:
	
	static int compare(const CurveSortRec & a, const CurveSortRec & b)
	{
#if defined(_DEBUG)
		u32 n0 = a.getNodeId();
		u32 n1 = b.getNodeId();
		u32 f0 = (n0 & cac_QuatMask) ? cef_Quaternion : 0;
		u32 f1 = (n1 & cac_QuatMask) ? cef_Quaternion : 0;
		CurveSampler s0 = CurveType::get(*a.pCurve)->getSampler(*a.pCurve,f0);
		CurveSampler s1 = CurveType::get(*b.pCurve)->getSampler(*b.pCurve,f1);
		ASSERT(s0 && s1);
		ASSERT(s0 == s1 || n0 != n1 || f0 != f1);
#endif
		if(a.nKeyVal < b.nKeyVal) 
		{
			return -1;
		}
		else if(a.nKeyVal == b.nKeyVal) 
		{
			return 0;
		}
		else 
		{
			return 1;
		}
	}
};
///////////////////////////////////////////////////////////////////////////////
void CurveAnimMaker::setNodeTypes(NodeProps * pprops,int * ptypes,int ntypes,bool brec,const NodeData * pnd)
{
	ASSERT(pnd);
	s32 inode = pnd->idxNode;
	ASSERT(inode >= 0 && inode < nNodes);
	pprops[inode].aCurveTypes.resize(ntypes);

	for(int i = 0; i < ntypes; i++)
	{
		pprops[inode].aCurveTypes[i] = CurveType::get(ptypes[i]);
		ASSERT(pprops[inode].aCurveTypes[i]);
	}

	if(brec)
	{
		for(int i = 0; i < pnd->nChildren; i++)
			setNodeTypes(pprops,ptypes,ntypes,brec,pnd->pChildren[i]);
	}
}
///////////////////////////////////////////////////////////////////////////////
bool CurveAnimMaker::loadSettingsFile(const MakeAnimSettings & mas)
{
	struct named_enum
	{
		const char * pName;
		u32          nEnum;
	};

	nNodes = mas.pNodeTree->noNodes();
	if(nNodes == 0)
	{
		char buf[256];
		mas.pNodeTree->getName(buf,sizeof(buf));
		LOG(lg_Error,"CurveAnimMaker::loadSettingsFile failed because \"%s\" has no nodes.\n",buf);
		return false;
	}

	FileLock fl;
	if (!gFileManager.lockFile(mas.szSettingsName, &fl))
	{
		LOG(lg_Error,"CurveAnimMaker::loadSettingsFile could not lock file \"%s\"\n", mas.szSettingsName);
		return false;
	}

	aPosProps.setData(new NodeProps[nNodes]);
	aRotProps.setData(new NodeProps[nNodes]);
	
	BufferTokenizer bt((const char *) fl.pData, fl.nLen);
	for(;;)
	{
		TokenRec tr_op, tr_name, tr_val;
		if(!bt.getCToken(&tr_op)) 
		{
			// All done.
			break;
		}

		char buf[256];
		if(!bt.getStringToken(&tr_name))
		{
			tr_op.translateStringToken(buf,sizeof(buf));
			LOG(lg_Error,"CurveAnimMaker::loadSettingsFile failed due to malformed token \"%s\".\n",buf);
			return false;
		}

		tr_name.translateStringToken(buf, sizeof(buf));

		bool bdef = stricmp(buf,CURVE_DEFAULT_TOKEN) == 0;
		const NodeData * pnd = bdef ? mas.pNodeTree->getNodeOrders()[0] : mas.pNodeTree->findNode(buf);
		s32 inode = pnd ? pnd->idxNode : -1;

		if( check_prefix_i_n(tr_op.pToken,CURVE_ERR_TOL_TOKEN,tr_op.nToken) || 
			check_prefix_i_n(tr_op.pToken,CURVE_CUT_TOL_TOKEN,tr_op.nToken) )
		{
			tr_op.translateStringToken(buf,sizeof(buf));
			if(!bt.getFloatToken(&tr_val))
			{
				LOG(lg_Error,"CurveAnimMaker::loadSettingsFile failed due to malformed token \"%s\".\n",buf);
				return false;
			}

			float val = (float)matof(tr_val.pToken, tr_val.nToken);
			if(val < 0)
			{
				LOG(lg_Error,"CurveAnimMaker::loadSettingsFile - Invalid value %f. (must be >= 0).\n", val);
				return false;
			}

			if(!pnd) 
				continue;

			bool bpos = mstristr(buf,CURVE_POS_SUFFIX) != NULL;
			bool bcut = check_prefix_i_n(tr_op.pToken,CURVE_CUT_TOL_TOKEN,tr_op.nToken);

			if(!bpos)
			{
				// Error tolerance was specified in radians, but we need to convert to quat error metric.
				val = CurveRadToQuatErrorMetric(Clamp(val,0.f,Pi));
			}
		
			if(mstristr(buf,CURVE_REC_SUFFIX) || bdef)
			{
				// Default is always considered recursive.
				if(bcut)
					setNodeProps<float,&NodeProps::fCutTol>(bpos ? &aPosProps[0] : &aRotProps[0],val,pnd);
				else
					setNodeProps<float,&NodeProps::fErrTol>(bpos ? &aPosProps[0] : &aRotProps[0],val,pnd);
			}
			else
			{
				if(bcut)
					bpos ? aPosProps[inode].fCutTol = val : aRotProps[inode].fCutTol = val;
				else
					bpos ? aPosProps[inode].fErrTol = val : aRotProps[inode].fErrTol = val;
			}
		}
		else if(check_prefix_i_n(tr_op.pToken,CURVE_DEGREE_TOKEN,tr_op.nToken))
		{
			tr_op.translateStringToken(buf,sizeof(buf));

			if(!bt.getIntToken(&tr_val))
			{
				LOG(lg_Error,"CurveAnimMaker::loadSettingsFile failed due to malformed token \"%s\".\n",buf);
				return false;
			}

			int val = matoi(tr_val.pToken, tr_val.nToken);
			if(val <= 0)
			{
				LOG(lg_Error,"CurveAnimMaker::loadSettingsFile - Invalid value %f. (must be > 0).\n", val);
				return false;
			}

			if(!pnd) 
				continue;
			
			bool bpos = mstristr(buf,CURVE_POS_SUFFIX) != NULL;
			if(mstristr(buf,CURVE_REC_SUFFIX) || bdef)
				setNodeProps<u32,&NodeProps::nDegree>(bpos ? &aPosProps[0] : &aRotProps[0],val,pnd);
			else
				bpos ? aPosProps[inode].nDegree = val : aRotProps[inode].nDegree = val;
		}
		else if(check_prefix_i_n(tr_op.pToken,CURVE_FLAGS_TOKEN,tr_op.nToken))
		{
			// Get full op name.
			tr_op.translateStringToken(buf,sizeof(buf));
			if(!bt.getStringToken(&tr_val))
			{
				LOG(lg_Error,"CurveAnimMaker::loadSettingsFile failed due to malformed token \"%s\".\n",buf);
				return false;
			}

			if(!pnd) 
				continue;

			char flags_buf[1024];
			tr_val.translateStringToken(flags_buf,sizeof(flags_buf));

			const named_enum flag_table[] = 
			{
				ADD_NAMED_ENUM(cmf_Constant),
				ADD_NAMED_ENUM(cmf_AllowCuts),
				ADD_NAMED_ENUM(cmf_AllowKnotReduction),
				ADD_NAMED_ENUM(cmf_ConstrainEndPos),
				ADD_NAMED_ENUM(cmf_ConstrainEndVel),
				ADD_NAMED_ENUM(cmf_ForceMatchingEndVel),
				ADD_NAMED_ENUM(cmf_DoublePrecisionSolver)
			};

			const int n = sizeof(flag_table) / sizeof(named_enum);
			COMPILE_TIME_ASSERT(cmf_NumMakerFlags == n);
			
			u32 flags = 0;
			BufferTokenizer bt_flags(flags_buf,strlen(flags_buf));
			for(;;)
			{
				TokenRec tr_flag;
				if(!bt_flags.getCToken(&tr_flag))
					break;
				
				int i;
				for(i = 0; i < n; i++)
				{
					if(compare_string_token(tr_flag,flag_table[i].pName))
					{
						flags |= flag_table[i].nEnum;
						break;
					}
				}

				if(i == n)
				{
					tr_flag.translateStringToken(buf,sizeof(buf));
					LOG(lg_Warning,"CurveAnimMaker::loadSettingsFile encountered unrecognized flag \"%s\".\n",buf);
				}
			}

			bool bpos = mstristr(buf,CURVE_POS_SUFFIX) != NULL;
			if(mstristr(buf,CURVE_REC_SUFFIX) || bdef)
				setNodeProps<u32,&NodeProps::nFlags>(bpos ? &aPosProps[0] : &aRotProps[0],flags,pnd);
			else
				bpos ? aPosProps[inode].nFlags = flags : aRotProps[inode].nFlags = flags;
		}
		else if(check_prefix_i_n(tr_op.pToken,CURVE_TYPES_TOKEN,tr_op.nToken))
		{
			tr_op.translateStringToken(buf,sizeof(buf));
			if(!bt.getStringToken(&tr_val))
			{
				LOG(lg_Error,"CurveAnimMaker::loadSettingsFile failed due to malformed token \"%s\".\n",buf);
				return false;
			}

			if(!pnd) 
				continue;

			const named_enum type_table[] = 
			{
				ADD_NAMED_ENUM(ct_Dv_F32),
				ADD_NAMED_ENUM(ct_D3_U16),
				ADD_NAMED_ENUM(ct_D3_U8),
				ADD_NAMED_ENUM(ct_D4n_U16),
				ADD_NAMED_ENUM(ct_D4n_U8),
				ADD_NAMED_ENUM(ct_D4_U16),
				ADD_NAMED_ENUM(ct_D4_U8),
				ADD_NAMED_ENUM(ct_D3_ConstF32),
				ADD_NAMED_ENUM(ct_D4_ConstF32)
			};

			const int n = sizeof(type_table) / sizeof(named_enum);
			COMPILE_TIME_ASSERT(ct_LastType + 1 == n);
			
			char types_buf[1024];
			tr_val.translateStringToken(types_buf,sizeof(types_buf));
			BufferTokenizer bt_flags(types_buf,strlen(types_buf));
			
			int types[ct_MaxTypes];
			int ntypes = 0;

			for(;;)
			{
				TokenRec tr_type;
				if(!bt_flags.getCToken(&tr_type))
					break;
				
				int i;
				for(i = 0; i < n; i++)
				{
					if(compare_string_token(tr_type,type_table[i].pName))
					{
						types[ntypes++] = type_table[i].nEnum;
						break;
					}
				}

				if(i == n)
				{
					tr_type.translateStringToken(buf,sizeof(buf));
					LOG(lg_Warning,"CurveAnimMaker::loadSettingsFile encountered unrecognized type \"%s\".\n",buf);
				}
			}

			bool bpos = mstristr(buf,CURVE_POS_SUFFIX) != NULL;
			setNodeTypes(bpos ? &aPosProps[0] : &aRotProps[0],types,ntypes,mstristr(buf,CURVE_REC_SUFFIX) || bdef,pnd);
		}
	}

	return true;
}
///////////////////////////////////////////////////////////////////////////////
bool CurveAnimMaker::createPosCurve(
						AssetFileResourcePtr       passet,
						AssetFileNodePtr           pnode,
						const NodeData *           pnode_data,
						int                        inode,
						MVector<float,POS_COUNT> & pos_samples,
						const MakeAnimSettings &   mas
						)
{
	float time_scale = pnode->getPositionEnd() - pnode->getPositionStart();
	pos_samples.resize(0);

	Point3 start, delta;
	for(int i = 0, n = pnode->getPositionSegmentCount(); i < n; i++)
	{
		float fpos = ((time_scale * float(i)) / float(n)) + pnode->getPositionStart();
		pnode->getPositionSegment(i, &start, &delta);

		if(mas.pBaseAnim)
		{
			QuatXForm xf_base, xf_final;
			mas.pBaseAnim->calcNTP(pnode_data->idxNode,&xf_base,min(mas.fMaxBaseAnimPos,fpos));
			pnode->calcXForm(&xf_final,fpos);
			xf_base.quatPart.mulInverse(xf_base.transPart,&delta);
			xf_final.quatPart.mul(delta,&xf_base.transPart);
			start -= xf_base.transPart;
		}

		pos_samples.push_back(start.getX());
		pos_samples.push_back(start.getY());
		pos_samples.push_back(start.getZ());
	}

	// Add last translational sample.
	pnode->getPositionSegment(pnode->getPositionSegmentCount()-1, &start, &delta);
	start += delta;

	if(mas.pBaseAnim)
	{
		QuatXForm qxf;
		mas.pBaseAnim->calcNTP(pnode_data->idxNode,&qxf,min(mas.fMaxBaseAnimPos,1.f));
		qxf.quatPart.mulInverse(qxf.transPart,&delta);

		QuatArc qa;
		pnode->getRotationSegment(pnode->getRotationSegmentCount() - 1,&qa);
		qa.calcArc(1.f,&qxf.quatPart);
		qxf.quatPart.mul(delta,&qxf.transPart);
		start -= qxf.transPart;
	}

	pos_samples.push_back(start.getX());
	pos_samples.push_back(start.getY());
	pos_samples.push_back(start.getZ());

#if defined(_DEBUG)
	if(mas.pBaseAnim)
	{
		// Verify Bt * B-1 * Fr + Dt = Ft for all translation keys.
		for(int i = 0, n = pos_samples.size() / 3; i < n; i++)
		{
			float t = ((time_scale * float(i)) / float(n - 1)) + pnode->getPositionStart();

			QuatXForm xf_base, xf_final;
			mas.pBaseAnim->calcNTP(pnode_data->idxNode,&xf_base,min(mas.fMaxBaseAnimPos,t));
			pnode->calcXForm(&xf_final,t);
			
			s32 is = i * 3;
			Point3 trans, delta;
			MulQuatInvQuat(xf_base.quatPart,xf_final.quatPart).mul(xf_base.transPart,&trans);
			delta.set(pos_samples[is+0],pos_samples[is+1],pos_samples[is+2]);
			trans += delta;

			if(!trans.isEqual(xf_final.transPart))
			{
				String node_name, asset_name;
				pnode->getName(&node_name);
				passet->getName(&asset_name);
				LOG(lg_Warning,"\"%s\" produced incorrect translation delta for node \"%s\".\n",(const char *)asset_name,(const char *)node_name);
			}
		}
	}
#endif

	bool bhit_tol = false;
	u32 curve_size;
	CurveData * pcd = curveMaker.makeCompressedCurve(
						&pos_samples[0],
						pos_samples.size() / 3, 
						3,
						aPosProps[inode].nDegree,
						aPosProps[inode].nFlags,
						0,
						aPosProps[inode].fErrTol,
						aPosProps[inode].fCutTol,
						&bhit_tol,
						&curve_size,
						aPosProps[inode].aCurveTypes.size() ? &aPosProps[inode].aCurveTypes[0] : NULL,
						aPosProps[inode].aCurveTypes.size()
						);
	if(!pcd)
	{
		String name;
		passet->getName(&name);
		LOG
		(
			lg_Error,
			"Node \"%s\" for animation \"%s\" failed to create translation curve.\n",
			pnode_data->pszName,name
		);
		return false;
	}

	const CurveType * ptype = CurveType::get(*pcd);
	ASSERT(ptype && ptype->getSize(*pcd) == curve_size);

	if(!bhit_tol)
	{
		String name;
		passet->getName(&name);

		CurveError ce;
		u32 nsamples = (u32)(pos_samples.size() / 3);
		ptype->calcError(*pcd,&pos_samples[0],nsamples,0.f,float(nsamples - 1),0,&ce); 
		LOG
		(
			lg_Warning,
			"Node \"%s\" for animation \"%s\" achieved linear error = %g, but tolerance = %g\n",
			pnode_data->pszName,(const char *)name,sqrtf(ce.fMaxErrSq),aPosProps[inode].fErrTol
		);
	}

	// This bit indicates that translation data was exported for this node.
	ASSERT(inode == pnode_data->idxNode);
	SetBitArrayBit(&aNodeUsage[0],inode * 4 + 1);
	ASSERT(CurveAnim::getNodeUsage(&aNodeUsage[0],inode) & nu_HasTransEntry);
	
	bool bref_curve = false;
	if((ptype->getFlags() & caf_Constant) || (curveMaker.getMakerFlags() & cmf_Constant))
	{
		float p[3];
		p[0] = pnode_data->qxfRefNTP.transPart.getX();
		p[1] = pnode_data->qxfRefNTP.transPart.getY();
		p[2] = pnode_data->qxfRefNTP.transPart.getZ();

		CurveEvalBuf buf;
		ptype->evaluate(*pcd,0.f,0,&buf);
		if(CurveCalcErrorSq(buf.aCurve,p,3,0) <= square(aPosProps[inode].fErrTol))
		{
			// Turns out curve was just the reference pose.
			bref_curve = true;
		}
	}

	if(bref_curve)
	{
		// Curve is the reference pose.
		SetBitArrayBit(&aNodeUsage[0],inode * 4 + 3);
		ASSERT(CurveAnim::getNodeUsage(&aNodeUsage[0],inode) & nu_RefTransEntry);
		nRefCount++;
	}
	else
	{
		// We're actually going to emit this curve.
		aNodeIds.push_back(inode);
		
		// Allocate space and copy finished curve.
		aCurveData.increment().setData(new byte[curve_size]);
		memcpy(&aCurveData.last()[0],pcd,curve_size);
		pcd = (CurveData *)&aCurveData.last()[0];
		aCurves.push_back(pcd);

		nTotalCurveSize += curve_size;
		ASSERT(CurveType::get(*pcd)->getSize(*pcd) == curve_size);
	}

	return true;
}
///////////////////////////////////////////////////////////////////////////////
bool CurveAnimMaker::createRotCurve(
						AssetFileResourcePtr       passet,
						AssetFileNodePtr           pnode,
						const NodeData *           pnode_data,
						int                        inode,
						MVector<float,ROT_COUNT> & rot_samples,
						const MakeAnimSettings &   mas
						)
{
	float time_scale = pnode->getPositionEnd() - pnode->getPositionStart();
	rot_samples.resize(0);

	QuatArc qa;
	for(int i = 0, n = pnode->getRotationSegmentCount(); i < n; i++)
	{
		pnode->getRotationSegment(i, &qa);
		float fpos = ((time_scale * float(i)) / float(n)) + pnode->getPositionStart();
		
		if(mas.pBaseAnim)
		{
			Quat qbase;
			mas.pBaseAnim->calcNTPQuat(pnode_data->idxNode,&qbase,min(mas.fMaxBaseAnimPos,fpos));
			qa.u = MulQuatInvQuat(qbase,qa.u).normalize();
		}
		
		rot_samples.push_back(qa.u.pt.getX());
		rot_samples.push_back(qa.u.pt.getY());
		rot_samples.push_back(qa.u.pt.getZ());
		rot_samples.push_back(qa.u.pt.getW());
	}

	pnode->getRotationSegment(pnode->getRotationSegmentCount() - 1, &qa);
	
	Quat qlast;
	qa.calcArc(1,&qlast);
	if(mas.pBaseAnim)
	{
		Quat qbase;
		mas.pBaseAnim->calcNTPQuat(pnode_data->idxNode,&qbase,min(mas.fMaxBaseAnimPos,1.f));
		qlast = MulQuatInvQuat(qbase,qlast).normalize();
	}

	rot_samples.push_back(qlast.pt.getX());
	rot_samples.push_back(qlast.pt.getY());
	rot_samples.push_back(qlast.pt.getZ());
	rot_samples.push_back(qlast.pt.getW());

#if defined(_DEBUG)
	if(mas.pBaseAnim)
	{
		// Verify base * delta = sample for all rotation samples.
		for(int i = 0, n = rot_samples.size() / 4; i < n; i++)
		{
			float t = ((time_scale * float(i)) / float(n - 1)) + pnode->getPositionStart();
			QuatXForm xf_final;
			pnode->calcXForm(&xf_final,t);
			
			Quat qrot;
			mas.pBaseAnim->calcNTPQuat(pnode_data->idxNode,&qrot,min(mas.fMaxBaseAnimPos,t));

			s32 is = i * 4;
			Quat qdelta;
			qdelta.pt.set(rot_samples[is+0],rot_samples[is+1],rot_samples[is+2],rot_samples[is+3]);
			qrot = MulQuatQuat(qrot,qdelta);

			if(!qrot.isEqual(xf_final.quatPart))
			{
				String node_name, asset_name;
				pnode->getName(&node_name);
				passet->getName(&asset_name);
				LOG(lg_Warning,"\"%s\" produced incorrect rotational delta for node \"%s\".\n",(const char *)asset_name,(const char *)node_name);
			}
		}
	}
#endif

	bool bhit_tol = false;
	u32 curve_size;
	CurveData * pcd = curveMaker.makeCompressedCurve(
						&rot_samples[0],
						rot_samples.size() / 4, 
						4,
						aRotProps[inode].nDegree,
						aRotProps[inode].nFlags,
						cef_Quaternion,
						aRotProps[inode].fErrTol,
						aRotProps[inode].fCutTol,
						&bhit_tol,
						&curve_size,
						aRotProps[inode].aCurveTypes.size() ? &aRotProps[inode].aCurveTypes[0] : NULL,
						aRotProps[inode].aCurveTypes.size()
						);
	if(!pcd)
	{
		String name;
		passet->getName(&name);
		LOG
		(
			lg_Error,
			"Node \"%s\" for animation \"%s\" failed to create rotation curve.\n",
			pnode_data->pszName,(const char *)name
		);
		return false;
	}

	const CurveType * ptype = CurveType::get(*pcd);
	ASSERT(ptype && ptype->getSize(*pcd) == curve_size);

	if(!bhit_tol)
	{
		String name;
		passet->getName(&name);

		CurveError ce;
		u32 nsamples = (u32)(rot_samples.size() / 4);
		ptype->calcError(*pcd,&rot_samples[0],nsamples,0.f,float(nsamples - 1),cef_Quaternion,&ce); 
		LOG
		(
			lg_Warning,
			"Node \"%s\" for animation \"%s\" achieved angular error of %g, but tolerance = %g\n",
			pnode_data->pszName,(const char *)name,sqrtf(ce.fMaxErrSq),aRotProps[inode].fErrTol
		);
	}

	// Set bit indicating that this node has rotational data.
	ASSERT(inode == pnode_data->idxNode);
	SetBitArrayBit(&aNodeUsage[0],inode * 4);
	ASSERT(CurveAnim::getNodeUsage(&aNodeUsage[0],inode) & nu_HasRotEntry);

	bool bref_curve = false;
	if((ptype->getFlags() & caf_Constant) || (curveMaker.getMakerFlags() & cmf_Constant))
	{
		float q[4];
		q[0] = pnode_data->qxfRefNTP.quatPart.pt.getX();
		q[1] = pnode_data->qxfRefNTP.quatPart.pt.getY();
		q[2] = pnode_data->qxfRefNTP.quatPart.pt.getZ();
		q[3] = pnode_data->qxfRefNTP.quatPart.pt.getW();

		CurveEvalBuf buf;
		ptype->evaluate(*pcd,0.f,cef_Quaternion,&buf);
		if(CurveCalcErrorSq(buf.aCurve,q,4,cef_Quaternion) <= square(aRotProps[inode].fErrTol))
		{
			// Turns out curve was just the reference pose.
			bref_curve = true;
		}
	}

	if(bref_curve)
	{
		// Curve is the reference pose.
		SetBitArrayBit(&aNodeUsage[0],inode * 4 + 2);
		ASSERT(CurveAnim::getNodeUsage(&aNodeUsage[0],inode) & nu_RefRotEntry);
		nRefCount++;
	}
	else
	{
		// We're actually going to emit this curve.
		aNodeIds.push_back(cac_QuatMask | inode);

		// Allocate space and copy finished curve.
		aCurveData.increment().setData(new byte[curve_size]);
		memcpy(&aCurveData.last()[0],pcd,curve_size);
		pcd = (CurveData *)&aCurveData.last()[0];
		aCurves.push_back(pcd);

		nTotalCurveSize += curve_size;
		ASSERT(CurveType::get(*pcd)->getSize(*pcd) == curve_size);
	}

	return true;
}
///////////////////////////////////////////////////////////////////////////////
void CurveAnimMaker::sortCurves()
{
	if(aCurves.size() <= 1)
		return;

	MVector<CurveSortRec> sort_recs;
	sort_recs.resize(aCurves.size());

	int i,n = aCurves.size();
	for(i = 0; i < n; i++)
	{
		sort_recs[i].init(aCurves[i],aNodeIds[i]);
	}

	quick_sort<CurveSortRec,CurveSortRecComparator>(&sort_recs[0],sort_recs.size());

	for(i = 0; i < n; i++)
	{
		aCurves[i]  = sort_recs[i].pCurve;
		aNodeIds[i] = sort_recs[i].getNodeId();
	}

#if defined(_DEBUG)
	static bool bdump_samplers = false;
	if(bdump_samplers)
	{
		for(i = 0; i < n; i++)
		{
			const CurveData & cd = *sort_recs[i].pCurve;
			u32 inode = sort_recs[i].getNodeId();
			u32 flags = (inode & cac_QuatMask) ? cef_Quaternion : 0;
			CurveSampler ps = CurveType::get(cd)->getSampler(cd,flags);
			ASSERT(ps);
			LOG(lg_Anim,"Sampler for node %d = 0x%Xh.\n",inode & cac_NodeMask,ps);
		}
	}
#endif
}
///////////////////////////////////////////////////////////////////////////////
bool CurveAnimMaker::init(AssetFileResourcePtr passet,const MakeAnimSettings & mas)
{
	if(!loadSettingsFile(mas))
		return false;

	u32 ndw = (nNodes * 4 + 31) >> 5;
	aNodeUsage.setData(new u32[ndw]);
	nNodeUsageSize = sizeof(u32) * ndw;
	memset(&aNodeUsage[0],0,nNodeUsageSize);

	// Probably we won't have this many curves, this is just the max.
	aCurves.reserve(nNodes * 2);
	aNodeIds.reserve(nNodes * 2);
	aCurveData.reserve(nNodes * 2);

	nRefCount = 0;
	nTotalCurveSize = 0;

	MVector<float,POS_COUNT> pos_samps;
	MVector<float,ROT_COUNT> rot_samps;
	pos_samps.reserve(POS_COUNT);
	rot_samps.reserve(ROT_COUNT);

	const NodeData ** ppnd = mas.pNodeTree->getNodeOrders();
	
	// First, create all position curves.
	int i;
	for(i = 0; i < nNodes; i++)
	{
		AssetFileNodePtr pnode = passet->getNode(ppnd[i]->pszName);
		if(!pnode || passet->getNodeTree()->isRootNode(pnode))
			continue;

		if(!createPosCurve(passet,pnode,ppnd[i],i,pos_samps,mas))
			return false;
	}

	// Then create all rotation curves.
	for(i = 0; i < nNodes; i++)
	{
		AssetFileNodePtr pnode = passet->getNode(ppnd[i]->pszName);
		if(!pnode || passet->getNodeTree()->isRootNode(pnode))
			continue;	

		if(!createRotCurve(passet,pnode,ppnd[i],i,rot_samps,mas))
			return false;
	}

	u32 type_sizes[ct_MaxTypes];
	u32 type_count[ct_MaxTypes];
	memset(type_sizes,0,sizeof(type_sizes));
	memset(type_count,0,sizeof(type_count));
	
	for(i = 0; i < aCurves.size(); i++)
	{
		const CurveData & cd = *aCurves[i];
		type_sizes[cd.nType] += CurveType::get(cd)->getSize(cd);
		type_count[cd.nType]++;
	}

	for(i = ct_FirstType; i <= ct_LastType; i++)
	{
		LOG(lg_Anim,"Type \"%s\": Count = %d, Size = %d.\n",CurveType::get(i)->getName(),type_count[i],type_sizes[i]);
	}

	// All curves created, now sort them by sampler.
	sortCurves();
	return true;
}
///////////////////////////////////////////////////////////////////////////////

#endif