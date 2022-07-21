#ifndef __CURVE_ANIM_H__
#define __CURVE_ANIM_H__

#include "NodeTreeAnim.h"
#include "Curve.h"
#include "CurveMaker.h"

BEGIN_PMATH_NAMESPACE

#if defined(PWTOOL)
class CurveAnimMaker;
#endif

DECLARE_LOG_CATEGORY(lg_Anim)

enum CurveAnimConsts
{
	cac_NodeMask = (1 << 15) - 1,
	cac_QuatMask = 1 << 15
};

///////////////////////////////////////////////////////////////////////////////
class CurveAnim : public NodeTreeAnim
{
protected:

	struct InstHeader
	{
		u8  nPosRef;
		u8  nRotRef;
		u16 nSamplers;
		u32 nTotSize;
	};

	u16 nCurves;
	u16 nNodes;
	u32 nMaxInstSize;

	union
	{
		u32 nCurvesOffset;
		const CurveData ** ppCurves;
	};

	// Translates curve indices into node indices.
	union
	{
		u32 nNodeIdsOffset;
		const u16 * pNodeIds;
	};

	union
	{
		u32 nNodeUsageOffset;
		const u32 * pNodeUsage;
	};

	u32 getNodeUsage(s32 inode) const;
	s32 getCurveIndex(u32 inode) const;

public:

    CurveAnim() {}
    virtual ~CurveAnim() {}

	virtual void  calcNTP(s32 idx_node, QuatXForm * pxf, float fpos) const;
	virtual void  calcNTPQuat(s32 idx_node, Quat * pquat, float fpos) const;
	virtual void  calcNTPTrans(s32 idx_node, Point3 * ppt, float fpos) const;
	
	virtual u32   nodeUsage(int idx_node) const;

	virtual int   getMaxInstSize() const;
	virtual int   getInstSize(int nnodes, const u32 * pnode_usage) const;
	virtual void  initInst(byte * pb, int nmax_nodes, const u32 * pnode_usage) const;
	virtual void  shutdownInst(byte * pb) const;
	virtual void  calcNTPInst(QuatXForm * pxf, float fpos, byte * pinst_data) const;
	
	virtual void  fixup(NodeTree *pnt);		

	static CurveAnim * initCurveAnimFromImage(byte *pdata, NodeTree *pnt);
	static u32 getNodeUsage(const u32 * pnode_usage, s32 inode);

#if defined(PWTOOL)
	static CurveAnim * makeCurveAnimFromAssetFile(AssetFileResourcePtr passet,const CurveAnimMaker & cam, NodeTree * pnt);
	static void fixEndian(byte * pdata);
	virtual void finalWritePrep(u32 id_node_tree); 
#endif
};
///////////////////////////////////////////////////////////////////////////////
inline u32 CurveAnim::getNodeUsage(s32 idx_node) const
{
	ASSERT_PRINTF(idx_node >= 0 && idx_node < pNodeTree->noNodes(), 
		"Animation \"%s\" uses a skeleton which only has %d bones.\n",debName(),pNodeTree->noNodes());
	
#if defined(PWTOOL)
	if(idx_node < 0 || idx_node >= pNodeTree->noNodes())
	{
		return nu_Unused;
	}
#endif

	return getNodeUsage(pNodeUsage,idx_node);
}
///////////////////////////////////////////////////////////////////////////////
inline u32 CurveAnim::getNodeUsage(const u32 * pnode_usage, s32 inode)
{
	u32 istart_bit = inode << 2;
	u32 idx = istart_bit >> 5;
	u32 dw = pnode_usage[idx];

	u32 offset = istart_bit - (idx << 5);

	dw >>= offset;
	dw = dw & nu_UsageMask;
	return dw;
}
///////////////////////////////////////////////////////////////////////////////
#if defined(PWTOOL)
class CurveAnimMaker
{
protected:

	typedef MVector<const CurveType *> CurveTypeArray;
	typedef MVectorCD< AutoArrayPtr<byte> > CurveDataArray;

	enum
	{
		POS_COUNT  = 2   * 1024,
		ROT_COUNT  = 4   * 1024
	};

	struct NodeProps
	{
		float          fErrTol;
		float          fCutTol;
		u32            nDegree;
		u32            nFlags;
		CurveTypeArray aCurveTypes;

		NodeProps()
		{
			fErrTol = 0.f;
			fCutTol = 0.f;
			nDegree = 3;
			nFlags  = 0;
		};
	};

	AutoArrayPtr<NodeProps> aPosProps;
	AutoArrayPtr<NodeProps> aRotProps;

	MVector<u32>            aNodeIds;       // List of nodes that actually need a curve.
	AutoArrayPtr<u32>       aNodeUsage;     // Bitfield indicating how each node is going to be used.
	CurveDataArray          aCurveData;     // MVector of arrays, one for each finished curve.
	MVector<CurveData *>    aCurves;        // Pointers to finished curves, memory is owned by aCurveData.  
	CurveMaker              curveMaker;     // This is more memory efficient than having one per node.

	u32                     nTotalCurveSize;
	u32                     nNodeUsageSize;
	u32                     nRefCount;
	int                     nNodes;

	bool loadSettingsFile(const MakeAnimSettings & mas);

	bool createRotCurve(
			AssetFileResourcePtr       passet,
			AssetFileNodePtr           pnode,
			const NodeData *           pnode_data,
			int                        inode,
			MVector<float,ROT_COUNT> & rot_samps,
			const MakeAnimSettings &   mas
			);

	bool createPosCurve(
			AssetFileResourcePtr       passet,
			AssetFileNodePtr           pnode,
			const NodeData *           pnode_data,
			int                        inode,
			MVector<float,POS_COUNT> & pos_samps,
			const MakeAnimSettings &   mas
			);

	void sortCurves();

	template<typename V, V NodeProps::* M> 
	void setNodeProps(NodeProps * pprops, V val, const NodeData * pnd);
	void setNodeTypes(NodeProps * pprops,int * ptypes,int ntypes,bool brec,const NodeData * pnd);
	
public:

	bool init(AssetFileResourcePtr passet,const MakeAnimSettings & settings);

	u32 getCurveSize()         const { return nTotalCurveSize; }
	u32 getCurveCount()        const { return aCurves.size(); }
	u32 getRefCount()          const { return nRefCount; }
	u32 getNodeId(s32 icurve)  const { return aNodeIds[icurve]; }
	u32 getNodeUsageSize()     const { return nNodeUsageSize; }

	const u32 * getNodeUsage() const { return &aNodeUsage[0]; }
	const CurveData * const * getCurves() const { return aCurves.size() ? &aCurves[0] : NULL; }
};
///////////////////////////////////////////////////////////////////////////////
template<typename V,V CurveAnimMaker::NodeProps::* M> 
void CurveAnimMaker::setNodeProps(NodeProps * pprops, V val, const NodeData * pnd)
{
	ASSERT(pnd);
	s32 inode = pnd->idxNode;
	ASSERT(inode >= 0 && inode < nNodes);
	tset_member<NodeProps,V,M>(pprops + inode,val);
	
	for(int i = 0; i < pnd->nChildren; i++)
		setNodeProps<V,M>(pprops,val,pnd->pChildren[i]);
}
///////////////////////////////////////////////////////////////////////////////
#endif

END_PMATH_NAMESPACE

#endif // __CURVE_ANIM_H__