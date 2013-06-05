#ifndef CreateBubbleSimMayaNode_H_
#define CreateBubbleSimMayaNode_H_

#include <maya/MIOStream.h>
#include <maya/MVector.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MPxFieldNode.h>
#include "BubbleGeniusCmd.h"

// macros
#define MAKE_INPUT(attr) \
	CHECK_MSTATUS(attr.setKeyable(true)); \
	CHECK_MSTATUS(attr.setStorable(true)); \
	CHECK_MSTATUS(attr.setReadable(true)); \
	CHECK_MSTATUS(attr.setWritable(true));
#define MAKE_OUTPUT(attr) \
	CHECK_MSTATUS(attr.setKeyable(false)); \
	CHECK_MSTATUS(attr.setStorable(false)); \
	CHECK_MSTATUS(attr.setReadable(true)); \
	CHECK_MSTATUS(attr.setWritable(false));
#define MAKE_ADDR(attr) \
	CHECK_MSTATUS(attr.setKeyable(false)); \
	CHECK_MSTATUS(attr.setStorable(false)); \
	CHECK_MSTATUS(attr.setReadable(true)); \
	CHECK_MSTATUS(attr.setWritable(false)); \
	CHECK_MSTATUS(attr.setHidden(true));

#define McheckErr(stat, msg)		\
	if ( MS::kSuccess != stat )		\
	{								\
	cerr << msg;				\
	return MS::kFailure;		\
	}

class BubbleSimMaya: public MPxFieldNode
{
public:
	BubbleSimMaya() {};
	virtual ~BubbleSimMaya() {};

	static void		*creator();
	static MStatus	initialize();

	// will compute output force.
	//
	virtual MStatus	compute( const MPlug& plug, MDataBlock& block );

	static MTypeId	id;

	static MObject inNumPoints;

/*
	static MObject inMinimumX;
	static MObject inMinimumY;
	static MObject inMinimumZ;
	static MObject inMinimum;

	static MObject inMaximumX;
	static MObject inMaximumY;
	static MObject inMaximumZ;
	static MObject inMaximum;*/
///////////////////////////////////////////////////////////////
	static MObject eTransX;
	static MObject eTransY;
	static MObject eTransZ;
	static MObject eTrans;

	static MObject eScaleX;
	static MObject eScaleY;
	static MObject eScaleZ;
	static MObject eScale;

	static MObject fResWidth;
	static MObject fResHeight;
	static MObject fResDepth;
	static MObject fResolution;

	static MObject  time;

	static MObject boundingBoxMinX;
	static MObject boundingBoxMinY;
	static MObject boundingBoxMinZ;
	static MObject boundingBoxMin;

	static MObject boundingBoxMaxX;
	static MObject boundingBoxMaxY;
	static MObject boundingBoxMaxZ;
	static MObject boundingBoxMax;

	static MObject boundingBoxSizeX;
	static MObject boundingBoxSizeY;
	static MObject boundingBoxSizeZ;
	static MObject boundingBoxSize;

	static MObject baseResolution;
//////////////////////////////////////////////////////////////
	static MObject outPoints;

	static MObject positionX;
	static MObject positionY;
	static MObject positionZ;
	static MObject stiffCoeff;
	static MObject weakCoeff;
	static MObject dampCoeffV;
	static MObject dampCoeffL;
	static MObject solidAdhesionCoeff;
	static MObject solidAttractionCoeff;
	static MObject burstingCoeff;
	static MObject wetness;
	static MObject liquidAdhesion;
	static MObject dragCoeff;

	static MObject emitterPositionX;
	static MObject emitterPositionY;
	static MObject emitterPositionZ;
	static MObject emitterPosition;

	static MObject bubbleNumber;
	static MObject bubbleSize;
};

#endif