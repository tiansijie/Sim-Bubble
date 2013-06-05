#ifndef _BubbleEmitterMayaNode
#define _BubbleEmitterMayaNode

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

#define McheckErr(stat,msg)                     \
        if ( MS::kSuccess != stat ) {   \
                cerr << msg;                            \
                return MS::kFailure;            \
        }

//
// Copyright (C) CIS660
// 
// File: BubbleSimMayaNode.h
//
// Dependency Graph Node: BubbleSimMaya
//
// Author: Maya Plug-in Wizard 2.0.1
// Revised by Isaac Peral
//

#include <maya/MPxNode.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MTypeId.h> 

 
class BubbleEmitterMaya : public MPxNode
{
public:
						BubbleEmitterMaya();
	virtual				~BubbleEmitterMaya(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();

public:

	static  MObject		input;		// Example input attribute
	static  MObject		output;		// Example output attribute


	static	MTypeId		id;
};

#endif
