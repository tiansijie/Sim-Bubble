//
// Copyright (C) CIS660
// 
// File: BubbleSimMayaNode.cpp
//
// Dependency Graph Node: BubbleSimMaya
//
// Author: Maya Plug-in Wizard 2.0.1
// Revised by Isaac Peral
//

#include "BubbleEmitterMayaNode.h"

#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>

#include <maya/MGlobal.h>

// You MUST change this to a unique value!!!  The id is a 32bit value used
// to identify this type of node in the binary file format.  
//
MTypeId     BubbleEmitterMaya::id( 0x00007 );

MObject     BubbleEmitterMaya::input;        
MObject     BubbleEmitterMaya::output;       

BubbleEmitterMaya::BubbleEmitterMaya() {}
BubbleEmitterMaya::~BubbleEmitterMaya() {}

MStatus BubbleEmitterMaya::compute( const MPlug& plug, MDataBlock& data )
//
//	Description:
//		This method computes the value of the given output plug based
//		on the values of the input attributes.
//
//	Arguments:
//		plug - the plug to compute
//		data - object that provides access to the attributes for this node
//
{
	MStatus returnStatus;
 
	if( plug == output )
	{
		MDataHandle inputData = data.inputValue( input, &returnStatus );

		if( returnStatus != MS::kSuccess )
			MGlobal::displayError( "Node BubbleSimMaya cannot get value\n" );
		else
		{
			float result = inputData.asFloat();

			MDataHandle outputHandle = data.outputValue( BubbleEmitterMaya::output );
			outputHandle.set( result );
			data.setClean(plug);
		}
	} else {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}

void* BubbleEmitterMaya::creator()
//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
{
	return new BubbleEmitterMaya();
}

MStatus BubbleEmitterMaya::initialize()
//
//	Description:
//		This method is called to create and initialize all of the attributes
//      and attribute dependencies for this node type.  This is only called 
//		once when the node type is registered with Maya.
//
//	Return Values:
//		MS::kSuccess
//		MS::kFailure
//		
{
	MFnNumericAttribute nAttr;
	MStatus				stat;

	input = nAttr.create( "input", "in", MFnNumericData::kFloat, 0.0 );
 	nAttr.setStorable(true);
 	nAttr.setKeyable(true);

	output = nAttr.create( "output", "out", MFnNumericData::kFloat, 0.0 );
	nAttr.setWritable(false);
	nAttr.setStorable(false);

	stat = addAttribute( input );
		if (!stat) { stat.perror("addAttribute"); return stat;}
	stat = addAttribute( output );
		if (!stat) { stat.perror("addAttribute"); return stat;}

	stat = attributeAffects( input, output );
		if (!stat) { stat.perror("attributeAffects"); return stat;}

	return MS::kSuccess;

}

