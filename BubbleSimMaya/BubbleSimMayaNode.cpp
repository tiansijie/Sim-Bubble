#include <maya/MIOStream.h>
#include <math.h>
#include <time.h>    


#include <BubbleSimMayaNode.h>

#include <maya/MTime.h>
#include <maya/MVectorArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MMatrix.h>
#include <maya/MGlobal.h>
#include <maya/MPlane.h>
#include <maya/MSelectionList.h>
#include <maya/MFnMesh.h>
#include <maya/MDagPath.h>

#include <maya/MFnDependencyNode.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnArrayAttrsData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFloatVectorArray.h>  
#include <maya/MArrayDataHandle.h>

#define MNoVersionString
#define MNoPluginEntry
#include <maya/MFnPlugin.h>

#include <../bubble/bubble_sim.h>
#include <../bubble/bubble.h>

#include "makelevelset3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

#define SimScale 1.0f

MTypeId BubbleSimMaya::id( 0x9 );


MObject BubbleSimMaya::inNumPoints;

/*
MObject BubbleSimMaya::inMinimumX;
MObject BubbleSimMaya::inMinimumY;
MObject BubbleSimMaya::inMinimumZ;
MObject BubbleSimMaya::inMinimum;

MObject BubbleSimMaya::inMaximumX;
MObject BubbleSimMaya::inMaximumY;
MObject BubbleSimMaya::inMaximumZ;
MObject BubbleSimMaya::inMaximum;*/

MObject BubbleSimMaya::outPoints;
/////////////////////////////////////////////
MObject BubbleSimMaya:: eTransX;
MObject BubbleSimMaya:: eTransY;
MObject BubbleSimMaya:: eTransZ;
MObject BubbleSimMaya:: eTrans;

MObject BubbleSimMaya:: eScaleX;
MObject BubbleSimMaya:: eScaleY;
MObject BubbleSimMaya:: eScaleZ;
MObject BubbleSimMaya:: eScale;

MObject BubbleSimMaya:: fResWidth;
MObject BubbleSimMaya:: fResHeight;
MObject BubbleSimMaya:: fResDepth;
MObject BubbleSimMaya:: fResolution;

MObject BubbleSimMaya:: time;

MObject BubbleSimMaya:: boundingBoxMinX;
MObject BubbleSimMaya:: boundingBoxMinY;
MObject BubbleSimMaya:: boundingBoxMinZ;
MObject BubbleSimMaya:: boundingBoxMin;

MObject BubbleSimMaya:: boundingBoxMaxX;
MObject BubbleSimMaya:: boundingBoxMaxY;
MObject BubbleSimMaya:: boundingBoxMaxZ;
MObject BubbleSimMaya:: boundingBoxMax;

MObject BubbleSimMaya:: boundingBoxSizeX;
MObject BubbleSimMaya:: boundingBoxSizeY;
MObject BubbleSimMaya:: boundingBoxSizeZ;
MObject BubbleSimMaya:: boundingBoxSize;

MObject BubbleSimMaya:: baseResolution;
/////////////////////////////////////////////
 MObject BubbleSimMaya::  positionX;
 MObject BubbleSimMaya::  positionY;
 MObject BubbleSimMaya::  positionZ;
 MObject BubbleSimMaya::  stiffCoeff;
 MObject BubbleSimMaya::  weakCoeff;
 MObject BubbleSimMaya::  dampCoeffV;
 MObject BubbleSimMaya::  dampCoeffL;
 MObject BubbleSimMaya::  solidAdhesionCoeff;
 MObject BubbleSimMaya::  solidAttractionCoeff;
 MObject BubbleSimMaya::  burstingCoeff;
 MObject BubbleSimMaya::  wetness;
 MObject BubbleSimMaya::  liquidAdhesion;
 MObject BubbleSimMaya::  dragCoeff;


 MObject BubbleSimMaya:: emitterPositionX;
 MObject BubbleSimMaya:: emitterPositionY;
 MObject BubbleSimMaya:: emitterPositionZ;
 MObject BubbleSimMaya:: emitterPosition;

 MObject BubbleSimMaya:: bubbleNumber;
 MObject BubbleSimMaya:: bubbleSize;


std::vector<std::pair<glm::vec3, glm::vec3>> solidWall;

void *BubbleSimMaya::creator()
{
	return new BubbleSimMaya;
}

Bubbles_sim myBubble;

MStatus BubbleSimMaya::initialize()
	//
	//	Descriptions:
	//		Initialize the node, attributes.
	//
{
	MStatus status;
	MFnTypedAttribute tAttr;
	MFnNumericAttribute	nAttr;
	MFnUnitAttribute uAttr;


	MSelectionList list;
	status = MGlobal::getActiveSelectionList(list);
	if(!status) {
		cout << "getActiveSelectionList FAILED\n";
		return( status );
	}

	MDagPath path;
	//list.length(status);
	//MObject mm;
	//status = list.getDependNode(0,mm);
	//if(!status) {
	//	cout << "getDependNode FAILED\n";
	//	return( status );
	//}

	/*for(int i = 0; i < list.length(); i++){
		status = list.getDagPath(i, path);
		if(!status) {
			cout << "getDagPath FAILED\n";
			return( status );
		}
		else if(status) {
			cout<<"success"<<endl;
			MFnMesh mplane(path);
			MPointArray planePoints;
			status  = mplane.getPoints(planePoints, MSpace::kWorld);
			if(!status){
				cout<<"getPoint falied"<<endl;
				return(status);
			}
			MFloatVectorArray  planeNormals;
			status = mplane.getNormals(planeNormals, MSpace::kWorld);
			if(!status){
				cout<<"getNormal falied"<<endl;
				return(status);
			}
			glm::vec3 p = glm::vec3(planePoints[0][0], planePoints[0][1], planePoints[0][2]);
			glm::vec3 n = glm::vec3(planeNormals[0][0], planeNormals[0][1], planeNormals[0][2]);
			n = glm::normalize(n);
			myBubble.addSoild(n,p);
			std::pair<glm::vec3, glm::vec3> solid;
			solid.first = p;
			solid.second = n;
			solidWall.push_back(solid);
			cout<<planePoints[0]<<endl;
			cout<<planeNormals[0]<<endl;
		}
	}*/

	
	//MStringArray kk;
	//list.getSelectionStrings(0,kk);
	//cout<<kk[0]<<endl;
	//MPlane m;

	/*myBubble.addSoild(glm::vec3(-1,0,0), glm::vec3(2.5,1,0) / SimScale);
	myBubble.addSoild(glm::vec3(1,0,0), glm::vec3(-2.5,1,0) / SimScale);
	myBubble.addSoild(glm::vec3(0,0,-1), glm::vec3(0,1,5) / SimScale);
	myBubble.addSoild(glm::vec3(0,0,1), glm::vec3(0,1,-5) / SimScale);*/

	time = uAttr.create( "time", "tm", MFnUnitAttribute::kTime, 0.0, &status );
	MAKE_INPUT(uAttr);


	inNumPoints = nAttr.create("numberPts", "np", MFnNumericData::kInt, 10);
	MAKE_INPUT(nAttr)

/*
	inMinimumX = nAttr.create("inMinX", "MiX", MFnNumericData::kFloat, -10.0);
	MAKE_INPUT(nAttr)
	inMinimumY = nAttr.create("inMinY", "MiY", MFnNumericData::kFloat, 0.0);
	MAKE_INPUT(nAttr)
	inMinimumZ = nAttr.create("inMinZ", "MiZ", MFnNumericData::kFloat, -10.0);
	MAKE_INPUT(nAttr)
	inMinimum = nAttr.create("inMin", "Min", inMinimumX, inMinimumY, inMinimumZ);
	MAKE_INPUT(nAttr)

	inMaximumX = nAttr.create("inMaxX", "MaX", MFnNumericData::kFloat, 10.0);
	MAKE_INPUT(nAttr)
	inMaximumY = nAttr.create("inMaxY", "MaY", MFnNumericData::kFloat, 10.0);
	MAKE_INPUT(nAttr)
	inMaximumZ = nAttr.create("inMaxZ", "MaZ", MFnNumericData::kFloat, 10.0);
	MAKE_INPUT(nAttr)
	inMaximum = nAttr.create("inMax", "Max", inMaximumX, inMaximumY, inMaximumZ);
	MAKE_INPUT(nAttr)*/

	outPoints = tAttr.create("outPoints", "op", MFnArrayAttrsData::kDynArrayAttrs);
	MAKE_OUTPUT(tAttr)
////////////////////////////////////////////////////////////////////////////////////
	eTransX = nAttr.create("emitterTransX", "etX", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	eTransY = nAttr.create("emitterTransY", "etY", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	eTransZ = nAttr.create("emitterTransZ", "etZ", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	eTrans = nAttr.create("emitterTrans", "et", eTransX, eTransY, eTransZ);
	MAKE_INPUT(nAttr)

	eScaleX = nAttr.create("emitterScaleX", "esX", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	eScaleY = nAttr.create("emitterScaleY", "esY", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	eScaleZ = nAttr.create("emitterScaleZ", "esZ", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	eScale = nAttr.create("emitterScale", "es", eScaleX, eScaleY, eScaleZ);
	MAKE_INPUT(nAttr)

    fResWidth = nAttr.create("fluidResoWidth", "frW", MFnNumericData::kInt);
	MAKE_INPUT(nAttr)
	fResHeight = nAttr.create("fluidResoHeight", "frH", MFnNumericData::kInt);
	MAKE_INPUT(nAttr)
	fResDepth = nAttr.create("fluidResoDepth", "frD", MFnNumericData::kInt);
	MAKE_INPUT(nAttr)
    fResolution = nAttr.create("fluidResolution", "fr", fResWidth, fResHeight, fResDepth);
	MAKE_INPUT(nAttr)

	boundingBoxMinX = nAttr.create("boundBoxMinX", "minX", MFnNumericData::kFloat, -10.0f);
	MAKE_INPUT(nAttr)
	boundingBoxMinY = nAttr.create("boundBoxMinY", "minY", MFnNumericData::kFloat, -10.0f);
	MAKE_INPUT(nAttr)
	boundingBoxMinZ = nAttr.create("boundBoxMinZ", "minZ", MFnNumericData::kFloat, -10.0f);
	MAKE_INPUT(nAttr)
	boundingBoxMin = nAttr.create("boundBoxMin", "bmin", boundingBoxMinX, boundingBoxMinY, boundingBoxMinZ);
	MAKE_INPUT(nAttr)

	boundingBoxMaxX = nAttr.create("boundBoxMaxX", "maxX", MFnNumericData::kFloat, 10.0f);
	MAKE_INPUT(nAttr)
	boundingBoxMaxY = nAttr.create("boundBoxMaxY", "maxY", MFnNumericData::kFloat, 10.0f);
	MAKE_INPUT(nAttr)
	boundingBoxMaxZ = nAttr.create("boundBoxMaxZ", "maxZ", MFnNumericData::kFloat, 10.0f);
	MAKE_INPUT(nAttr)
	boundingBoxMax = nAttr.create("boundBoxMax", "bmax", boundingBoxMaxX, boundingBoxMaxY, boundingBoxMaxZ);
	MAKE_INPUT(nAttr)

	boundingBoxSizeX = nAttr.create("boundBoxSizeX", "szX", MFnNumericData::kFloat, 20.0f);
	MAKE_INPUT(nAttr)
	boundingBoxSizeY = nAttr.create("boundBoxSizeY", "szY", MFnNumericData::kFloat, 20.0f);
	MAKE_INPUT(nAttr)
	boundingBoxSizeZ = nAttr.create("boundBoxSizeZ", "szZ", MFnNumericData::kFloat, 20.0f);
	MAKE_INPUT(nAttr)
	boundingBoxSize = nAttr.create("boundBoxSize", "bsz", boundingBoxSizeX, boundingBoxSizeY, boundingBoxSizeZ);
	MAKE_INPUT(nAttr)

	baseResolution = nAttr.create("baseReso", "bRes", MFnNumericData::kInt, 200);
	MAKE_INPUT(nAttr)
///////////////////////////////////////////////////////////////////////////////////////////


	////bubble attribute
	positionX = nAttr.create("positionX", "px", MFnNumericData::kFloat, 0.0f);
	MAKE_INPUT(nAttr)
	positionY = nAttr.create("positionY", "py", MFnNumericData::kFloat, 0.0f);
	MAKE_INPUT(nAttr)
	positionZ = nAttr.create("positionZ", "pz", MFnNumericData::kFloat, 0.0f);
	MAKE_INPUT(nAttr)

	stiffCoeff = nAttr.create("stiffCoeff", "sf", MFnNumericData::kFloat, 2.0f);
	MAKE_INPUT(nAttr)
	weakCoeff = nAttr.create("weakCoeff", "wt", MFnNumericData::kFloat, 0.6f);
	MAKE_INPUT(nAttr)
	dampCoeffV = nAttr.create("dampCoeffV", "dpv", MFnNumericData::kFloat, 0.9f);
	MAKE_INPUT(nAttr)
	dampCoeffL = nAttr.create("dampCoeffL", "dpl", MFnNumericData::kFloat, 0.9f);
	MAKE_INPUT(nAttr)
	solidAdhesionCoeff = nAttr.create("solidAdhesionCoeff", "sd", MFnNumericData::kFloat, 15.0f);
	MAKE_INPUT(nAttr)
	solidAttractionCoeff = nAttr.create("solidAttractionCoeff", "st", MFnNumericData::kFloat, 6.0f);
	MAKE_INPUT(nAttr)
	burstingCoeff = nAttr.create("burstingCoeff", "bs", MFnNumericData::kFloat, 0.02f);
	MAKE_INPUT(nAttr)
	wetness = nAttr.create("wetness", "ws", MFnNumericData::kDouble, 0.0);
	MAKE_INPUT(nAttr)
	liquidAdhesion = nAttr.create("liquidAdhesion", "lad", MFnNumericData::kFloat, 15.0f);
	MAKE_INPUT(nAttr)
	dragCoeff = nAttr.create("dragCoeff", "drg", MFnNumericData::kFloat, 0.5f);
	MAKE_INPUT(nAttr)


	bubbleNumber = nAttr.create("bubbleNumber", "bn", MFnNumericData::kInt, 500);
	MAKE_INPUT(nAttr)
	bubbleSize = nAttr.create("bubbleSize", "bus", MFnNumericData::kFloat, 1.0f);
	MAKE_INPUT(nAttr)
		
	emitterPositionX = nAttr.create("emitterPositionX", "epX", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	emitterPositionY = nAttr.create("emitterPositionY", "epY", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	emitterPositionZ = nAttr.create("emitterPositionZ", "epZ", MFnNumericData::kFloat);
	MAKE_INPUT(nAttr)
	emitterPosition = nAttr.create("emitterPosition", "ep", emitterPositionX, emitterPositionY, emitterPositionZ);
	MAKE_INPUT(nAttr)
	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder(true);

	////bubble attribute
	status = addAttribute(BubbleSimMaya::positionX);
	McheckErr(status, "ERROR adding positionX attribute\n");
	status = addAttribute(BubbleSimMaya::positionY);
	McheckErr(status, "ERROR adding positionY attribute\n");
	status = addAttribute(BubbleSimMaya::positionZ);
	McheckErr(status, "ERROR adding positionZ attribute\n");
	status = addAttribute(BubbleSimMaya::stiffCoeff);
	McheckErr(status, "ERROR adding stiffCoeff attribute\n");
	status = addAttribute(BubbleSimMaya::weakCoeff);
	McheckErr(status, "ERROR adding weakCoeff attribute\n");
	status = addAttribute(BubbleSimMaya::dampCoeffL);
	McheckErr(status, "ERROR adding dampCoeffL attribute\n");
	status = addAttribute(BubbleSimMaya::dampCoeffV);
	McheckErr(status, "ERROR adding dampCoeffV attribute\n");
	status = addAttribute(BubbleSimMaya::solidAdhesionCoeff);
	McheckErr(status, "ERROR adding solidAdhesionCoeff attribute\n");
	status = addAttribute(BubbleSimMaya::solidAttractionCoeff);
	McheckErr(status, "ERROR adding solidAttractionCoeff attribute\n");
	status = addAttribute(BubbleSimMaya::burstingCoeff);
	McheckErr(status, "ERROR adding burstingCoeff attribute\n");
	status = addAttribute(BubbleSimMaya::wetness);
	McheckErr(status, "ERROR adding wetness attribute\n");
	status = addAttribute(BubbleSimMaya::liquidAdhesion);
	McheckErr(status, "ERROR adding liquidAdhesion attribute\n");
	status = addAttribute(BubbleSimMaya::dragCoeff);
	McheckErr(status, "ERROR adding dragCoeff attribute\n");

	status = addAttribute(BubbleSimMaya::bubbleNumber);
	McheckErr(status, "ERROR adding bubbleNumber attribute\n");
	status = addAttribute(BubbleSimMaya::bubbleSize);
	McheckErr(status, "ERROR adding bubbleSize attribute\n");


    status = addAttribute(BubbleSimMaya::inNumPoints);
	McheckErr(status, "ERROR adding inNumPoints attribute\n");
/*
	status = addAttribute(BubbleSimMaya::inMinimum);
	McheckErr(status, "ERROR adding inMinimum attribute\n");
	status = addAttribute(BubbleSimMaya::inMaximum);   
	McheckErr(status, "ERROR adding inMaximum attribute\n");*/
	status = addAttribute(BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR adding outPoints attribute\n");


	status = addAttribute(BubbleSimMaya:: emitterPosition);
	McheckErr(status, "ERROR adding eTrans attribute\n");
//////////////////////////////////////////////////////////////////////////////////////////////
	status = addAttribute(BubbleSimMaya:: eTrans);
	McheckErr(status, "ERROR adding eTrans attribute\n");
	status = addAttribute(BubbleSimMaya:: eScale);   
	McheckErr(status, "ERROR adding eScale attribute\n");
	status = addAttribute(BubbleSimMaya:: fResolution);
	McheckErr(status, "ERROR adding fResolution attribute\n");
	status = addAttribute(BubbleSimMaya:: time);
	McheckErr(status, "ERROR adding time attribute\n");

	status = addAttribute(BubbleSimMaya::boundingBoxMin);
	McheckErr(status, "ERROR adding boundingBoxMin attribute\n");
	status = addAttribute(BubbleSimMaya:: boundingBoxMax);
	McheckErr(status, "ERROR adding boundingBoxMax attribute\n");
	status = addAttribute(BubbleSimMaya:: boundingBoxSize);
	McheckErr(status, "ERROR adding boundingBoxSize attribute\n");
	status = addAttribute(BubbleSimMaya:: baseResolution);
	McheckErr(status, "ERROR adding baseResolution attribute\n");
/////////////////////////////////////////////////////////////////////////////////////////////

	status = attributeAffects(BubbleSimMaya::inNumPoints, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
/*
	status = attributeAffects(BubbleSimMaya::inMinimum, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::inMaximum, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");*/
	status = attributeAffects(MPxFieldNode::mInputData, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::time, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");

	////bubble attribute
	status = attributeAffects(BubbleSimMaya::positionX, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::positionY, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::positionZ, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::stiffCoeff, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::weakCoeff, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::dampCoeffL, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::dampCoeffV, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::solidAdhesionCoeff, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::solidAttractionCoeff, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::burstingCoeff, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::wetness, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::liquidAdhesion, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::dragCoeff, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");


	status = attributeAffects(BubbleSimMaya::bubbleNumber, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::bubbleSize, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");


	status = attributeAffects(BubbleSimMaya::emitterPosition, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");

	status = attributeAffects(BubbleSimMaya::boundingBoxMin, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::boundingBoxMax, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::boundingBoxSize, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");
	status = attributeAffects(BubbleSimMaya::baseResolution, BubbleSimMaya::outPoints);
	McheckErr(status, "ERROR in attributeAffects\n");


    MGlobal::displayInfo( "BubbleSim Node initialization!\n" );
	return( MS::kSuccess );
}

int countT = 0;
int oldIteration = -1;

MStatus BubbleSimMaya::compute(const MPlug& plug, MDataBlock& block)
	//
	//	Descriptions:
	//		compute updated bubble positions and velocities
	//
{
	MStatus status;

	//cout<<count<<endl;
	//count++;

	if( !(plug == outPoints) )
		return( MS::kUnknownParameter );

	MDataHandle timeData = block.inputValue( time, &status ); 
	McheckErr(status, "Error getting time data handle\n");
	int iteration = timeData.asInt();


	MDataHandle  pxData = block.inputValue(positionX, &status);  
	McheckErr(status, "Error getting px data handle\n");
	float px = pxData.asFloat();
	MDataHandle  pyData = block.inputValue(positionY, &status);             
	McheckErr(status, "Error getting py data handle\n");
	float py = pyData.asFloat();
	MDataHandle  pzData = block.inputValue(positionZ, &status);             
	McheckErr(status, "Error getting pz data handle\n");
	float pz = pzData.asFloat();
	MDataHandle  stiffData = block.inputValue(stiffCoeff, &status);             
	float stiff = stiffData.asFloat();
	MDataHandle  weakData = block.inputValue(weakCoeff, &status);             
	float weak = weakData.asFloat();
	MDataHandle  damplData = block.inputValue(dampCoeffL, &status);             
	float dampl = damplData.asFloat();
	MDataHandle  dampvData = block.inputValue(dampCoeffV, &status);             
	float dampv = dampvData.asFloat();
	MDataHandle  solidadData = block.inputValue(solidAdhesionCoeff, &status);             
	float solidad = solidadData.asFloat();
	MDataHandle  solidatData = block.inputValue(solidAttractionCoeff, &status);             
	float solidat = solidatData.asFloat();
	MDataHandle  burstData = block.inputValue(burstingCoeff, &status);             
	float burst = burstData.asFloat();
	MDataHandle  wetnessData = block.inputValue(wetness, &status);             
	double wet = wetnessData.asDouble();
	MDataHandle  LAdhesionData = block.inputValue(liquidAdhesion, &status);             
	float lad = LAdhesionData.asFloat();
	MDataHandle  DragData = block.inputValue(dragCoeff, &status);             
	float drg = DragData.asFloat();

	MDataHandle  numberData = block.inputValue(bubbleNumber, &status);             
	int bNumber = numberData.asInt();
	MDataHandle  bsizeData = block.inputValue(bubbleSize, &status);             
	float bSize = bsizeData.asFloat();


	// obtain container information
	MDataHandle boundingBoxMinHandle = block.inputValue(boundingBoxMin, &status);
	float containerMinX = boundingBoxMinHandle.child(boundingBoxMinX).asFloat();
	float containerMinY = boundingBoxMinHandle.child(boundingBoxMinY).asFloat();
	float containerMinZ = boundingBoxMinHandle.child(boundingBoxMinZ).asFloat();

	MDataHandle boundingBoxMaxHandle = block.inputValue(boundingBoxMax, &status);
	float containerMaxX = boundingBoxMaxHandle.child(boundingBoxMaxX).asFloat();
	float containerMaxY = boundingBoxMaxHandle.child(boundingBoxMaxY).asFloat();
	float containerMaxZ = boundingBoxMaxHandle.child(boundingBoxMaxZ).asFloat();



	MDataHandle boundingBoxSizeHandle = block.inputValue(boundingBoxSize, &status);
	float containerSizeX = boundingBoxSizeHandle.child(boundingBoxSizeX).asFloat();
	float containerSizeY = boundingBoxSizeHandle.child(boundingBoxSizeY).asFloat();
	float containerSizeZ = boundingBoxSizeHandle.child(boundingBoxSizeZ).asFloat();

	float maxSize = max(containerSizeX, max(containerSizeY, containerSizeZ));

	MDataHandle baseResoHandle = block.inputValue(baseResolution, &status);
	int baseResoGrid = baseResoHandle.asInt();
	float gridSize = maxSize / baseResoGrid;
	//cout<<"containerMin "<<containerMinX<<" "<<containerMinY<<" "<<containerMinZ<<endl;
	//cout<<"containerMax "<<containerMaxX<<" "<<containerMaxY<<" "<<containerMaxZ<<endl;
	//cout<<"containerSize "<<containerSizeX<<" "<<containerSizeY<<" "<<containerSizeZ<<endl;
	//cout<<"baseReso "<<baseResoGrid<<endl;

	MArrayDataHandle emitterArrayHandle = block.inputArrayValue(emitterPosition, &status);
	McheckErr(status, "arrayHandle construction failed\n");
	unsigned count = emitterArrayHandle.elementCount();
	float Posx, Posy, Posz;
	std::vector<glm::vec3> emitter;
	for( int i = 0; i < count; i++)
	{
		emitterArrayHandle.jumpToElement(i);
		MDataHandle eHandle = emitterArrayHandle.inputValue(&status).child(emitterPositionX);
		Posx = eHandle.asFloat();
		eHandle = emitterArrayHandle.inputValue(&status).child(emitterPositionY);
		Posy = eHandle.asFloat();
		eHandle = emitterArrayHandle.inputValue(&status).child(emitterPositionZ);
		Posz = eHandle.asFloat();
		emitter.push_back(glm::vec3(Posx, Posy, Posz));
		//cout<<"emitter Pos"<<Posx<<" "<<Posy<<" "<<Posz<<endl;
		
	}	
	// Get fluid velocity field
	MArrayDataHandle hInputArray = block.outputArrayValue( mInputData, &status );
	McheckErr(status,"ERROR in hInputArray = block.outputArrayValue().\n");

	MDataHandle hCompond = hInputArray.inputValue( &status );
	McheckErr(status, "ERROR in hCompond = hInputArray.inputValue\n");

	MDataHandle hVelocity = hCompond.child( mInputVelocities );
	MObject dVelocity = hVelocity.data();
	MFnVectorArrayData fnVelocity( dVelocity );
	MVectorArray velocities = fnVelocity.array( &status );
	McheckErr(status, "ERROR in fnVelocity.array(), not find velocities.\n");

	//for(int i = 0; i< velocities.length(); i++){
	//	cout<<"x "<<velocities[i].x <<"y "<< velocities[i].y<<"z "<<velocities[i].z<<endl;
	//}
	//cout<<"I am out"<<endl;
	
	// Updata bubbles	


	//MGlobal::executeCommand("select -r polySurface2; ");
	//MGlobal::executeCommand("file -f -options \"v=0;\"  -typ \"mayaBinary\" -o \"E:/Documents/maya/projects/Project3_Character/scenes/660.mb\";addRecentFile(\"E:/Documents/maya/projects/Project3_Character/scenes/660.mb\", \"mayaBinary\"); ");

	
	myBubble.setParameter(stiff, weak, dampv/100000.0, dampl/100000.0, solidad, solidat, burst, wet, lad, drg);
	//cout<<stiff<<" "<<weak<<" "<<dampl<<" "<<dampv<<" "<<solidad<<" "<<solidat<<" "<<burst<<endl;
	//cout<<"iteration "<<iteration<<endl;
	//for(int i = 0; i < (iteration/250.0f); i++){
	//cout<<"bubble Number "<<bNumber<<endl;
	//countT++;
	if(oldIteration != iteration){
		
		MObject pluginObj = MFnPlugin::findPlugin("BubbleSimMaya");
		MFnPlugin plugin(pluginObj);
		MString pathMEL = plugin.loadPath() + "/waterMesh.obj";
		//MGlobal::executeCommand("string $objPath = " + pathMEL + ";");
		
		MGlobal::executeCommand("select surfaceMesh; ");	
		//MGlobal::executeCommand("select -r polySurface2; ");
		MGlobal::executeCommand("file -force -options \"groups=0;ptgroups=0;materials=0;smoothing=0;normals=0\" -typ \"OBJexport\" -pr -es\"" + pathMEL + "\";");
		MGlobal::executeCommand("select BubbleEmitter0; ");
		cout<<"iteration is "<<iteration<<endl;
		if(iteration %20==0 && myBubble.size() < bNumber){
			for(int i = 0; i < emitter.size(); i++)
				myBubble.insertBubble( glm::vec3(emitter[i][0] +(rand() % 2 /10.0f), emitter[i][1], emitter[i][2] + (rand() % 2/10.f)) / SimScale, glm::vec3(0,0,0), (bSize/10.f + (rand() % (int)(bSize*2>1?2:bSize*2))/10.f) / SimScale);		
		}		


		//cal the signed distance to surface
		MString plugM = plugin.loadPath();
		std::string loadPath = plugM.asChar();
		cout<<loadPath<<endl;
		std::string filename(loadPath + "/waterMesh.obj");//the obj file
		cout<<"complete0"<<endl;
		float dx;
		dx = gridSize;//grid size

		int padding;
		padding = 3;

		//if(padding < 1) padding = 1;

		Vec3f min_box(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()), 
			max_box(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());
		cout<<"complete1"<<endl;
		std::ifstream infile(filename);
		if(!infile) {
			std::cerr << "Failed to open. Terminating.\n";
			//exit(-1);
		}
		cout<<"complete2"<<endl;
		int ignored_lines = 0;
		std::string line;
		std::vector<Vec3f> vertList;
		std::vector<Vec3ui> faceList;
		while(!infile.eof()) {
			std::getline(infile, line);
			if(line.substr(0,1) == std::string("v")) {
				std::stringstream data(line);
				char c;
				Vec3f point;
				data >> c >> point[0] >> point[1] >> point[2];
				vertList.push_back(point);
				update_minmax(point, min_box, max_box);
			}
			else if(line.substr(0,1) == std::string("f")) {
				std::stringstream data(line);
				char c;
				int v0,v1,v2;
				data >> c >> v0 >> v1 >> v2;
				faceList.push_back(Vec3ui(v0-1,v1-1,v2-1));
			}
			else {
				++ignored_lines; 
			}
		}
		cout<<"complete3"<<endl;
		infile.close();
		glm::vec3 minSolid, maxSolid;
		minSolid[0] = min_box[0]; minSolid[1] = min_box[1]; minSolid[2] = min_box[2];
		maxSolid[0] = max_box[0]; maxSolid[1] = max_box[1]; maxSolid[2] = max_box[2];

		if(countT == 0){			
			myBubble.addSoild(glm::vec3(-1,0,0), glm::vec3(maxSolid[0], maxSolid[1], maxSolid[2]));
			myBubble.addSoild(glm::vec3(1,0,0), glm::vec3(minSolid[0], minSolid[1], minSolid[2]));
			myBubble.addSoild(glm::vec3(0,0,-1), glm::vec3(maxSolid[0], maxSolid[1], maxSolid[2]));
			myBubble.addSoild(glm::vec3(0,0,1), glm::vec3(minSolid[0], minSolid[1], minSolid[2]));
		}

		countT++;

		cout<<"complete4"<<endl;
		Vec3f unit(1,1,1);
		min_box -= padding*dx*unit;
		max_box += padding*dx*unit;
		Vec3ui sizes = Vec3ui((max_box - min_box)/dx);

		//cout<<"make level"<<endl;
		Array3f phi_grid;
		make_level_set3(faceList, vertList, min_box, dx, sizes[0], sizes[1], sizes[2], phi_grid);//This is how they make the level set
		cout<<"complete5"<<endl;
		int dimensionX = phi_grid.ni;
		int dimensionY = phi_grid.nj;
		int dimensionZ = phi_grid.nk;

		float minboxX = min_box[0];
		float minboxY = min_box[1];
		float minboxZ = min_box[2];
		cout<<"minBox "<<minboxX<<" "<<minboxY<<" "<<minboxZ<<endl;

		for(int i = 0; i < myBubble.size(); i++){
			
			glm::vec3 bP = myBubble.getPosition(i);
			bP = bP * SimScale;
			cout<<"bp "<<bP[0]<<" "<<bP[1]<<" "<<bP[2]<<"                ";			

			float dtx = bP[0] - minboxX;float dty = bP[1] - minboxY;float dtz = bP[2] - minboxZ;
			
			int ii = dtx/dx; int jj = dty/dx; int kk = dtz/dx;

			if(ii <= 0 || jj <= 0 || kk <= 0 || ii >= sizes[0]-1 || jj >= sizes[1]-1 || kk >= sizes[2]-1){
				myBubble.setBubbleType(i,10.f);//set it as air bubble
				cout<<"out of bound"<<endl;				
				continue;
			}

			int indexSignedD = ii + jj * dimensionX + kk * dimensionY * dimensionX;
			int indexXneighbour = ii-1 + jj * dimensionX + kk * dimensionY * dimensionX;
			int indexYneighbour = ii + (jj-1) * dimensionX + kk * dimensionY * dimensionX;
			int indexZneighbour = ii+ jj * dimensionX + (kk-1) * dimensionY * dimensionX;
			glm::vec3 surfaceN;
			if(ii <= 0 || jj <= 0 || kk <= 0){
				cout<<"Wo Diu"<<endl;
				surfaceN = glm::vec3(0,1,0);
			}
			else{
				surfaceN[0] = phi_grid.a[indexSignedD] - phi_grid.a[indexXneighbour];
				surfaceN[1] = phi_grid.a[indexSignedD] - phi_grid.a[indexYneighbour];
				surfaceN[2] = phi_grid.a[indexSignedD] - phi_grid.a[indexZneighbour];
				surfaceN /= dx;
			}
			//cout<<"is here?"<<endl;
			myBubble.setBubbleType(i,phi_grid.a[indexSignedD]/SimScale);//set bubble type	

			myBubble.setSurfaceNormal(i,surfaceN);//set surface normal
			//cal velocity
			dtx = bP[0] - containerMinX; dty = bP[1] - containerMinY; dtz = bP[2] - containerMinZ;
			ii = dtx/dx; jj = dty/dx; kk = dtz/dx;
			indexSignedD = ii + jj * containerSizeX + kk * containerSizeY * containerSizeX;
			//cout<<"Or is here?"<<endl;
			glm::vec3 liquidVel = glm::vec3(0,0,0);
			if(indexSignedD < velocities.length()){
				liquidVel[0] = velocities[indexSignedD].x; liquidVel[1] = velocities[indexSignedD].y; liquidVel[2] = velocities[indexSignedD].z;
			}
			else
				cout<<"velocity is out!!"<<endl;
			myBubble.setLiquidVel(i,liquidVel);	
		}
		 
		myBubble.update(1.f/192.0f);
		


		for(int i = 0; i < myBubble.size(); i++){

			glm::vec3 bP = myBubble.getPosition(i);
			bP = bP * SimScale;
			cout<<"bp "<<bP[0]<<" "<<bP[1]<<" "<<bP[2]<<"                ";			

			float dtx = bP[0] - minboxX;float dty = bP[1] - minboxY;float dtz = bP[2] - minboxZ;

			int ii = dtx/dx; int jj = dty/dx; int kk = dtz/dx;

			if(ii <= 0 || jj <= 0 || kk <= 0 || ii >= sizes[0]-1 || jj >= sizes[1]-1 || kk >= sizes[2]-1){
				myBubble.setBubbleType(i,10.f);//set it as air bubble
				cout<<"out of bound"<<endl;
				continue;
			}

			int indexSignedD = ii + jj * dimensionX + kk * dimensionY * dimensionX;
			int indexXneighbour = ii-1 + jj * dimensionX + kk * dimensionY * dimensionX;
			int indexYneighbour = ii + (jj-1) * dimensionX + kk * dimensionY * dimensionX;
			int indexZneighbour = ii+ jj * dimensionX + (kk-1) * dimensionY * dimensionX;
			glm::vec3 surfaceN;
			if(ii <= 0 || jj <= 0 || kk <= 0){
				cout<<"Wo Diu"<<endl;
				surfaceN = glm::vec3(0,1,0);
			}
			else{
				surfaceN[0] = phi_grid.a[indexSignedD] - phi_grid.a[indexXneighbour];
				surfaceN[1] = phi_grid.a[indexSignedD] - phi_grid.a[indexYneighbour];
				surfaceN[2] = phi_grid.a[indexSignedD] - phi_grid.a[indexZneighbour];
				surfaceN /= dx;
			}
			//cout<<"is here?"<<endl;
			myBubble.setBubbleType(i,phi_grid.a[indexSignedD]/SimScale);//set bubble type	

			myBubble.setSurfaceNormal(i,surfaceN);//set surface normal
			//cal velocity
			dtx = bP[0] - containerMinX; dty = bP[1] - containerMinY; dtz = bP[2] - containerMinZ;
			ii = dtx/dx; jj = dty/dx; kk = dtz/dx;
			indexSignedD = ii + jj * containerSizeX + kk * containerSizeY * containerSizeX;
			//cout<<"Or is here?"<<endl;
			glm::vec3 liquidVel = glm::vec3(0,0,0);
			if(indexSignedD < velocities.length()){
				liquidVel[0] = velocities[indexSignedD].x; liquidVel[1] = velocities[indexSignedD].y; liquidVel[2] = velocities[indexSignedD].z;
			}
			else
				cout<<"velocity is out!!"<<endl;
			myBubble.setLiquidVel(i,liquidVel);	
		}

		myBubble.update(1.f/192.0f);




		oldIteration = iteration;
	}


	//}
//cout<<"x "<<px<<" y "<<py<<" z "<<pz<<endl;0

	/*5.0f + rand() % 12*/

	/*MDataHandle  numData= block.inputValue(inNumPoints, &status);             
	int numPts = numData.asInt();
	MDataHandle  MiXData = block.inputValue(inMinimumX, &status);
	float MinimumX = MiXData.asFloat();
	MDataHandle  MiYData = block.inputValue(inMinimumY, &status);
	float MinimumY = MiYData.asFloat();
	MDataHandle  MiZData = block.inputValue(inMinimumZ, &status);
	float MinimumZ = MiZData.asFloat();

	MDataHandle  MaXData = block.inputValue(inMaximumX, &status);
	float MaximumX = MaXData.asFloat();
	MDataHandle  MaYData = block.inputValue(inMaximumY, &status);
	float MaximumY = MaYData.asFloat();
	MDataHandle  MaZData = block.inputValue(inMaximumZ, &status);
	float MaximumZ = MaZData.asFloat();*/

	MDataHandle  pointsData = block.outputValue(outPoints, &status);
	MFnArrayAttrsData pointsAAD;           
	MObject pointsObject = pointsAAD.create();   

	MVectorArray positionArray = pointsAAD.vectorArray("position");
	MVectorArray scaleArray = pointsAAD.vectorArray("scale");
	MDoubleArray idArray = pointsAAD.doubleArray("id");

//	srand (time(NULL));
	float x, y, z;
	float scaleX, scaleY, scaleZ;
	for(int i=0; i <myBubble.bubbleArr.size(); i++)
	{
		x = myBubble.bubbleArr[i].center[0] * SimScale;
		y = myBubble.bubbleArr[i].center[1] * SimScale;
		z = myBubble.bubbleArr[i].center[2] * SimScale;
		scaleX = myBubble.bubbleArr[i].radius * SimScale;
		scaleY = myBubble.bubbleArr[i].radius * SimScale;
		scaleZ = myBubble.bubbleArr[i].radius * SimScale;
		MVector point(x, y, z);
		MVector scale(scaleX, scaleY, scaleZ);
		positionArray.append(point);
		scaleArray.append(scale);
		idArray.append(i);
	}

	pointsData.set(pointsObject);
	block.setClean( plug );
	


	return( MS::kSuccess );
}

