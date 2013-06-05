//
// Copyright (C) CIS660
// 
// File: pluginMain.cpp
//
// Author: Maya Plug-in Wizard 2.0.1
// Revised by Isaac Peral
//

#include "BubbleSimMayaNode.h"
#include "BubbleEmitterMayaNode.h"
#include "BubbleGeniusCmd.h"

#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <cstdio>

MStatus initializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is loaded into Maya.  It 
//		registers all of the services that this plug-in provides with 
//		Maya.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "CIS660", "2012", "Any");

	status = plugin.registerNode( "BubbleSimMaya", BubbleSimMaya::id, BubbleSimMaya::creator,
								  BubbleSimMaya::initialize, MPxNode::kFieldNode );
	if (!status) {
		status.perror("Error registerNode BubbleSimMaya");
		return status;
	}

	status = plugin.registerNode( "BubbleEmitterMaya", BubbleEmitterMaya::id, BubbleEmitterMaya::creator,
		BubbleEmitterMaya::initialize );
	if (!status) {
		status.perror("Error registerNode BubbleEmitterMaya");
		return status;
	}

	status = plugin.registerCommand( "BubbleGeniusCmd", BubbleGeniusCmd::creator, BubbleGeniusCmd::newSyntax );
	if (!status) {
		status.perror("registerCommand");
		return status;
	}

	MGlobal::executeCommand("source \"" + plugin.loadPath() + "/BubbleSimUI.mel\"");
	status = plugin.registerUI("createBubbleSimMenu", "deleteBubbleSimMenu");

	return status;
}

MStatus uninitializePlugin( MObject obj)
//
//	Description:
//		this method is called when the plug-in is unloaded from Maya. It 
//		deregisters all of the services that it was providing.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{
	MStatus   status;
	MFnPlugin plugin( obj );

	status = plugin.deregisterNode( BubbleSimMaya::id );
	if (!status) {
		status.perror("Error deregisterNode BubbleSimMaya");
		return status;
	}

	status = plugin.deregisterNode( BubbleEmitterMaya::id );
	if (!status) {
		status.perror("Error deregisterNode BubbleEmitterMaya");
		return status;
	}

	status = plugin.deregisterCommand( "BubbleGeniusCmd" );
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	/*MObject pluginObj = MFnPlugin::findPlugin("BubbleSimMaya");
	MFnPlugin plugin(pluginObj);*/
	MString pathMEL = plugin.loadPath() + "/waterMesh.obj";
	
	MString command1 = "sysFile -delete \"" + pathMEL +"\";";
	
	//MGlobal::executeCommand("string $objPath = `getenv \"MAYA_LOCATION\"`+ \"/bin/waterMesh.obj\"; ");
	//MGlobal::executeCommand("string $objPath = " + plugin.loadPath() + "/waterMesh.obj\"; ");
	MGlobal::executeCommand(command1);
	
	/*if( remove( "waterMesh.obj" ) != 0 )
		std::cout<< "Error deleting temporary fluid mesh file" <<std::endl;
	else
		std::cout<<"Temporary fluid mesh file successfully deleted" <<std::endl;
*/

	return status;
}
