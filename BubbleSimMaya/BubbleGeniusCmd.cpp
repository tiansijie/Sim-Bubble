#include "BubbleGeniusCmd.h"

#include <maya/MGlobal.h>
#include <list>
#include <maya/MArgList.h>
#include <maya/MSyntax.h>
#include <maya/MString.h>
#include <maya/MArgDatabase.h>


BubbleGeniusCmd::BubbleGeniusCmd() : MPxCommand()
{
}

BubbleGeniusCmd::~BubbleGeniusCmd() 
{
}

const char *postionXFlag = "-px", *postionXLongFlag = "-postionx";
const char *postionYFlag = "-py", *postionYLongFlag = "-postiony";
const char *postionZFlag = "-pz", *postionZLongFlag = "-postionz";

const char *wetnessFlag = "-w", *wetnessLongFlag = "-wetness";
const char *stiffnessFlag = "-s", *stiffnessLongFlag = "-stiffness";
const char *dragFlag = "-d", *dragLongFlag = "-drag";
const char *liquidAdFlag = "-ld", *liquidLongFlag = "-liquidAd";
const char *bubbleAtFlag = "-lt", *bubbleLongFlag = "-bubbleAt";
const char *solidAdFlag = "-sd", *solidAdLongFlag = "-solidAd";
const char *solidAtFlag = "-st", *solidAtLongFlag = "-solidAt";
const char *dampingFlag = "-dp", *dampingLongFlag = "-damping";
const char *volumeCFlag = "-vc", *volumeCLongFlag = "-volumeC";
const char *burstingSFlag = "-bs", *burstingSLongFlag = "-burstingS";


MSyntax BubbleGeniusCmd::newSyntax(){
	MSyntax syntax;
	syntax.addFlag(postionXFlag, postionXLongFlag, MSyntax::kDouble);
	syntax.addFlag(postionYFlag, postionYLongFlag, MSyntax::kDouble);
	syntax.addFlag(postionZFlag, postionZLongFlag, MSyntax::kDouble);

	syntax.addFlag (wetnessFlag, wetnessLongFlag, MSyntax::kDouble);
	syntax.addFlag (stiffnessFlag, stiffnessLongFlag, MSyntax::kDouble);
	syntax.addFlag (dragFlag, dragLongFlag, MSyntax::kDouble);
	syntax.addFlag (liquidAdFlag, liquidLongFlag, MSyntax::kDouble);
	syntax.addFlag (bubbleAtFlag, bubbleLongFlag, MSyntax::kDouble);
	syntax.addFlag (solidAdFlag, solidAdLongFlag, MSyntax::kDouble);
	syntax.addFlag (solidAtFlag, solidAtLongFlag, MSyntax::kDouble);
	syntax.addFlag (dampingFlag, dampingLongFlag, MSyntax::kDouble);
	syntax.addFlag (volumeCFlag, volumeCLongFlag, MSyntax::kDouble);
	syntax.addFlag (burstingSFlag, burstingSLongFlag, MSyntax::kDouble);
	return syntax;
}

MStatus BubbleGeniusCmd::doIt( const MArgList& args )
{

	/*double stepsize,degree = 90;
	MString grammar = "hi";*/
	double px,py,pz;
	double wetness, stiffness, drag, liquidAd, bubbleAt, solidAd, solidAt, damping, volumeC, burstingS;
	int iter;

	MArgDatabase argData(syntax(), args);

	if(argData.isFlagSet(postionXFlag))
		argData.getFlagArgument(postionXFlag,0,px);
	if(argData.isFlagSet(postionYFlag))
		argData.getFlagArgument(postionYFlag,0,py);
	if(argData.isFlagSet(postionZFlag))
		argData.getFlagArgument(postionZFlag,0,pz);


	if(argData.isFlagSet(wetnessFlag))
		argData.getFlagArgument(wetnessFlag,0,wetness);

	if(argData.isFlagSet(stiffnessFlag))
		argData.getFlagArgument(stiffnessFlag,0,stiffness);

	if(argData.isFlagSet(dragFlag))
		argData.getFlagArgument(dragFlag,0,drag);

	if(argData.isFlagSet(liquidAdFlag))
		argData.getFlagArgument(liquidAdFlag,0,liquidAd);

	if(argData.isFlagSet(bubbleAtFlag))
		argData.getFlagArgument(bubbleAtFlag,0,bubbleAt);

	if(argData.isFlagSet(solidAdFlag))
		argData.getFlagArgument(solidAdFlag,0,solidAd);

	if(argData.isFlagSet(solidAtFlag))
		argData.getFlagArgument(solidAtFlag,0,solidAt);

	if(argData.isFlagSet(dampingFlag))
		argData.getFlagArgument(dampingFlag,0,damping);

	if(argData.isFlagSet(volumeCFlag))
		argData.getFlagArgument(volumeCFlag,0,volumeC);

	if(argData.isFlagSet(burstingSFlag))
		argData.getFlagArgument(burstingSFlag,0,burstingS);


	//LSystem system;
	//std::string grammarString = grammar.asChar();
	/*system.loadProgramFromString(grammarString);
	system.setDefaultAngle(degree);
	system.setDefaultStep(stepsize);*/
	
	wetness, stiffness, drag, liquidAd, bubbleAt, solidAd, solidAt, damping, volumeC, burstingS;
	cout<<"Implement Me!"<<endl;
	cout<<"wetness "<<wetness <<" stiffness "<<stiffness<<" drag "<< drag<<" liquidAd "<<liquidAd<<" bubbleAt "<<bubbleAt<<" solidAd "<<solidAd
		<<" solidAt "<<solidAt<<" solidAt "<<solidAt<<" damping "<<damping<<" volumeC "<<volumeC<<" burstingS "<<burstingS<<endl;

	cout<<"position X "<<px<<" position Y "<<py<<" position Z "<<pz<<endl;

	EpositionX = px; EpositionY = py; EpositionZ = pz;

	Ewetness = wetness; Estiffness = stiffness; Edrag = drag; EliquidAd = liquidAd; EbubbleAt = bubbleAt; EsolidAd = solidAd; EsolidAt = solidAt; Edamping = damping;
	EvolumeC = volumeC; EburstingS = burstingS;

	return MStatus::kSuccess;
}

