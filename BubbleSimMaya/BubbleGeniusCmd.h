#ifndef CreateBubbleGeniusCmd_H_
#define CreateBubbleGeniusCmd_H_

#include <maya/MPxCommand.h>
#include <string>
#include "extern.h"

class BubbleGeniusCmd : public MPxCommand
{
public:
	BubbleGeniusCmd();
	virtual ~BubbleGeniusCmd();
	static void* creator() { return new BubbleGeniusCmd(); }
	MStatus doIt( const MArgList& args );
	static MSyntax newSyntax();

};

#endif