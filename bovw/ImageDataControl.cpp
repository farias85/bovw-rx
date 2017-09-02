#include "ImageDataControl.h"

namespace img
{
	ImageDataControl::ImageDataControl(long startByte, int columnLength, int zeroPadding)
	{
		this->startByte = startByte;
		this->columnLength = columnLength;
		this->zeroPadding = zeroPadding;
	}

	ImageDataControl::~ImageDataControl()
	{
	
	}
}
