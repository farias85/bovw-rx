/**
* Licensed under the Apache License, Version 2.0 (the "License"); you may not
* use this file except in compliance with the License. You may obtain a copy of
* the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
* License for the specific language governing permissions and limitations under
* the License.
*
* Created by Felipe Rodriguez Arias <ucifarias@gmail.com>.
*/

#include "Util.h"

namespace bovw
{
	Util::Util()
	{
	}

	Util::~Util()
	{
	}

	float* Util::loadFloatArray(char* dir, int dataSize)
	{
		float value = 0;
		int valueCount = 0;
		float* data = new float[dataSize];

		std::fstream stream(dir, std::ios_base::in);
		while (stream >> value)
		{
			data[valueCount] = value;
			valueCount++;
		}
		stream.close();

		return data;
	}

	int* Util::loadIntegerArray(char* dir, int imgvecSize)
	{
		float value = 0;
		int valueCount = 0;
		int* data = new int[imgvecSize];

		std::fstream stream(dir, std::ios_base::in);
		while (stream >> value)
		{
			data[valueCount] = value;
			valueCount++;
		}
		stream.close();

		return data;
	}

	char* Util::pathTo(const char *src)
	{
		const char* path = "C:/_old/bovw/";

		int length = std::strlen(path);
		char* str = new char[length + strlen(src)];
		
		std::strcpy(str, path);
		std::strcpy(str + length, src);
		
		return str;
	}
}
