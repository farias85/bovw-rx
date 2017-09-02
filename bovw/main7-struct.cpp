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

#include <algorithm>
#include <iostream>
#include <vector>

namespace main7
{

	class A{
	public:
		int i;
	};

	class A2 : A{
	};

	struct A3 :A{
	};


	struct abc{
		int i;
	};

	struct abc2 :abc{
	};

	class abc3 :abc{
	};

	
	int main ()
	{
		printf("Este es el main %d \n", 7);

		abc2 objabc;
		objabc.i = 10;
		printf("Struct: %d \n", objabc.i);

		A3 ob;
		ob.i = 10;
		printf("Class: %f \n", (double)ob.i);

		//A2 obja; //privately inherited
		//obja.i = 10;

		//abc3 obss;
		//obss.i = 10;

		system("pause");
		return 1;
	}

}