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

// even_lambda.cpp
// compile with: cl /EHsc /nologo /W4 /MTd
#include <algorithm>
#include <iostream>
#include <vector>

using namespace std;

namespace main8 
{
	int main()
	{
		// Create a vector object that contains 10 elements.
		vector<int> v;
		for (int i = 1; i < 10; ++i) {
			v.push_back(i);
		}

		// Count the number of even numbers in the vector by 
		// using the for_each function and a lambda.
		int evenCount = 0;
		for_each(v.begin(), v.end(), [&evenCount](int n) {
			cout << n;
			if (n % 2 == 0) {
				cout << " is even " << endl;
				++evenCount;
			}
			else {
				cout << " is odd " << endl;
			}
		});

		// Print the count of even numbers to the console.
		cout << "There are " << evenCount
			<< " even numbers in the vector." << endl;

		return 1;
	}
}