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

#include <vl/homkermap.h>
#include <vl/mathop.h>
#include <vl/stringop.h>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

namespace main10
{

	void main()
	{
		char* dir = "C:/histogram-neg.txt";
		fstream myfile(dir, ios_base::in);

		float value;
		int valueCount = 0;
		vl_size const histSize = 1000;
		float* histogram = new float[histSize];

		while (myfile >> value)
		{
			histogram[valueCount] = value;
			valueCount++;
		}
		myfile.close();

		vl_size const kernelSize = 5;
		double gamma = 1.0;
		vl_size order = 2;
		double period = -1;
		vl_size psiStride = 1;

		VlHomogeneousKernelMap * hom = vl_homogeneouskernelmap_new(VlHomogeneousKernelChi2,
			gamma, order, period, VlHomogeneousKernelMapWindowRectangular);

		float* psi = new float[kernelSize];
		vl_size const resultSize = kernelSize * histSize;
		float* resultHistogram = new float[resultSize];

		for (size_t i = 0; i < histSize; i++)
		{
			//pass data througt the kernel
			vl_homogeneouskernelmap_evaluate_f(hom, psi, psiStride, histogram[i]);
			for (size_t j = 0; j < kernelSize; j++)
				resultHistogram[(i * kernelSize) + j] = psi[j];
		}

		//for (int i = 0; i < resultSize; i++)
		//	cout << "Valor " << i << ": " << result[i] << endl;

		vl_homogeneouskernelmap_delete(hom);

		char* dirW = "C:/weight.svm";
		fstream svmfile(dirW, ios_base::in);

		valueCount = 0;
		vl_size const weightSize = 5000;
		float* svm_weight = new float[weightSize];

		while (svmfile >> value)
		{
			svm_weight[valueCount] = value;
			valueCount++;
		}
		svmfile.close();

		float svm_bias = -1.081984805533894f;
		float svm_threshold = -1.0167452f;

		//compute linear proyection
		double score = 0;
		for (size_t i = 0; i < weightSize; i++)
			score += svm_weight[i] * resultHistogram[i];
		score += svm_bias;

		int classT = 2 * (score > svm_threshold) - 1;

		cout << "Score: " << score << endl;
		cout << "Tipo: " << classT << endl;
	}

}