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

#include <stdio.h>
#include <vl/svm.h>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace main2
{
	int main()
	{
		printf("Este es el main %d \n", 2);

		cv::Mat image = cv::imread("C:/1.png");
		cv::namedWindow("My Image");
		cv::imshow("My Image", image);
		cv::waitKey(5000);

		vl_size const numData = 4;
		vl_size const dimension = 2;

		double x[dimension * numData] = {
			0.0, -0.5,
			0.6, -0.3,
			0.0, 0.5,
			0.6, 0.0 };

		double y[numData] = { 1, 1, -1, 1 };
		double lambda = 0.01;

		//double *const model;

		double bias;
		VlSvm * svm = vl_svm_new(VlSvmSolverSgd,
			x, dimension, numData,
			y,
			lambda);
		vl_svm_train(svm);
		//vl_svmdataset_
		//vl_homogeneouskernelmap

		const double * model = vl_svm_get_model(svm);
		bias = vl_svm_get_bias(svm);

		printf("model w = [ %f , %f ] , bias b = %f \n",
			model[0],
			model[1],
			bias);

		VlSvmStatistics const* st = vl_svm_get_statistics(svm);

		VL_PRINT("Hello:");
		std::cout << st->elapsedTime << std::endl;

		vl_svm_delete(svm);

		system("pause");
		return 0;
	}
}