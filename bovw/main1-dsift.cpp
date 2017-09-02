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

extern "C" {
	#include <vl/generic.h>
}


#include <vl/dsift.h>
#include <vl/kmeans.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <iostream>

namespace main1
{
	using namespace cv;

	int main () {
		printf("Este es el main %d \n", 1);

		VL_PRINT("Hello world!\n");

		Mat image = imread("C:/2.jpg");

		// create filter
		VlDsiftFilter *vlf = vl_dsift_new_basic(image.rows, image.cols, 1, 3);

		// transform image in cv::Mat to float vector
		std::vector<float> imgvec;
		for (int i = 0; i < image.rows; ++i){
			for (int j = 0; j < image.cols; ++j){
				imgvec.push_back(image.at<unsigned char>(i, j) / 255.0f);
			}
		}
		
		// call processing function of vl
		vl_dsift_process(vlf, &imgvec[0]);

		// echo number of keypoints found
		int numKeyPoints = vl_dsift_get_keypoint_num(vlf);
		std::cout << "Points: " << numKeyPoints << std::endl;

		// Extract keypoints
		VlDsiftKeypoint const *vlkeypoints;
		vlkeypoints = vl_dsift_get_keypoints(vlf);

		for (int i = 0; i < 5; i++) {
			std::cout << "X4: " << vlkeypoints[i].x << std::endl;
			std::cout << "Y4: " << vlkeypoints[i].y << std::endl;
			std::cout << std::endl;
		}
		
		return 0;
	}

}