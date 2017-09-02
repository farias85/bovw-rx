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

//#include <opencv\cxcore.h>
//#include <opencv\highgui.h>

#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\features2d\features2d.hpp>

#include <vector>
#include <iostream>

using namespace cv;
using namespace std;

namespace main5
{
	int main() {

		printf("Este es el main %d \n", 5);

		Mat image = imread("C:/2.jpg", CV_LOAD_IMAGE_COLOR);

		//IplImage* img = cvLoadImage("C:/2.jpg");
		//Mat image = Mat(img);

		if (!image.data)
		{
			cout << "Could not open or find the image" << std::endl;
			waitKey(5000);
			return -1;
		}

		int rows = image.rows;
		int cols = image.cols;

		namedWindow("Display window", CV_WINDOW_AUTOSIZE);
		imshow("Display window", image);

		cvtColor(image, image, CV_BGR2HSV);

		vector<Mat> layers;
		split(image, layers);

		//cvSmooth

		imshow("layer 0", layers[0]);
		waitKey();
		//equalizeHist(layers[0], layers[0]);
		imshow("layer 0", layers[0]);
		waitKey();

		imshow("layer 0", layers[1]);
		waitKey();
		//equalizeHist(layers[1], layers[1]);
		imshow("layer 0", layers[1]);
		waitKey();

		imshow("layer 0", layers[2]);
		waitKey();
		//equalizeHist(layers[2], layers[2]);
		imshow("layer 0", layers[2]);
		waitKey();

		merge(layers, image);
		cvtColor(image, image, CV_HSV2BGR);

		imshow("Final", image);
		waitKey();

		//system("pause");
	}

}