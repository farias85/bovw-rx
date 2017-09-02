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

#include "main12-completo.cpp"
#include "BoVW.h";

using namespace std;
using namespace cv;
using namespace img;

#include <vl/float.th>
#include <vl/host.h>

#include <string.h>

void all()
{
	size_t length = 145;
	string dir = "C:/_old/por_clases/pistola+bg/format/";

	//size_t length = 146;
	//string dir = "C:/_old/por_clases/pistola_nobg/format/";
	
	//size_t length = 171;
	//string dir = "C:/_old/por_clases/revolver+bg/format/";

	//size_t length = 171;
	//string dir = "C:/_old/por_clases/revolver_nobg/format/";

	//size_t length = 8706;
	//string dir = "C:/_old/por_clases/no_pistola2/format/";

	bovw::BoVW b;
	int target = 0;

	for (int i = 1; i <= length; i++)
	{
		string full = dir + "pistola" + to_string(i) + ".png";
		//string full = dir + "revolver" + to_string(i) + ".png";
		//string full = dir + "no_pistola" + to_string(i) + ".png";

		Mat image = imread(full.data(), CV_LOAD_IMAGE_COLOR);
		int classT = b.computeAlgorithm(image);

		if (classT == 1)
			target++;

		cout << i << " => " << classT << endl;
		//imshow(full, image);
	}

	cout << "Porciento de acierto " << (double)target / (double)length * 100 << "%" << endl;

	waitKey();
}

void one()
{
	Mat image = imread("C:/_old/bovw_matrix/pistola_1.png", CV_LOAD_IMAGE_COLOR);
	//Mat image = imread("C:/_old/bovw_matrix/pistola_43.png", CV_LOAD_IMAGE_COLOR);
	//Mat image = imread("C:/_old/bovw_matrix/pistola_55.png", CV_LOAD_IMAGE_COLOR);

	//Mat image = imread("C:/_old/bovw_matrix/no_pistola_3.png", CV_LOAD_IMAGE_COLOR);
	//Mat image = imread("C:/_old/bovw_matrix/no_pistola_4.png", CV_LOAD_IMAGE_COLOR);
	//Mat image = imread("C:/_old/bovw_matrix/no_pistola_10.png", CV_LOAD_IMAGE_COLOR);

	double t1 = omp_get_wtime();
	bovw::BoVW b;
	double t2 = omp_get_wtime();
	double t = t2 - t1;
	std::cout << "Time - bovw::BoVW b;!!! - " << t << std::endl << std::endl;

	t1 = omp_get_wtime();
	int classT = b.computeAlgorithm(image);
	t2 = omp_get_wtime();
	t = t2 - t1;
	std::cout << "CLASS -> " << classT << std::endl;
	std::cout << "Time - b.computeAlgorithm(image);!!! - " << t << std::endl << std::endl;
}

void main(int argc, const char * argv[])
{
	//main0::main(); 
	//main1::main();
	//main2::main();
	//main3::main();
	//main4::main();
	//main5::main();
	//main6::main(argc, argv);
	//main7::main(); 
	//main8::main();
	//main9::main(argc, argv);
	//main10::main();
	//main11::main();
	//main12::main(argc, argv);   	

	one();

	//double t1 = omp_get_wtime();
	//all();
	//double t2 = omp_get_wtime();
	//double t = t2 - t1;
	//std::cout << "Time - all();!!! - " << t << std::endl << std::endl;
	

	std::system("pause");
}

