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

#include <vl/kmeans.h>
#include <iostream>
#include <vl/kdtree.h>
#include <vl/ikmeans.h>

#include <algorithm>
#include <iterator>
#include <fstream>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <vl/mathop.h>
#define FLOAT_MAX 3.40e38

using namespace std;
using namespace cv;

namespace main9
{
	float* vl_alldist(size_t dimension, float* inX, size_t numDataX, float* inY, size_t numDataY)
	{
		vl_bool autoComparison = VL_FALSE;
		VlVectorComparisonType comparisonType = VlDistanceL2;

		float* result = new float[numDataY * numDataX];

		VlFloatVectorComparisonFunction f = vl_get_vector_comparison_function_f(comparisonType);

		vl_eval_vector_comparison_on_all_pairs_f(result,
			dimension,
			inX, numDataX,
			inY, numDataY,
			f);

		return result;
	}

	float* computeHistogram(int* clusterAssignments, int numDataX, int numCenters)
	{
		float* histogram = new float[numCenters];

		for (int i = 0; i < numCenters; i++)
			histogram[i] = 0;

		for (int i = 0; i < numDataX; i++)
			histogram[clusterAssignments[i]]++;

		for (int i = 0; i < numCenters; i++)
			histogram[i] /= (numDataX * 4);

		return histogram;
	}

	void kmeans_histogram()
	{
		double energy;

		vl_size const numData = 5000;
		vl_size const dimension = 2;

		//Cargando el archivo generado en matlab
		char* dir = "C:/data.txt";
		fstream myfile(dir, ios_base::in);

		float a;
		int count = 0;
		float* data = new float[numData * dimension];
		while (myfile >> a)
		{
			//printf("Valor: %f \n", a);
			data[count] = a;
			count++;
		}
		myfile.close();

		// Use float data and the L2 distance for clustering
		VlKMeans * kmeans = vl_kmeans_new(VL_TYPE_FLOAT, VlDistanceL2);

		// Use Elkan algorithm
		vl_kmeans_set_algorithm(kmeans, VlKMeansLloyd);

		// Initialize the cluster centers by randomly sampling the data
		//vl_kmeans_init_centers_with_rand_data (kmeans, data, dimension, numData, numCenters) ;
		int numCenters = 3;
		//vl_kmeans_init_centers_with_rand_data(kmeans, data, dimension, numData, numCenters);
		vl_kmeans_set_centers(kmeans, data, dimension, numCenters);

		// Run at most 100 iterations of cluster refinement using Lloyd algorithm
		vl_kmeans_set_max_num_iterations(kmeans, 100);
		vl_kmeans_refine_centers(kmeans, data, numData);

		// Obtain the energy of the solution
		energy = vl_kmeans_get_energy(kmeans);
		// Obtain the cluster centers
		float *centers = (float*)vl_kmeans_get_centers(kmeans);

		int i, j;
		//for (i = 0; i < numCenters; i++) {
		//	printf("center # %d:\n", i);
		//	for (j = 0; j < dimension; j++) {
		//		printf("    coord1[%d] = %f\n", j, centers[dimension * i + j]);
		//	}
		//}

		//Cantidad de elementos del vector a comparar
		size_t const numDataX = 9;
		//Vector a comparar
		float x[] = {
			3.9f, 3.9f,
			-1.1f, -1.2f,
			-5.9f, 3.9f,
			-5.9f, 3.9f,
			3.9f, 3.9f,
			-1.1f, -1.2f,
			-5.9f, 3.9f,
			-1.1f, -1.2f,
			-1.1f, -1.2f,
		};

		float* result = vl_alldist(dimension, x, numDataX, centers, numCenters);

		int* assignments = new int[numDataX];

		//muestra los valores agrupados por cada elemento espacial de entrada
		for (i = 0; i < numDataX; i++) {
			printf("valores # %d:\n", i);
			float minValue = FLOAT_MAX;
			int cluster = -1;
			for (j = 0; j < numCenters; j++) {
				printf("    result[%d] = %f\n", j, result[numDataX * j + i]);
				if (result[numDataX * j + i] < minValue)
				{
					minValue = result[numDataX * j + i];
					cluster = j;
				}
			}
			assignments[i] = cluster;
		}

		for (int i = 0; i < numDataX; i++)
		{
			cout << "Asignaciones: " << assignments[i] << endl;
		}

		auto hist = computeHistogram(assignments, numDataX, numCenters);
		double sum = 0;
		for (int i = 0; i < numCenters; i++)
		{
			printf("Hist [%d]: -> %f \n", i, hist[i]);
			sum += hist[i];
		}

		printf("Sum: -> %f \n", sum);

		vl_kmeans_delete(kmeans);
	}


	void main(int argc, const char * argv[])
	{
		kmeans_histogram();
	}

	void ikmeans5()
	{
		int row = 255;
		int col = 255;
		Mat show1 = Mat::zeros(row, col, CV_8UC3);
		Mat show2 = show1.clone();

		int data_num = 5000;
		int data_dim = 2;

		vl_uint8 * data = new vl_uint8[data_dim * data_num];

		for (int i = 0; i < data_num; i++)
		{
			vl_uint8 x = data[i * data_dim] = rand() % col;
			vl_uint8 y = data[i * data_dim + 1] = rand() % row;
			circle(show1, Point(x, y), 2, Scalar(255, 255, 255));
		}

		namedWindow("random_points");
		imshow("random_points", show1);
		waitKey(0);

		VlIKMFilt * kmeans = vl_ikm_new(VL_IKM_ELKAN);
		vl_uint k = 3;
		vl_ikm_init_rand(kmeans, data_dim, k);
		vl_ikm_train(kmeans, data, data_num);
		vl_uint * label = new vl_uint[data_num];
		vl_ikm_push(kmeans, label, data, data_num);

		for (int i = 0; i < data_num; i++)
		{
			vl_uint8 x = data[i * data_dim];
			vl_uint8 y = data[i * data_dim + 1];
			switch (label[i])
			{
			case 0:
				circle(show2, Point(x, y), 2, Scalar(255, 0, 0));
				break;
			case 1:
				circle(show2, Point(x, y), 2, Scalar(0, 255, 0));
				break;
			case 2:
				circle(show2, Point(x, y), 2, Scalar(0, 0, 255));
				break;
			}
		}

		const auto *centers = vl_ikm_get_centers(kmeans);
		circle(show2, Point(centers[0], centers[1]), 4, Scalar(255, 255, 0), 4);
		circle(show2, Point(centers[2], centers[3]), 4, Scalar(255, 255, 0), 4);
		circle(show2, Point(centers[4], centers[5]), 4, Scalar(255, 255, 0), 4);

		namedWindow("vlikmeans_result");
		imshow("vlikmeans_result", show2);
		waitKey(0);

		centers = NULL;
		vl_ikm_delete(kmeans);
		delete[] label;
		label = NULL;
		delete[]data;
		data = NULL;
	}

	void kmeans6()
	{
		int data_num = 10;
		int data_dim = 2;
		int k = 2;

		float *data = new float[data_dim * data_num];

		cout << "Points to clustering: " << endl;
		for (int i = 0; i < data_num; i++)
		{
			data[i * data_dim] = (float)rand() / 3.0;
			data[i * data_dim + 1] = (float)rand() / 3.0;
			cout << data[i * data_dim] << "\t" << data[i * data_dim + 1] << endl;
		}

		float * init_centers = new float[data_dim * k];
		cout << "Initial centers: " << endl;
		for (int i = 0; i < k; i++)
		{
			init_centers[i * data_dim] = (float)rand() / 3.0;
			init_centers[i * data_dim + 1] = (float)rand() / 3.0;
			cout << init_centers[i * data_dim] << "\t" << init_centers[i * data_dim + 1] << endl;
		}

		VlKMeans * fkmeans = vl_kmeans_new(VL_TYPE_FLOAT, VlDistanceL2);
		vl_kmeans_set_algorithm(fkmeans, VlKMeansElkan);

		// vl_kmeans_init_centers_plus_plus(fkmeans, data, data_dim, data_num, k);  
		vl_kmeans_set_centers(fkmeans, (void *)init_centers, data_dim, k);
		vl_kmeans_cluster(fkmeans, data, data_dim, data_num, k);
		// vl_kmeans_set_max_num_iterations(fkmeans, 100);  
		// vl_kmeans_refine_centers(fkmeans, data, data_num);  
		// vl_kmeans_cluster(fkmeans, data, data_dim, data_num, k);  

		const float * centers = (float *)vl_kmeans_get_centers(fkmeans);

		cout << "Clustering Centers: " << endl;
		for (int i = 0; i < k; i++)
		{
			cout << centers[i * data_dim] << "\t" << centers[i * data_dim + 1] << endl;
		}
	}

	void kmeans4()
	{
		/*initialize data point*/
		int row = 255;
		int col = 255;
		Mat show = Mat::zeros(row, col, CV_8UC3);
		Mat show2 = show.clone();

		int data_num = 200;
		int data_dim = 2;
		vl_uint8 *data = new vl_uint8[data_num * data_dim];

		for (int i = 0; i < data_num; ++i)
		{
			vl_uint8 x = data[i*data_dim] = rand() % col;
			vl_uint8 y = data[i*data_dim + 1] = rand() % row;
			circle(show, Point(x, y), 2, Scalar(255, 255, 255));
		}

		VlIKMFilt *kmeans = vl_ikm_new(VL_IKM_ELKAN);
		vl_uint K = 3;
		vl_ikm_init_rand(kmeans, data_dim, K);
		vl_ikm_train(kmeans, data, data_num);

		vl_uint * label = new vl_uint[data_num];

		vl_ikm_push(kmeans, label, data, data_num);

		for (int i = 0; i < data_num; ++i)
		{
			vl_uint8 x = data[i*data_dim];
			vl_uint8 y = data[i*data_dim + 1];
			switch (label[i])
			{
			case 0:
				circle(show2, Point(x, y), 2, Scalar(255, 0, 0));
				break;
			case 1:
				circle(show2, Point(x, y), 2, Scalar(0, 255, 0));
				break;
			case 2:
				circle(show2, Point(x, y), 2, Scalar(0, 0, 255));
				break;
			}
		}

		imshow("Imge", show);
		waitKey();

		imshow("Imge2", show2);
		waitKey();

		//imwrite("show.jpg", show);
		//imwrite("show2.jpg", show2);

		vl_ikm_delete(kmeans);

		delete[]label;
		label = NULL;

		delete[]data;
		data = NULL;
	}

	void kmeans3()
	{
		// data generation
		float data[200];
		RNG rng(time(0));
		for (int i = 0; i < 50; i++)
		{
			data[i * 2] = rng.uniform(20, 230);
			data[i * 2 + 1] = rng.uniform(271, 480);
			data[100 + i * 2] = rng.uniform(271, 480);
			data[100 + i * 2 + 1] = rng.uniform(20, 230);
		}

		// kmeans initialization
		VlKMeans *km = vl_kmeans_new(VL_TYPE_FLOAT, VlDistanceL2);

		// parameter setting
		vl_kmeans_set_algorithm(km, VlKMeansAlgorithm::VlKMeansLloyd);
		vl_kmeans_set_initialization(km, VlKMeansInitialization::VlKMeansPlusPlus);

		// kmeans clustering
		vl_kmeans_cluster(km, data, 2, 100, 2);

		// get centers
		const float *centers = (const float *)vl_kmeans_get_centers(km);

		// quanization
		vl_uint32 assignments[100];
		float distances[100];
		vl_kmeans_quantize(km, assignments, distances, data, 100);

		// result visualization
		Mat image(500, 500, CV_8UC3);
		circle(image, Point(centers[0], centers[1]), 3, Scalar(255, 0, 0), 3);
		circle(image, Point(centers[2], centers[3]), 3, Scalar(0, 0, 255), 3);

		for (int i = 0; i < 100; i++)
		{
			if (assignments[i] == 0)
				circle(image, Point(data[i * 2], data[i * 2 + 1]), 1, Scalar(255, 0, 0), 2);
			else
				circle(image, Point(data[i * 2], data[i * 2 + 1]), 1, Scalar(0, 0, 255), 2);
		}

		//cvShowImage("Image", &(IplImage)image);
		imshow("Imge", image);
		waitKey();

		// termination
		vl_kmeans_delete(km);
	}

}