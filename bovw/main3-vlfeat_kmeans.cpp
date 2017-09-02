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

extern "C" {
	#include <vl/generic.h>
}
namespace main3
{
	int main()
	{
		printf("Este es el main %d \n", 3);

		double energy;
		
		vl_size const numData = 10;
		vl_size const dimension = 2;

		double data[dimension * numData] = {
			0.0, -0.5,
			0.6, -0.3,
			0.6, -0.3,
			0.6, -0.3,
			0.0, 0.5,
			0.0, 0.5,
			0.0, 0.5,
			0.0, 0.5,
			0.0, 0.5,
			0.6, 0.0 };


		// Use float data and the L2 distance for clustering
		VlKMeans * kmeans = vl_kmeans_new(VL_TYPE_FLOAT, VlDistanceL2);

		// Use Elkan algorithm
		vl_kmeans_set_algorithm(kmeans, VlKMeansElkan);

		// Initialize the cluster centers by randomly sampling the data
		//vl_kmeans_init_centers_with_rand_data (kmeans, data, dimension, numData, numCenters) ;
		int numCenters = 5;
		vl_kmeans_init_centers_with_rand_data(kmeans, data, dimension, numData, numCenters);

		// Run at most 100 iterations of cluster refinement using Lloyd algorithm
		vl_kmeans_set_max_num_iterations(kmeans, 100);
		vl_kmeans_refine_centers(kmeans, data, numData);

		// Obtain the energy of the solution
		energy = vl_kmeans_get_energy(kmeans);
		// Obtain the cluster centers
		const float *centers = (float*)vl_kmeans_get_centers(kmeans);

		int i, j;
		for (i = 0; i < numCenters; i++) {
			printf("center # %d:\n", i);
			for (j = 0; j < dimension; j++) {
				printf("    coord1[%d] = %f\n", j, centers[dimension * i + j]);
			}
		}

		return 0;
	}
}