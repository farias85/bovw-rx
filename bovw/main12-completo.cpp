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

#include "main0-img.cpp"
#include "main1-dsift.cpp"
#include "main2-svm.cpp"
#include "main3-vlfeat_kmeans.cpp"
#include "main4-vlfeat_sift.cpp"
#include "main5-split.cpp"
#include "main6-save-load-array.cpp"
#include "main7-struct.cpp"
#include "main8-lambda.cpp"
#include "main9-compute-histogram.cpp"
#include "main10-homkernelmap.cpp"
#include "main11-merge-arrays.cpp"

using namespace std;
using namespace cv;
using namespace img;

#include <vl/imopv.h>
#include <vl/dsift.h>
#include <vl/kdtree.h>

#include <assert.h>
#include <string.h>

#define FLOAT_MAX 3.40e38

namespace main12
{

	float* loaddata(char* dir, int imgvecSize)
	{
		float value = 0;
		int valueCount = 0;
		float* p1Data = new float[imgvecSize];

		fstream p1(dir, ios_base::in);
		while (p1 >> value)
		{
			p1Data[valueCount] = value;
			valueCount++;
		}
		p1.close();

		return p1Data;
	}

	int* loaddataint(char* dir, int imgvecSize)
	{
		float value = 0;
		int valueCount = 0;
		int* p1Data = new int[imgvecSize];

		fstream p1(dir, ios_base::in);
		while (p1 >> value)
		{
			p1Data[valueCount] = value;
			valueCount++;
		}
		p1.close();

		return p1Data;
	}

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

	float* computeHistogramvlsize(vl_size* clusterAssignments, int numDataX, int numCenters)
	{
		float* histogram = new float[numCenters];

		for (int i = 0; i < numCenters; i++)
			histogram[i] = 0;

		for (int i = 0; i < numDataX; i++)
			histogram[clusterAssignments[i]]++;

		for (int i = 0; i < numCenters; i++)
			histogram[i] /= (double)(numDataX * 4);

		return histogram;
	}

	float* computeHistogram(int* clusterAssignments, int numDataX, int numCenters)
	{
		float* histogram = new float[numCenters];

		for (int i = 0; i < numCenters; i++)
			histogram[i] = 0;

		for (int i = 0; i < numDataX; i++)
			histogram[clusterAssignments[i]]++;

		for (int i = 0; i < numCenters; i++)
			histogram[i] /= (double)(numDataX * 4);

		return histogram;
	}

	void computeSVM(float* histogram, int numCenters)
	{
		vl_size const kernelSize = 5;
		double gamma = 1.0;
		vl_size order = 2;
		double period = -1;
		vl_size psiStride = 1;

		VlHomogeneousKernelMap * hom = vl_homogeneouskernelmap_new(VlHomogeneousKernelChi2,
			gamma, order, period, VlHomogeneousKernelMapWindowRectangular);

		float* psi = new float[kernelSize];
		vl_size const resultSize = kernelSize * numCenters;
		float* resultHistogram = new float[resultSize];

		std::cout << "vl_homogeneouskernelmap_evaluate_f" << endl << endl;
		for (size_t i = 0; i < numCenters; i++)
		{
			//pass data througt the kernel
			vl_homogeneouskernelmap_evaluate_f(hom, psi, psiStride, histogram[i]);
			//copy result
			for (size_t j = 0; j < kernelSize; j++)
				resultHistogram[(i * kernelSize) + j] = psi[j];
		}

		//for (int i = 0; i < resultSize; i++)
		//	cout << "Valor " << i << ": " << result[i] << endl;

		vl_homogeneouskernelmap_delete(hom);

		std::cout << "svm_weight" << endl << endl;
		vl_size const weightSize = 5000;
		float* svm_weight = loaddata("C:/_old/bovw/weight.svm", weightSize);

		float svm_bias = -1.081984805533894f;
		float svm_threshold = -1.0167452f;

		//compute linear proyection
		double score = 0;
		for (size_t i = 0; i < weightSize; i++)
			score += svm_weight[i] * resultHistogram[i];
		score += svm_bias;

		int classT = 2 * (score > svm_threshold) - 1;

		std::cout << "Score: " << score << endl << endl;
		std::cout << "Tipo: " << classT << endl << endl;
	}

	const float* computeDSift2(float* imgvec, int image_cols, int image_rows, int binSizeActual, int maxValueBinSize,
		int* numKeyPoints, int* descrSize)
	{
		VlDsiftFilter *dsift = vl_dsift_new_basic(image_cols, image_rows, 4, binSizeActual);
		//VlDsiftFilter *dsift = vl_dsift_new(image_cols, image_rows);

		//set geometry
		//VlDsiftDescriptorGeometry* geom = new VlDsiftDescriptorGeometry();
		//geom->numBinX = 4;
		//geom->numBinY = 4;
		//geom->numBinT = 8;
		//geom->binSizeX = binSizeActual;
		//geom->binSizeY = binSizeActual;
		//vl_dsift_set_geometry(dsift, geom);

		//set step
		//vl_size step = 4;
		//vl_dsift_set_steps(dsift, step, step);

		//set flatwindows
		//bool useFlatWindow = true;
		//vl_dsift_set_flat_window(dsift, useFlatWindow);

		//set windos size
		vl_dsift_set_window_size(dsift, 1.5f);

		//set the bounds
		int minX, minY, maxX, maxY;
		vl_dsift_get_bounds(dsift, &minX, &minY, &maxX, &maxY);
		int off = floor(1 + 3.0f / 2 * (maxValueBinSize - binSizeActual));
		vl_dsift_set_bounds(dsift, off, off, maxX, maxY);

		////smooth image with gaussian filter
		//int magnif = 6;
		//double sigma = (double)binSizes[0] / (double)magnif;
		//vl_size smoothedStride = 1;
		//vl_size stride = 1;
		//float* smoothed = new float[imgvecSize];
		//vl_imsmooth_f(smoothed, smoothedStride, imgvec, image_cols, image_rows, stride, sigma, sigma);

		//vl_dsift_process(dsift, smoothed);
		vl_dsift_process(dsift, imgvec);

		// echo number of keypoints found
		*numKeyPoints = vl_dsift_get_keypoint_num(dsift);
		std::cout << "Points: " << *numKeyPoints << endl << endl;

		// Extract keypoints
		VlDsiftKeypoint const *vlkeypoints;
		vlkeypoints = vl_dsift_get_keypoints(dsift);

		//descriptors dimention
		*descrSize = vl_dsift_get_descriptor_size(dsift);
		//cout << "descrSize: " << descrSize << endl;

		//extract descriptors
		const float* descriptors = vl_dsift_get_descriptors(dsift);

		//for (int i = 0; i < 1; i++) {
		//	printf("center # %d:\n", i);
		//	for (int j = 0; j < descrSize; j++) {
		//		printf("    coord [%d] = %f\n", j, descriptors[descrSize * i + j]);
		//	}
		//}

		return descriptors;
	}

	float* computeDSift(float* imgvec, int image_cols, int image_rows, int binSizeActual, int maxValueBinSize,
		int* numKeyPoints, int* descrSize)
	{
		//VlDsiftFilter *dsift = vl_dsift_new_basic(image.cols, image.rows, 4, binSizeActual);
		VlDsiftFilter *dsift = vl_dsift_new(image_cols, image_rows);

		//set geometry
		VlDsiftDescriptorGeometry* geom = new VlDsiftDescriptorGeometry();
		geom->numBinX = 4;
		geom->numBinY = 4;
		geom->numBinT = 8;
		geom->binSizeX = binSizeActual;
		geom->binSizeY = binSizeActual;
		vl_dsift_set_geometry(dsift, geom);

		//set step
		vl_size step = 4;
		vl_dsift_set_steps(dsift, step, step);

		//set flatwindows
		bool useFlatWindow = true;
		vl_dsift_set_flat_window(dsift, useFlatWindow);

		//set windos size
		vl_dsift_set_window_size(dsift, 1.5f);

		//set the bounds
		int minX, minY, maxX, maxY;
		vl_dsift_get_bounds(dsift, &minX, &minY, &maxX, &maxY);
		int off = floor(1 + 3.0f / 2 * (maxValueBinSize - binSizeActual)) - 1;
		vl_dsift_set_bounds(dsift, off, off, maxX, maxY);

		////smooth image with gaussian filter
		//int magnif = 6;
		//double sigma = (double)binSizes[0] / (double)magnif;
		//vl_size smoothedStride = 1;
		//vl_size stride = 1;
		//float* smoothed = new float[imgvecSize];
		//vl_imsmooth_f(smoothed, smoothedStride, imgvec, image_cols, image_rows, stride, sigma, sigma);

		//vl_dsift_process(dsift, smoothed);
		vl_dsift_process(dsift, imgvec);

		// echo number of keypoints found
		*numKeyPoints = vl_dsift_get_keypoint_num(dsift);
		std::cout << "Points: " << *numKeyPoints << endl << endl;

		// Extract keypoints
		VlDsiftKeypoint const *vlkeypoints;
		vlkeypoints = vl_dsift_get_keypoints(dsift);

		//descriptors dimention
		*descrSize = vl_dsift_get_descriptor_size(dsift);
		//cout << "descrSize: " << descrSize << endl;

		//extract descriptors
		const float* descriptors = vl_dsift_get_descriptors(dsift);

		//for (int i = 0; i < 1; i++) {
		//	printf("center # %d:\n", i);
		//	for (int j = 0; j < descrSize; j++) {
		//		printf("    coord [%d] = %f\n", j, descriptors[descrSize * i + j]);
		//	}
		//}

		float contrastthreshold = 0.005f;
		float* descr = new float[*descrSize * *numKeyPoints];
		for (int i = 0; i < *numKeyPoints; i++) {
			//printf("center # %d:\n", i);
			for (int j = 0; j < *descrSize; j++) {
				//printf("    coord [%d] = %f\n", j, descriptors[descrSize * i + j] * 512);			
				descr[*descrSize * i + j] = vlkeypoints[i].norm < contrastthreshold
					? 0
					: descriptors[*descrSize * i + j] * 512;
				//descr[*descrSize * i + j] = 0;
			}
		}

		//for (int i = 0; i < 1; i++) {
		//	printf("center # %d:\n", i);
		//	for (int j = 0; j < *descrSize; j++) {
		//		printf("    coord [%d] = %f\n", j, descriptors[*descrSize * i + j]);
		//		printf("    coord [%d] = %f\n", j, descr[*descrSize * i + j]);
		//	}
		//}

		return descr;
	}

	void salvar(char* dir, float* arr, size_t size)
	{
		ofstream fs(dir);

		for (size_t i = 0; i < size; i++)
			fs << arr[i] << endl;

		fs.close();
	}

	void salvar2d(char* dir, float* arr, size_t numData, size_t dimension)
	{
		ofstream fs(dir);

		for (int i = 0; i < numData; i++) {
			for (int j = 0; j < dimension; j++) {
				fs << arr[dimension * i + j] << " ";
			}
			fs << endl;
		}
		fs.close();
	}

	void salvar2dint(char* dir, int* arr, size_t numData, size_t dimension)
	{
		ofstream fs(dir);

		for (int i = 0; i < numData; i++) {
			for (int j = 0; j < dimension; j++) {
				fs << arr[dimension * i + j] << " ";
			}
			fs << endl;
		}
		fs.close();
	}

	void prueba_dsift();
	void main_pincha_ok();
	void main_completo();
	void mergecopy();

	void restore_parent_recursively(VlKDTree * tree, int nodeIndex, int * numNodesToVisit)
	{
		VlKDTreeNode * node = tree->nodes + nodeIndex;
		int lowerChild = node->lowerChild;
		int upperChild = node->upperChild;

		if (*numNodesToVisit == 0) {
			printf("FOREST.TREES has an inconsitsent tree structure.");
		}

		*numNodesToVisit -= 1;

		if (lowerChild >= 0) {
			VlKDTreeNode * child = tree->nodes + lowerChild;
			child->parent = nodeIndex;
			restore_parent_recursively(tree, lowerChild, numNodesToVisit);
		}
		if (upperChild >= 0) {
			VlKDTreeNode * child = tree->nodes + upperChild;
			child->parent = nodeIndex;
			restore_parent_recursively(tree, upperChild, numNodesToVisit);
		}
	}

	VlKDForest* new_kdforest_from_array(float* centers, vl_size numCenters, vl_size dimension)
	{
		VlKDForest* forest;
		VlVectorComparisonType distance = VlDistanceL2;
		vl_size numTrees = 1;

		forest = vl_kdforest_new(VL_TYPE_FLOAT, dimension, numTrees, distance);
		forest->numData = numCenters;
		forest->trees = (VlKDTree**)vl_malloc(sizeof(VlKDTree*)* numTrees);
		forest->data = centers;

		vl_size maxNumNodes = 0;

		vl_uindex ti;
		for (ti = 0; ti < numTrees; ++ti) {
			VlKDTree * tree = (VlKDTree*)vl_malloc(sizeof(VlKDTree));

			//cantidad de valores de los nodos del kdtree
			vl_size numUsedNodes = 1999;
			maxNumNodes += numUsedNodes;

			tree->numAllocatedNodes = numUsedNodes;
			tree->numUsedNodes = numUsedNodes;
			tree->nodes = (VlKDTreeNode*)vl_malloc(sizeof(VlKDTreeNode)* numUsedNodes);
			tree->dataIndex = (VlKDTreeDataIndexEntry*)vl_malloc(sizeof(VlKDTreeDataIndexEntry)* numCenters);

			int* lowerChild = loaddataint("C:/_old/bovw/vocabulary.kdtree.trees.nodes.lowerChild", numUsedNodes);
			int* upperChild = loaddataint("C:/_old/bovw/vocabulary.kdtree.trees.nodes.upperChild", numUsedNodes);
			int* splitDimension = loaddataint("C:/_old/bovw/vocabulary.kdtree.trees.nodes.splitDimension", numUsedNodes);

			float* splitThreshold = loaddata("C:/_old/bovw/vocabulary.kdtree.trees.nodes.splitThreshold", numUsedNodes);
			float* lowerBound = loaddata("C:/_old/bovw/vocabulary.kdtree.trees.nodes.lowerBound", numUsedNodes);
			float* upperBound = loaddata("C:/_old/bovw/vocabulary.kdtree.trees.nodes.upperBound", numUsedNodes);

			{
				vl_uindex ni;
				for (ni = 0; ni < numUsedNodes; ++ni) {
					vl_int32 lc = lowerChild[ni];
					vl_int32 uc = upperChild[ni];
					vl_uint32 d = splitDimension[ni];

					tree->nodes[ni].parent = 0;
					tree->nodes[ni].upperChild = (uc >= 1) ? uc - 1 : uc;
					tree->nodes[ni].lowerChild = (lc >= 1) ? lc - 1 : lc;
					tree->nodes[ni].splitDimension = d - 1;
					tree->nodes[ni].splitThreshold = splitThreshold[ni];
					tree->nodes[ni].lowerBound = lowerBound[ni];
					tree->nodes[ni].upperBound = upperBound[ni];
				}
			}

			{
				vl_uindex di;
				int* dataIndex = loaddataint("C:/_old/bovw/vocabulary.kdtree.trees.dataIndex", numCenters);
				for (di = 0; di < numCenters; ++di) {
					tree->dataIndex[di].index = dataIndex[di] - 1;
				}
			}

			{
				int numNodesToVisit = tree->numUsedNodes;
				restore_parent_recursively(tree, 0, &numNodesToVisit);
				if (numNodesToVisit != 0) {
					printf("TREE has an inconsitsent tree structure.");
				}
			}

			forest->trees[ti] = tree;
		}

		forest->maxNumNodes = maxNumNodes;

		return forest;
	}

	float* normHist(vector<float*> histVec, int numCenters)
	{
		//Construyendo el histograma
		float* histogram = new float[numCenters];
		for (size_t e = 0; e < numCenters; e++)
			histogram[e] = 0;


		for (size_t e = 0; e < histVec.size(); e++)
		for (size_t y = 0; y < numCenters; y++)
			histogram[y] += histVec[e][y];

		for (size_t e = 0; e < numCenters; e++)
			histogram[e] /= 4;


		float sum = 0;
		for (int i = 0; i < histVec.size(); i++)
			sum += histogram[i];

		std::cout << "Sum de histograma: " << sum << endl << endl;

		return histogram;
	}

	void main(int argc, const char * argv[])
	{
		//std::cout << "loaddata - words" << endl << endl;
		//int numCenters = 1000;
		//int dimension = 384;
		////float* centers = loaddata("C:/_old/bovw/words.txt", numCenters * dimension);
		//float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * dimension);

		////numQuerys hace referencia a la cantidad de datos (numData) que se pretenden comparar
		////En el caso de descriptor1.array hay numData = 1, dimension 384, por eso numQueries = 1
		//vl_size numQueries = 911;
		//float* query = loaddata("C:/_old/bovw/descriptor.all.array.tran", dimension * numQueries);

		//VlKDForest* forest = new_kdforest_from_array(centers, numCenters, dimension);
		//vl_size numNeighbors = 1;
		//unsigned int numComparisons = 0;

		//unsigned int maxNumComparisons = 15;
		//vl_kdforest_set_max_num_comparisons(forest, maxNumComparisons);

		//vl_uint32* index = new vl_uint32[numNeighbors * numQueries];
		//float* distance = new float[numNeighbors * numQueries];

		//numComparisons = vl_kdforest_query_with_array(forest, index, numNeighbors, numQueries, distance, query);

		////for (size_t i = 0; i < numQueries; i++)
		////{
		////	cout << i << " ";
		////	cout << "Index: " << index[i] << endl;
		////	cout << i << " ";
		////	cout << "Dist: " << distance[i] << endl;
		////}

		//vl_kdforest_delete(forest);

		//auto histogram = computeHistogramvlsize(index, numQueries, numCenters);

		//salvar("C:/_old/bovw/histogram.last.c++", histogram, numCenters);

		//prueba_dsift();
		//main_pincha_ok();
		main_completo();
		//mergecopy();

		////std::vector<std::vector<float> > ldescriptors = leftImage->descriptors;
		////std::vector<std::vector<float> > rdescriptors = rightImage->descriptors;

		///* KDTree, L1 comparison metric, dimension 128, 1 tree, L1 metric */
		//VlKDForest* forest = vl_kdforest_new(VL_TYPE_FLOAT, 384, 1, VlDistanceL2);

		///* Build the tree from the left descriptors */
		//vl_kdforest_build(forest, numCenters, centers);

		///* Searcher object */
		//VlKDForestSearcher* searcher = vl_kdforest_new_searcher(forest);

		//const vl_size numVecinos = 1000;
		//VlKDForestNeighbor neighbours[numVecinos];

		//int nvisited = vl_kdforestsearcher_query(searcher, neighbours, numVecinos, descriptors384);
		//cout << "Visited: " << nvisited << endl;
		//
		//float minValue = FLOAT_MAX;
		//int cluster = -1;
		//for (int i = 0; i < numVecinos; i++) {			
		//	if (neighbours[i].distance < minValue)
		//	{
		//		cluster = i;
		//		minValue = neighbours[i].distance;
		//	}
		//	cout << i << ": " << neighbours[i].distance << endl;
		//}

		//cout << "Cluster: " << cluster << endl;
		//
		//vl_kdforest_delete(forest);
		//delete[] centers;
		//delete[] descriptors384;

		std::system("pause");
	}

	/*	se comprueban los valores por defecto de configuracion del dsift
	se pudo constatar que los valores del vector de 128 elementos
	son 512 veces menores que los que entrega matlab. Es por eso que se multiplican
	los valores del descriptor por 512. (ver xq sucede esto)
	De igual forma los puntos del frame quedan desplazados por valores que van desde 1 hasta 4 o un poco mas

	La llamada en c++
	VlDsiftFilter *dsift = vl_dsift_new_basic(image_cols, image_rows, 1, 3);
	vl_dsift_process(dsift, tdense);

	es equivalente a [f, d] = vl_dsift(tdense, 'floatdescriptors') ; en matlab
	se obtienen la misma cantidad de puntos
	*/
	void prueba_dsift()
	{
		vl_size const image_rows = 15;
		vl_size const image_cols = 12;

		vl_size imgvecSize = image_rows * image_cols;
		float* p1Data = loaddata("C:/_old/bovw/tdense.im", imgvecSize);

		VlDsiftFilter *dsift = vl_dsift_new_basic(image_cols, image_rows, 1, 3);
		//VlDsiftFilter *dsift = vl_dsift_new(image_cols, image_rows);

		//const VlDsiftDescriptorGeometry* geom = vl_dsift_get_geometry(dsift);
		//cout << geom->binSizeX << endl; //5
		//cout << geom->binSizeY << endl; //5
		//cout << geom->numBinT << endl;  //8
		//cout << geom->numBinX << endl;  //4
		//cout << geom->numBinY << endl;  //4

		//int stepX = 0, stepY = 0;
		//vl_dsift_get_steps(dsift, &stepX, &stepY);
		//cout << stepX << endl; //5
		//cout << stepY << endl; //5

		//int minX, minY, maxX, maxY;
		//vl_dsift_get_bounds(dsift, &minX, &minY, &maxX, &maxY);
		//cout << minX << endl; //0
		//cout << minY << endl; //0
		//cout << maxX << endl; //image_cols - 1
		//cout << maxY << endl; //image_rows - 1

		//cout << vl_dsift_get_flat_window(dsift) << endl; //false -> 0;
		//
		//cout << vl_dsift_get_window_size(dsift) << endl; //2


		vl_dsift_process(dsift, p1Data);

		// echo number of keypoints found
		int numKeyPoints = vl_dsift_get_keypoint_num(dsift);
		std::cout << "Points: " << numKeyPoints << endl << endl;

		// Extract keypoints
		VlDsiftKeypoint const *vlkeypoints;
		vlkeypoints = vl_dsift_get_keypoints(dsift);

		//descriptors dimention
		int descrSize = vl_dsift_get_descriptor_size(dsift);
		//cout << "descrSize: " << descrSize << endl;

		//extract descriptors
		const float* descriptors = vl_dsift_get_descriptors(dsift);

		for (int i = 0; i < 1; i++) {
			printf("center # %d:\n", i);
			for (int j = 0; j < descrSize; j++) {
				printf("    coord [%d] = %f\n", j, descriptors[descrSize * i + j] * 512);
			}
		}

		for (int i = 0; i < numKeyPoints; i++) {
			std::cout << i + 1 << ": ";
			std::cout << "X: " << vlkeypoints[i].x << std::endl;
			std::cout << i + 1 << ": ";
			std::cout << "Y: " << vlkeypoints[i].y << std::endl;
			std::cout << i + 1 << ": ";
			std::cout << "Y: " << vlkeypoints[i].norm << std::endl;
			std::cout << std::endl;
		}
	}

	void main_completo()
	{
		std::cout << "loaddata - words" << endl << endl;
		int numCenters = 1000;
		float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * 384);
		//float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * arr3Dimension);

		//dimensiones de imagen positiva
		//vl_size const image_rows = 79;
		//vl_size const image_cols = 90;

		//vl_size imgvecSize = image_rows * image_cols;

		//float* p1Data = loaddata("C:/_old/bovw/ims.im1", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/ims.im2", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/ims.im3", imgvecSize);

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3", imgvecSize);

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1.tran", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2.tran", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3.tran", imgvecSize);

		//dimensiones de imagen negativa
		vl_size const image_rows = 95;
		vl_size const image_cols = 95;

		vl_size imgvecSize = image_rows * image_cols;

		float* p1Data = loaddata("C:/_old/bovw/pistola.im1.no", imgvecSize);
		float* p2Data = loaddata("C:/_old/bovw/pistola.im2.no", imgvecSize);
		float* p3Data = loaddata("C:/_old/bovw/pistola.im3.no", imgvecSize);

		vector<float*> layer(0);
		layer.push_back(p1Data);
		layer.push_back(p2Data);
		layer.push_back(p3Data);

		vl_size maxValueBinSize = 10;
		vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };

		int numKeyPoints = 0;
		int descrSize = 0;

		vector<float*> histVec(0);

		for (size_t z = 0; z < _countof(binSizes); z++)
		{
			vl_size binSizeActual = binSizes[z];
			vector<float*> dsiftVec(0);

			for (size_t i = 0; i < layer.size(); i++)
			{
				float* dsiftDescriptors = computeDSift(layer[i], image_cols, image_rows,
					binSizeActual, maxValueBinSize, &numKeyPoints, &descrSize);

				//if (z == 0 && i == 0)
				//	salvar2d("C:/_old/bovw/descriptors.all.array.c++", dsiftDescriptors, numKeyPoints, descrSize);

				dsiftVec.push_back(dsiftDescriptors);
			}

			int phowDescriptorSize = 3 * descrSize;
			float* phowDescriptors = new float[numKeyPoints * phowDescriptorSize];

			int rr = 0;
			for (int i = 0; i < numKeyPoints; i++) {

				for (size_t f = 0, rr = 0; f < dsiftVec.size(); f++, rr += descrSize)
				{
					std::copy(dsiftVec[f] + (descrSize * i), dsiftVec[f] + descrSize * (i + 1),
						phowDescriptors + (phowDescriptorSize * i) + rr);
				}
			}

			//if (z == 0)
			//	salvar2d("C:/_old/bovw/descriptors.all.array.c++", phowDescriptors, numKeyPoints, phowDescriptorSize);

			dsiftVec.clear();

			VlKDForest* forest = new_kdforest_from_array(centers, numCenters, phowDescriptorSize);
			vl_size numNeighbors = 1;

			unsigned int maxNumComparisons = 15;
			vl_kdforest_set_max_num_comparisons(forest, maxNumComparisons);

			vl_uint32* index = new vl_uint32[numNeighbors * numKeyPoints];
			float* distance = new float[numNeighbors * numKeyPoints];

			unsigned int numComparisons = vl_kdforest_query_with_array(forest, index, numNeighbors,
				numKeyPoints, distance, phowDescriptors);
			vl_kdforest_delete(forest);

			auto hist = computeHistogramvlsize(index, numKeyPoints, numCenters);

			histVec.push_back(hist);
		}

		auto histogram = normHist(histVec, numCenters);
		//histogram = loaddata("C:/_old/bovw/histogram.txt", numCenters);
		salvar("C:/_old/bovw/histogram.c++", histogram, numCenters);

		computeSVM(histogram, numCenters);
	}

	void main_no_picha_casi_ok()
	{
		std::cout << "loaddata - words" << endl << endl;
		int numCenters = 1000;
		float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * 384);
		//float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * arr3Dimension);

		//dimensiones de imagen positiva
		vl_size const image_rows = 79;
		vl_size const image_cols = 90;

		//dimensiones de imagen negativa
		//vl_size const image_rows = 95;
		//vl_size const image_cols = 95;

		vl_size imgvecSize = image_rows * image_cols;

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3", imgvecSize);

		float* p1Data = loaddata("C:/_old/bovw/ims.im1", imgvecSize);
		float* p2Data = loaddata("C:/_old/bovw/ims.im2", imgvecSize);
		float* p3Data = loaddata("C:/_old/bovw/ims.im3", imgvecSize);

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1.tran", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2.tran", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3.tran", imgvecSize);

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1.no", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2.no", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3.no", imgvecSize);

		vector<float*> layer(0);
		layer.push_back(p1Data);
		layer.push_back(p2Data);
		layer.push_back(p3Data);

		vl_size maxValueBinSize = 10;
		vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };

		int numKeyPoints = 0;
		int descrSize = 0;

		//vector<float*> histXbinSize(0);

		int sumKeyPoints = 0;
		vector<int> sizeVec(0);

		int phowDescriptorSize = 384;
		vector<float*> phowVec(0);

		for (size_t z = 0; z < _countof(binSizes); z++)
		{
			vl_size binSizeActual = binSizes[0];
			vector<float*> dsiftVec(0);

			for (size_t i = 0; i < layer.size(); i++)
			{
				float* dsiftDescriptor = computeDSift(layer[i], image_cols, image_rows,
					binSizeActual, maxValueBinSize, &numKeyPoints, &descrSize);

				dsiftVec.push_back(dsiftDescriptor);

				salvar2d("C:/_old/bovw/descriptors.all.array.c++", dsiftDescriptor, numKeyPoints, descrSize);
			}

			sumKeyPoints += numKeyPoints;
			float* phowDescriptor = new float[numKeyPoints * phowDescriptorSize];

			for (int i = 0; i < numKeyPoints; i++) {

				std::copy(dsiftVec[0] + (descrSize * i), dsiftVec[0] + descrSize * (i + 1),
					phowDescriptor + (phowDescriptorSize * i));

				std::copy(dsiftVec[1] + (descrSize * i), dsiftVec[1] + descrSize * (i + 1),
					phowDescriptor + (phowDescriptorSize * i) + descrSize);

				std::copy(dsiftVec[2] + (descrSize * i), dsiftVec[2] + descrSize * (i + 1),
					phowDescriptor + (phowDescriptorSize * i) + 2 * descrSize);
			}

			dsiftVec.clear();
			phowVec.push_back(phowDescriptor);
			sizeVec.push_back(numKeyPoints * phowDescriptorSize);
		}

		//float* descriptors = new float[sumKeyPoints * phowDescriptorSize];
		//copy(phowVec[0], phowVec[0] + sizeVec[0], descriptors);
		//copy(phowVec[1], phowVec[1] + sizeVec[1], descriptors + sizeVec[0]);
		//copy(phowVec[2], phowVec[2] + sizeVec[2], descriptors + sizeVec[0] + sizeVec[1]);
		//copy(phowVec[3], phowVec[3] + sizeVec[3], descriptors + sizeVec[0] + sizeVec[1] + sizeVec[2]);

		////float* descriptors = phowVec[0];
		////salvar2d("C:/_old/bovw/descriptors.all.array.c++", descriptors, 270, phowDescriptorSize);
		////salvar2d("C:/_old/bovw/descriptors.all.array.c++", descriptors, sumKeyPoints, phowDescriptorSize);

		//VlKDForest* forest = new_kdforest_from_array(centers, numCenters, phowDescriptorSize);
		//vl_size numNeighbors = 1;
		//unsigned int numComparisons = 0;

		//unsigned int maxNumComparisons = 15;
		//vl_kdforest_set_max_num_comparisons(forest, maxNumComparisons);

		//vl_uint32* index = new vl_uint32[numNeighbors * sumKeyPoints];
		//float* distance = new float[numNeighbors * sumKeyPoints];

		//numComparisons = vl_kdforest_query_with_array(forest, index, numNeighbors,
		//	sumKeyPoints, distance, descriptors);

		//vl_kdforest_delete(forest);

		//auto histogram = computeHistogramvlsize(index, numKeyPoints, numCenters);

		////histXbinSize.push_back(hist);

		////Construyendo el histograma
		////float* histogram = new float[numCenters];
		////for (size_t e = 0; e < numCenters; e++)
		////{
		////	histogram[e] = 0;
		////}

		////for (size_t e = 0; e < histXbinSize.size(); e++)
		////{
		////	for (size_t y = 0; y < numCenters; y++)
		////	{
		////		histogram[y] += histXbinSize[e][y];
		////	}
		////}

		////for (size_t e = 0; e < numCenters; e++)
		////{
		////	histogram[e] /= 4;
		////}

		//float sum = 0;
		////vector<float> h(0);	
		//for (int i = 0; i < numCenters; i++)
		//{
		//	//printf("Center [%d]: -> %f \n", i, hist[i]);
		//	//h.push_back(hist[i]);
		//	sum += histogram[i];
		//}

		//std::cout << "Sum de histograma: " << sum << endl << endl;

		////histogram = loaddata("C:/_old/bovw/histogram.txt", numCenters);
		//salvar("C:/_old/bovw/histogram.c++", histogram, numCenters);
		//
		//computeSVM(histogram, numCenters);
	}

	void main_pincha_ok()
	{
		std::cout << "loaddata - words" << endl << endl;
		int numCenters = 1000;
		float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * 384);
		//float* centers = loaddata("C:/_old/bovw/words.tran", numCenters * arr3Dimension);

		//dimensiones de imagen positiva
		vl_size const image_rows = 79;
		vl_size const image_cols = 90;

		//dimensiones de imagen negativa
		//vl_size const image_rows = 95;
		//vl_size const image_cols = 95;

		vl_size imgvecSize = image_rows * image_cols;

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3", imgvecSize);

		float* p1Data = loaddata("C:/_old/bovw/ims.im1", imgvecSize);
		float* p2Data = loaddata("C:/_old/bovw/ims.im2", imgvecSize);
		float* p3Data = loaddata("C:/_old/bovw/ims.im3", imgvecSize);

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1.tran", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2.tran", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3.tran", imgvecSize);

		//float* p1Data = loaddata("C:/_old/bovw/pistola.im1.no", imgvecSize);
		//float* p2Data = loaddata("C:/_old/bovw/pistola.im2.no", imgvecSize);
		//float* p3Data = loaddata("C:/_old/bovw/pistola.im3.no", imgvecSize);

		vector<float*> layer(0);
		layer.push_back(p1Data);
		layer.push_back(p2Data);
		layer.push_back(p3Data);

		vl_size maxValueBinSize = 10;
		vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };

		int numKeyPoints = 0;
		int descrSize = 0;

		vector<float*> histVec(0);

		for (size_t z = 0; z < _countof(binSizes); z++)
		{
			vl_size binSizeActual = binSizes[z];
			vector<float*> dsiftVec(0);

			for (size_t i = 0; i < layer.size(); i++)
			{
				float* dsiftDescriptors = computeDSift(layer[i], image_cols, image_rows,
					binSizeActual, maxValueBinSize, &numKeyPoints, &descrSize);

				dsiftVec.push_back(dsiftDescriptors);
			}

			int phowDescriptorSize = 3 * descrSize;
			float* phowDescriptors = new float[numKeyPoints * phowDescriptorSize];

			int rr = 0;
			for (int i = 0; i < numKeyPoints; i++) {

				for (size_t f = 0, rr = 0; f < dsiftVec.size(); f++, rr += descrSize)
				{
					std::copy(dsiftVec[f] + (descrSize * i), dsiftVec[f] + descrSize * (i + 1),
						phowDescriptors + (phowDescriptorSize * i) + rr);
				}

			}

			dsiftVec.clear();

			VlKDForest* forest = new_kdforest_from_array(centers, numCenters, phowDescriptorSize);
			vl_size numNeighbors = 1;

			unsigned int maxNumComparisons = 15;
			vl_kdforest_set_max_num_comparisons(forest, maxNumComparisons);

			vl_uint32* index = new vl_uint32[numNeighbors * numKeyPoints];
			float* distance = new float[numNeighbors * numKeyPoints];

			unsigned int numComparisons = vl_kdforest_query_with_array(forest, index, numNeighbors,
				numKeyPoints, distance, phowDescriptors);
			vl_kdforest_delete(forest);

			auto hist = computeHistogramvlsize(index, numKeyPoints, numCenters);

			histVec.push_back(hist);
		}

		//Construyendo el histograma
		float* histogram = new float[numCenters];
		for (size_t e = 0; e < numCenters; e++)
		{
			histogram[e] = 0;
		}

		for (size_t e = 0; e < histVec.size(); e++)
		{
			for (size_t y = 0; y < numCenters; y++)
			{
				histogram[y] += histVec[e][y];
			}
		}

		for (size_t e = 0; e < numCenters; e++)
		{
			histogram[e] /= 4;
		}

		float sum = 0;
		//vector<float> h(0);	
		for (int i = 0; i < numCenters; i++)
		{
			//printf("Center [%d]: -> %f \n", i, hist[i]);
			//h.push_back(hist[i]);
			sum += histogram[i];
		}

		std::cout << "Sum de histograma: " << sum << endl << endl;

		//histogram = loaddata("C:/_old/bovw/histogram.txt", numCenters);
		salvar("C:/_old/bovw/histogram.c++", histogram, numCenters);

		computeSVM(histogram, numCenters);
	}

	void mergecopy()
	{
		//main sin opencv2.2, subir a la versión 2.4.8 compatible con vs2013
		int const numData = 3;
		int const dimension = 2;

		int arr0[numData * dimension] = {
			1, 1,
			2, 2,
			3, 3
		};

		int arr1[numData * dimension] = {
			4, 4,
			5, 5,
			6, 6
		};

		//int arr2[numData * dimension] = {
		//	7, 7,
		//	8, 8,
		//	9, 9
		//};

		int* arr2 = new int[numData * dimension];
		arr2[0] = 7;
		arr2[1] = 7;
		arr2[2] = 8;
		arr2[3] = 8;
		arr2[4] = 9;
		arr2[5] = 9;


		int const arr3Dimension = 3 * dimension;
		int arr3Size = numData * arr3Dimension;
		int* arr3 = new int[arr3Size];

		for (int i = 0; i < numData; i++) {

			copy(arr0 + (dimension * i), arr0 + dimension * (i + 1), arr3 + (arr3Dimension * i));

			copy(arr1 + (dimension * i), arr1 + dimension * (i + 1), arr3 + (arr3Dimension * i) + dimension);

			copy(arr2 + (dimension * i), arr2 + dimension * (i + 1), arr3 + (arr3Dimension * i) + 2 * dimension);
		}

		delete[] arr2;

		int numDataArr4 = 2;
		int arr4Size = numDataArr4 * arr3Dimension;
		int* arr4 = new int[arr4Size];
		//data1
		arr4[0] = 2; arr4[1] = 2; arr4[2] = 5; arr4[3] = 5;	arr4[4] = 8; arr4[5] = 8;
		//data2
		arr4[6] = 3; arr4[7] = 3; arr4[8] = 6; arr4[9] = 6;	arr4[10] = 9; arr4[11] = 9;

		int numDataResult = numData + numDataArr4;
		int resultSize = arr3Size + arr4Size;
		int* result = new int[resultSize];

		copy(arr3, arr3 + arr3Size, result);
		copy(arr4, arr4 + arr4Size, result + arr3Size);

		ofstream out("C:/_old/bovw/arrResult.array");
		for (int i = 0; i < numDataResult; i++) {
			printf("center # %d:\n", i);
			for (int j = 0; j < arr3Dimension; j++) {
				printf("    coord [%d] = %d\n", j, result[arr3Dimension * i + j]);
				out << result[arr3Dimension * i + j] << " ";
			}
		}

		out.close();

		salvar2dint("C:/_old/bovw/arrResult.array.2d", result, numDataResult, arr3Dimension);

		//delete[] arr0;
		//delete[] arr1;
	}

	void main_original(int argc, const char * argv[])
	{
		//	//main1::main();
		//	//main2::main();
		//	//main3::main(); 
		//	//main4::main();
		//	//main5::main();
		//	//main6::main(argc, argv);
		//	//main7::main();
		//	//main8::main();	
		//	//main9::main(argc, argv);	
		//	//main10::main();	
		//
		//	Mat load = imread("C:/_old/bovw/pistola.png", CV_LOAD_IMAGE_COLOR);
		//	cout << "W: " << load.cols << endl;
		//	cout << "H: " << load.rows << endl;
		//	cout << "C: " << load.channels() << endl;
		//
		//	cvtColor(load, load, CV_BGR2HSV);
		//	Mat layers[3];
		//	split(load, layers);
		//
		//	Mat image = layers[2];
		//	//Mat image = load;
		//
		//	// transform image in cv::Mat to float vector
		//	vl_size imgvecSize = image.rows * image.cols;
		//	float* imgvec = new float[imgvecSize];
		//	int k = 0;
		//	for (int i = 0; i < image.rows; ++i)
		//	{
		//		for (int j = 0; j < image.cols; ++j)
		//		{
		//			imgvec[k] = image.at<unsigned char>(i, j) / 255.0f;
		//			k++;
		//		}
		//	}
		//
		//	vl_size maxValueBinSize = 10;
		//	vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };
		//
		//	//VlDsiftFilter *dsift = vl_dsift_new_basic(image.cols, image.rows, 4, binSizes[0]);
		//	VlDsiftFilter *dsift = vl_dsift_new(image.cols, image.rows);
		//
		//	//set geometry
		//	VlDsiftDescriptorGeometry* geom = new VlDsiftDescriptorGeometry();
		//	geom->numBinX = 4;
		//	geom->numBinY = 4;
		//	geom->numBinT = 8;
		//	geom->binSizeX = binSizes[0];
		//	geom->binSizeY = binSizes[0];
		//	vl_dsift_set_geometry(dsift, geom);
		//
		//	//set step
		//	vl_size step = 4;
		//	vl_dsift_set_steps(dsift, step, step);
		//
		//	//set flatwindows
		//	bool useFlatWindow = true;
		//	vl_dsift_set_flat_window(dsift, useFlatWindow);
		//
		//	//set windos size
		//	vl_dsift_set_window_size(dsift, 1.5f);
		//
		//	//set the bounds
		//	int minX, minY, maxX, maxY;
		//	vl_dsift_get_bounds(dsift, &minX, &minY, &maxX, &maxY);
		//	int off = floor(1 + 3.0f / 2 * (maxValueBinSize - binSizes[0]));
		//	vl_dsift_set_bounds(dsift, off, off, maxX, maxY);
		//
		//	//smooth image with gaussian filter
		//	int magnif = 6;
		//	double sigma = (double)binSizes[0] / (double)magnif;
		//	vl_size smoothedStride = 1;
		//	vl_size stride = 1;
		//	float* smoothed = new float[imgvecSize];
		//	vl_imsmooth_f(smoothed, smoothedStride, imgvec, image.cols, image.rows, stride, sigma, sigma);	
		//	
		//	//vl_dsift_process(dsift, smoothed);
		//	vl_dsift_process(dsift, imgvec);	
		//	
		//	//namedWindow("AUX");
		//	//imshow("AUX", load);
		//	//waitKey();	
		//	
		//	//vl_size M_ = (image.rows - 1) / step + 1;
		//	//vl_size N_ = (image.cols - 1) / step + 1;
		//
		//	//M_ = M_ % 2 == 0 ? M_ + 1: M_;
		//	//N_ = N_ % 2 == 0 ? N_ + 1 : N_;
		//
		//	//GaussianBlur(load, load, Size(M_, N_), sigma);
		//
		///*	namedWindow("AUX");
		//	imshow("AUX", load);
		//	waitKey();
		//
		//	Mat aux(image.rows, image.cols, CV_32FC1, imgvec);
		//	namedWindow("AUX");
		//	imshow("AUX", aux);
		//	waitKey();	*/	
		//
		//	// echo number of keypoints found
		//	int numKeyPoints = vl_dsift_get_keypoint_num(dsift);
		//	cout << "Points: " << numKeyPoints << endl << endl;
		//
		//	// Extract keypoints
		//	VlDsiftKeypoint const *vlkeypoints;
		//	vlkeypoints = vl_dsift_get_keypoints(dsift);
		//
		//	//for (int i = 0; i < 5; i++) {
		//	//	cout << "X: " << vlkeypoints[i].x << endl;
		//	//	cout << "Y: " << vlkeypoints[i].y << endl;
		//	//	cout << "S: " << vlkeypoints[i].s << endl;
		//	//	cout << std::endl;
		//	//}
		//
		//	//descriptors dimention
		//	int descrSize = vl_dsift_get_descriptor_size(dsift);
		//	cout << "descrSize: " << descrSize << endl;
		//
		//	//extract descriptors
		//	const float* descriptors = vl_dsift_get_descriptors(dsift);
		//	
		//	for (int i = 0; i < 1; i++) {
		//		printf("center # %d:\n", i);
		//		for (int j = 0; j < 128; j++) {
		//			printf("    coord [%d] = %f\n", j, descriptors[128 * i + j]);
		//		}
		//	}
		//
		//	//for (int i = 0; i < numKeyPoints; i++){
		//	//	circle(image, Point(vlkeypoints[i].x, vlkeypoints[i].y),
		//	//		vlkeypoints[i].s,
		//	//		Scalar(255, 0, 0, 0));
		//	//}
		//
		//	//namedWindow("My Image");
		//	//imshow("My Image", image);
		//	//waitKey();
		//
		//	//system("pause");
	}

	void main_otro_dsift()
	{
		////main1::main();
		////main2::main();
		////main3::main(); 
		////main4::main();
		////main5::main();
		////main6::main(argc, argv);
		////main7::main();
		////main8::main();	
		////main9::main(argc, argv);	
		////main10::main();	

		//Mat image = imread("C:/_old/bovw/pistola.png");
		//cout << "W: " << image.cols << endl;
		//cout << "H: " << image.rows << endl;
		//cout << "C: " << image.channels() << endl;

		//// transform image in cv::Mat to float vector
		//vl_size imgvecSize = image.rows * image.cols;
		//float* imgvec = new float[imgvecSize];
		//int k = 0;
		//for (int i = 0; i < image.rows; ++i)
		//{
		//	for (int j = 0; j < image.cols; ++j)
		//	{
		//		imgvec[k] = image.at<unsigned char>(i, j) / 255.0f;
		//		k++;
		//	}
		//}

		//vl_size maxValueBinSize = 10;
		//vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };

		//VlDsiftFilter *vlf = vl_dsift_new_basic(image.cols, image.rows, 4, binSizes[0]);
		////cout << "Windows size: " << vl_dsift_get_window_size(vlf) << endl;
		//vl_dsift_set_window_size(vlf, 1.5f);

		////set the bounds
		//int* minX = new int[1]; int* minY = new int[1];
		//int* maxX = new int[1]; int* maxY = new int[1];
		//vl_dsift_get_bounds(vlf, minX, minY, maxX, maxY);
		//int off = floor(1 + 3.0f / 2 * (maxValueBinSize - binSizes[0]));
		//vl_dsift_set_bounds(vlf, off, off, *maxX, *maxY);

		////smooth image with gaussian filter
		//int magnif = 6;
		//double sigma = (double)binSizes[0] / (double)magnif;
		//vl_size smoothedStride = 1;
		//vl_size stride = 1;
		//float* smoothed = new float[imgvecSize];
		//vl_imsmooth_f(smoothed, smoothedStride, imgvec, image.cols, image.rows, stride, sigma, sigma);

		//vl_dsift_process(vlf, imgvec);
		////vl_dsift_process(vlf, smoothed);

		//// echo number of keypoints found
		//int numKeyPoints = vl_dsift_get_keypoint_num(vlf);
		//cout << "Points: " << numKeyPoints << endl << endl;

		//// Extract keypoints
		//VlDsiftKeypoint const *vlkeypoints;
		//vlkeypoints = vl_dsift_get_keypoints(vlf);

		////for (int i = 0; i < 5; i++) {
		////	cout << "X: " << vlkeypoints[i].x << endl;
		////	cout << "Y: " << vlkeypoints[i].y << endl;
		////	cout << "S: " << vlkeypoints[i].s << endl;
		////	cout << std::endl;
		////}

		////descriptors dimention
		//int dim = vl_dsift_get_descriptor_size(vlf);
		//cout << "Dim: " << dim << endl;

		////extract descriptors
		//const float* descriptors = vl_dsift_get_descriptors(vlf);

		////auto geom = vl_dsift_get_geometry(vlf);

		//int i, j;
		//for (i = 0; i < 1; i++) {
		//	printf("center # %d:\n", i);
		//	for (j = 0; j < 128; j++) {
		//		printf("    coord [%d] = %f\n", j, descriptors[128 * i + j]);
		//	}
		//}

		//for (int i = 0; i < numKeyPoints; i++){
		//	circle(image, Point(vlkeypoints[i].x, vlkeypoints[i].y),
		//		vlkeypoints[i].s,
		//		Scalar(255, 0, 0, 0));
		//}

		//namedWindow("My Image");
		//imshow("My Image", image);
		//waitKey();

		//system("pause");
	}
}