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

#include "BoVW.h"

namespace bovw
{

	BoVW::BoVW()
	{
		std::cout << "Loading Words!!!" << std::endl;

		centers = Util::loadFloatArray(Util::pathTo("words.tran"), numCenters * phowDescriptorSize);

		lowerChild = Util::loadIntegerArray(Util::pathTo("vocabulary.kdtree.trees.nodes.lowerChild"), numUsedNodes);
		upperChild = Util::loadIntegerArray(Util::pathTo("vocabulary.kdtree.trees.nodes.upperChild"), numUsedNodes);
		splitDimension = Util::loadIntegerArray(Util::pathTo("vocabulary.kdtree.trees.nodes.splitDimension"), numUsedNodes);

		splitThreshold = Util::loadFloatArray(Util::pathTo("vocabulary.kdtree.trees.nodes.splitThreshold"), numUsedNodes);
		lowerBound = Util::loadFloatArray(Util::pathTo("vocabulary.kdtree.trees.nodes.lowerBound"), numUsedNodes);
		upperBound = Util::loadFloatArray(Util::pathTo("vocabulary.kdtree.trees.nodes.upperBound"), numUsedNodes);

		dataIndex = Util::loadIntegerArray(Util::pathTo("vocabulary.kdtree.trees.dataIndex"), numCenters);

		svm_weight = Util::loadFloatArray(Util::pathTo("weight.svm"), weightSize);

		forest = newKdforestFromArray();
		vl_kdforest_set_max_num_comparisons(forest, maxNumComparisons);

		hom = vl_homogeneouskernelmap_new(VlHomogeneousKernelChi2,
			gamma, order, period, VlHomogeneousKernelMapWindowRectangular);
	}

	BoVW::~BoVW()
	{
		delete[] lowerChild;
		delete[] upperChild;
		delete[] splitDimension;

		delete[] splitThreshold;
		delete[] lowerBound;
		delete[] upperBound;

		delete[] dataIndex;

		delete[] svm_weight;

		vl_kdforest_delete(forest);
		vl_homogeneouskernelmap_delete(hom);
	}

	int BoVW::computeAlgorithm2(cv::Mat roi)
	{
		vl_size maxValueBinSize = 10;
		vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };

		int magnif = 6;
		double sigma = 0.8888888f;
		vl_size stride = 4;
		cv::GaussianBlur(roi, roi, cv::Size(magnif + 1, magnif + 1), sigma, sigma, stride);

		cv::cvtColor(roi, roi, CV_BGR2HSV);

		std::vector<cv::Mat> layers;
		cv::split(roi, layers);

		vl_size const image_rows = roi.rows;
		vl_size const image_cols = roi.cols;

		float* p1Data = new float[image_rows * image_cols];
		float* p2Data = new float[image_rows * image_cols];
		float* p3Data = new float[image_rows * image_cols];

		int i, j;
#pragma omp parallel for private(i, j) shared(p1Data, p2Data, p3Data, layers) schedule(runtime)
		for (i = 0; i < image_rows; ++i){
			for (j = 0; j < image_cols; ++j){
				p1Data[image_cols * i + j] = layers[0].at<unsigned char>(i, j);
				p2Data[image_cols * i + j] = layers[1].at<unsigned char>(i, j);
				p3Data[image_cols * i + j] = layers[2].at<unsigned char>(i, j);
			}
		}
		layers.clear();

		size_t channelsCount = 3;
		float** channels = new float*[channelsCount];
		channels[0] = p1Data;
		channels[1] = p2Data;
		channels[2] = p3Data;

		int numKeyPoints = 0;
		int descrSize = 0;
		int binsCount = _countof(binSizes);

		float** histogramsArray = new float*[binsCount];
		int z;

#pragma omp parallel for private(z, numKeyPoints, descrSize) shared(histogramsArray) schedule(runtime)
		for (z = 0; z < binsCount; z++)
		{
			vl_size binSizeActual = binSizes[z];
			float** dsiftPtr = new float*[channelsCount];

			for (int i = 0; i < channelsCount; i++)
			{
				float* dsiftDescriptors = computeDenseSift(channels[i], image_cols, image_rows,
					binSizeActual, maxValueBinSize, &numKeyPoints, &descrSize);

				dsiftPtr[i] = dsiftDescriptors;
			}

			float* phowDescriptors = new float[numKeyPoints * phowDescriptorSize];

			int rr = 0;
			for (int i = 0; i < numKeyPoints; i++)
			{
				for (size_t f = 0, rr = 0; f < channelsCount; f++, rr += descrSize)
				{
					std::copy(dsiftPtr[f] + (descrSize * i),
						dsiftPtr[f] + descrSize * (i + 1),
						phowDescriptors + (phowDescriptorSize * i) + rr);
				}
			}
			deletePtrArray(dsiftPtr, channelsCount);

			vl_uint32* index = new vl_uint32[numNeighbors * numKeyPoints];
			float* distance = new float[numNeighbors * numKeyPoints];

			unsigned int numComparisons = vl_kdforest_query_with_array(forest, index, numNeighbors,
				numKeyPoints, distance, phowDescriptors);

			auto hist = computeHistogram(index, numKeyPoints);
			histogramsArray[z] = hist;

			delete[] index;
			delete[] distance;
			delete[] phowDescriptors;
		}
		deletePtrArray(channels, channelsCount);

		auto histogram = normalizeHistogram(histogramsArray, binsCount);
		deletePtrArray(histogramsArray, binsCount);

		int classT = computeSVMClassifier(histogram);
		delete[] histogram;

		return classT;
	}

	int BoVW::computeAlgorithm(cv::Mat roi)
	{
		vl_size maxValueBinSize = 10;
		vl_size binSizes[] = { 4, 6, 8, maxValueBinSize };

		int magnif = 6;
		double sigma = 0.8888888f;
		vl_size stride = 4;
		cv::GaussianBlur(roi, roi, cv::Size(magnif + 1, magnif + 1), sigma, sigma, stride);

		cv::cvtColor(roi, roi, CV_BGR2HSV);

		std::vector<cv::Mat> layers;
		cv::split(roi, layers);

		vl_size const image_rows = roi.rows;
		vl_size const image_cols = roi.cols;

		float* p1Data = new float[image_rows * image_cols];
		float* p2Data = new float[image_rows * image_cols];
		float* p3Data = new float[image_rows * image_cols];

		int i, j;
#pragma omp parallel for private(i, j) shared(p1Data, p2Data, p3Data, layers) schedule(runtime)
		for (i = 0; i < image_rows; ++i){
			for (j = 0; j < image_cols; ++j){
				p1Data[image_cols * i + j] = layers[0].at<unsigned char>(i, j);
				p2Data[image_cols * i + j] = layers[1].at<unsigned char>(i, j);
				p3Data[image_cols * i + j] = layers[2].at<unsigned char>(i, j);
			}
		}
		layers.clear();

		size_t channelsCount = 3;
		float** channels = new float*[channelsCount];
		channels[0] = p1Data;
		channels[1] = p2Data;
		channels[2] = p3Data;

		int numKeyPoints = 0;
		int descrSize = 0;
		int binsCount = _countof(binSizes);

		float** histogramsArray = new float*[binsCount];
		int z;

#pragma omp parallel for private(z, numKeyPoints, descrSize) shared(histogramsArray) schedule(runtime)
		for (z = 0; z < binsCount; z++)
		{
			vl_size binSizeActual = binSizes[z];
			float** dsiftPtr = new float*[channelsCount];

			for (int i = 0; i < channelsCount; i++)
			{
				float* dsiftDescriptors = computeDenseSift(channels[i], image_cols, image_rows,
					binSizeActual, maxValueBinSize, &numKeyPoints, &descrSize);

				dsiftPtr[i] = dsiftDescriptors;
			}

			float* phowDescriptors = new float[numKeyPoints * phowDescriptorSize];

			int rr = 0;
			for (int i = 0; i < numKeyPoints; i++)
			{
				for (size_t f = 0, rr = 0; f < channelsCount; f++, rr += descrSize)
				{
					std::copy(dsiftPtr[f] + (descrSize * i),
						dsiftPtr[f] + descrSize * (i + 1),
						phowDescriptors + (phowDescriptorSize * i) + rr);
				}
			}
			deletePtrArray(dsiftPtr, channelsCount);

			vl_uint32* index = new vl_uint32[numNeighbors * numKeyPoints];
			float* distance = new float[numNeighbors * numKeyPoints];

			unsigned int numComparisons = vl_kdforest_query_with_array(forest, index, numNeighbors,
				numKeyPoints, distance, phowDescriptors);

			auto hist = computeHistogram(index, numKeyPoints);
			histogramsArray[z] = hist;

			delete[] index;
			delete[] distance;
			delete[] phowDescriptors;
		}
		deletePtrArray(channels, channelsCount);

		auto histogram = normalizeHistogram(histogramsArray, binsCount);
		deletePtrArray(histogramsArray, binsCount);

		int classT = computeSVMClassifier(histogram);
		delete[] histogram;

		return classT;
	}

	float* BoVW::computeDenseSift(float* imgvec, int image_cols, int image_rows, int binSizeActual,
		int maxValueBinSize, int* numKeyPoints, int* descrSize)
	{
		//VlDsiftFilter *dsift = vl_dsift_new_basic(image.cols, image.rows, 4, binSizeActual);
		VlDsiftFilter* dsift = vl_dsift_new(image_cols, image_rows);

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

		//vl_dsift_process(dsift, smoothed);
		vl_dsift_process(dsift, imgvec);

		// echo number of keypoints found
		*numKeyPoints = vl_dsift_get_keypoint_num(dsift);

		// Extract keypoints
		VlDsiftKeypoint const *vlkeypoints;
		vlkeypoints = vl_dsift_get_keypoints(dsift);

		//descriptors dimention
		*descrSize = vl_dsift_get_descriptor_size(dsift);
		//cout << "descrSize: " << descrSize << endl;

		//extract descriptors
		const float* descriptors = vl_dsift_get_descriptors(dsift);

		float contrastthreshold = 0.005f;
		float* descr = new float[*descrSize * *numKeyPoints];
		for (int i = 0; i < *numKeyPoints; i++) {
			for (int j = 0; j < *descrSize; j++) {
				descr[*descrSize * i + j] = vlkeypoints[i].norm < contrastthreshold
					? 0
					: descriptors[*descrSize * i + j] * 512;
			}
		}

		vl_dsift_delete(dsift);
		delete geom;

		return descr;
	}

	VlKDForest* BoVW::newKdforestFromArray()
	{
		VlKDForest* forest;
		VlVectorComparisonType distance = VlDistanceL2;
		vl_size numTrees = 1;

		forest = vl_kdforest_new(VL_TYPE_FLOAT, phowDescriptorSize, numTrees, distance);
		forest->numData = numCenters;
		forest->trees = (VlKDTree**)vl_malloc(sizeof(VlKDTree*)* numTrees);
		forest->data = centers;

		vl_size maxNumNodes = 0;

		vl_uindex ti;
		for (ti = 0; ti < numTrees; ++ti) {
			VlKDTree * tree = (VlKDTree*)vl_malloc(sizeof(VlKDTree));

			maxNumNodes += numUsedNodes;

			tree->numAllocatedNodes = numUsedNodes;
			tree->numUsedNodes = numUsedNodes;
			tree->nodes = (VlKDTreeNode*)vl_malloc(sizeof(VlKDTreeNode)* numUsedNodes);
			tree->dataIndex = (VlKDTreeDataIndexEntry*)vl_malloc(sizeof(VlKDTreeDataIndexEntry)* numCenters);

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
				for (di = 0; di < numCenters; ++di) {
					tree->dataIndex[di].index = dataIndex[di] - 1;
				}
			}

			{
				int numNodesToVisit = tree->numUsedNodes;
				restoreParentRecursively(tree, 0, &numNodesToVisit);
				if (numNodesToVisit != 0) {
					printf("TREE has an inconsitsent tree structure.");
				}
			}

			forest->trees[ti] = tree;
		}

		forest->maxNumNodes = maxNumNodes;

		return forest;
	}

	void BoVW::restoreParentRecursively(VlKDTree *tree, int nodeIndex, int *numNodesToVisit)
	{
		VlKDTreeNode *node = tree->nodes + nodeIndex;
		int lowerChild = node->lowerChild;
		int upperChild = node->upperChild;

		if (*numNodesToVisit == 0) {
			printf("FOREST.TREES has an inconsitsent tree structure.");
		}

		*numNodesToVisit -= 1;

		if (lowerChild >= 0) {
			VlKDTreeNode *child = tree->nodes + lowerChild;
			child->parent = nodeIndex;
			restoreParentRecursively(tree, lowerChild, numNodesToVisit);
		}
		if (upperChild >= 0) {
			VlKDTreeNode *child = tree->nodes + upperChild;
			child->parent = nodeIndex;
			restoreParentRecursively(tree, upperChild, numNodesToVisit);
		}
	}

	float* BoVW::computeHistogram(vl_size* clusterAssignments, int numDataX)
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

	float* BoVW::normalizeHistogram(float** histogramsArray, size_t histogramCount)
	{
		//Construyendo el histograma
		float* histogram = new float[numCenters];
		for (size_t e = 0; e < numCenters; e++)
			histogram[e] = 0;

		for (size_t e = 0; e < histogramCount; e++)
		for (size_t y = 0; y < numCenters; y++)
			histogram[y] += histogramsArray[e][y];

		for (size_t e = 0; e < numCenters; e++)
			histogram[e] /= 4;

		float sum = 0;
		for (int i = 0; i < histogramCount; i++)
			sum += histogram[i];

		//std::cout << "Sum de histograma: " << sum << std::endl << std::endl;

		return histogram;
	}

	int BoVW::computeSVMClassifier(float* histogram)
	{
		float* psi = new float[kernelSize];
		vl_size const resultSize = kernelSize * numCenters;
		float* resultHistogram = new float[resultSize];

		//std::cout << "vl_homogeneouskernelmap_evaluate_f" << std::endl << std::endl;
		for (size_t i = 0; i < numCenters; i++)
		{
			//pass data througt the kernel
			vl_homogeneouskernelmap_evaluate_f(hom, psi, psiStride, histogram[i]);
			//copy result
			for (size_t j = 0; j < kernelSize; j++)
				resultHistogram[(i * kernelSize) + j] = psi[j];
		}

		//std::cout << "svm_weight" << std::endl << std::endl;
		float svm_bias = -1.081984805533894f;

		float svm_threshold = -1.0167452f;  //88.9655% //no -> 95.2676%
		//float svm_threshold = -1.0137452f; //88.2259% //no -> 95.6697%

		//compute linear proyection
		double score = 0;
		for (size_t i = 0; i < weightSize; i++)
			score += svm_weight[i] * resultHistogram[i];
		score += svm_bias;

		int classT = 2 * (score > svm_threshold) - 1;

		//std::cout << "Score: " << score << std::endl << std::endl;
		std::cout << "Score: " << score << " => ";
		//std::cout << "Tipo: " << classT << std::endl << std::endl;

		delete[] psi;
		delete[] resultHistogram;

		return classT;
	}

	void BoVW::deletePtrArray(float** arr, size_t arrCount)
	{
		for (size_t p = 0; p < arrCount; p++)
			delete[] arr[p];
		delete[] arr;
	}
}