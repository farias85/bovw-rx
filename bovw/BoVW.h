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

#pragma once

#include <vl/dsift.h>
#include <vl/imopv.h>
#include <vl/kdtree.h>
#include <vl/svm.h>

#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>

#include <math.h>
#include <iostream>
#include <vector>
#include <omp.h>

#include "Util.h"

namespace bovw
{
	class BoVW
	{
	private:

		/**
		* cantidad de centros del vocabulario.
		*/
		vl_size numCenters = 1000;		

		/**
		* cantidad de elementos de un vector de características phow.
		*/
		vl_size phowDescriptorSize = 384;

		/**
		* centros del vocabulario. centerSize = numCenters * phowSize
		*/
		float* centers;

		//---------------------------------SVM-----------------------------------------------------
		/**
		* pesos de la svm
		*/
		float* svm_weight;

		/**
		* cantidad de pesos de la svm. Para el Kernel VlHomogeneousKernelChi2
		* se configura con un kernelSize = 5, entonces por cada centro del vocabulario (numCenters)
		* se mapean 5 valores en una dimension superior. weightSize = numCenters * kernelSize
		*/
		vl_size const weightSize = 5000;

		/**
		* tamaño del kernel Chi2
		*/
		vl_size const kernelSize = 5;

		/**
		* ganma
		*/
		double gamma = 1.0;

		/**
		* orden
		*/
		vl_size order = 2;

		/**
		* período
		*/
		double period = -1;

		/**
		* paso de la transformacion
		*/
		vl_size psiStride = 1;

		/**
		* Kernel Map para la transformación de los espacios multidimensionales de la SVM
		*/
		VlHomogeneousKernelMap* hom;

		//---------------------------------Árbol de Búsqueda-----------------------------------------------------
		/**
		* arbol de busqueda
		*/
		VlKDForest* forest;

		/**
		* cantidad de visitantes a chequear por cada nodo en el arbol de busquedas
		*/
		vl_size numNeighbors = 1;

		/**
		* número maximo de comparaciones en un camino del arbol
		*/
		unsigned int maxNumComparisons = 15;

		/**
		* cantidad de valores de los nodos del kdtree.
		*/
		vl_size numUsedNodes = 1999;

		/**
		* vocabulary.kdtree.trees.nodes.lowerChild
		*/
		int* lowerChild;

		/**
		* vocabulary.kdtree.trees.nodes.upperChild
		*/
		int* upperChild;

		/**
		* vocabulary.kdtree.trees.nodes.splitDimension
		*/
		int* splitDimension;

		/**
		* vocabulary.kdtree.trees.nodes.splitThreshold
		*/
		float* splitThreshold;

		/**
		* vocabulary.kdtree.trees.nodes.lowerBound
		*/
		float* lowerBound;

		/**
		* vocabulary.kdtree.trees.nodes.upperBound
		*/
		float* upperBound;

		/**
		* vocabulary.kdtree.trees.dataIndex
		*/
		int* dataIndex;

	public:
		BoVW();
		~BoVW();

		/**
		* Ejecuta el algoritmo BoVW.
		* @param roi Región de Interés donde se ejecutará el algoritmo. Se debe cargar un roi de tres canales RGB.		
		* @return Devuelve 1 si es pistola, -1 en caso contrario.
		*/
		int computeAlgorithm(cv::Mat roi);
		int computeAlgorithm2(cv::Mat roi);

	private:
		/**
		* Ejecuta el algoritmo Dense Sift configurado para ejcutarse en los 3 canales de la imagen.
		* @param imgvec Los datos de la imagen.
		* @param image_cols Cantidad de columnas.
		* @param image_rows Cantidad de filas.
		* @param binSizeActual Tamanno actual del espaciado de la maya de los descriptores phow.
		* @param maxValueBinSize Tamanno mayor del espaciado de la maya de los descriptores phow.
		* @param numKeyPoints Parámetro de salida, numero de puntos de interés encontrados en la maya.
		* @param descrSize Parámetro de salida, numero de elementos del descriptor.
		* @return Descriptores Dense Sift encontrados.
		*/
		float* computeDenseSift(float* imgvec, int image_cols, int image_rows, int binSizeActual,
			int maxValueBinSize, int* numKeyPoints, int* descrSize);

		/**
		* Crea un nuevo KdForest, para acelerar las busquedas en la asignación de una palabra visual
		* a un cluster específico, con lo cual se crear el histograma de palabras visuales				
		* @return Universo de árboles de búsquedas configurado.
		*/
		VlKDForest* newKdforestFromArray();

		/**
		* Restaura el nodo padre luego de una búsqueda por un camino sin salida satisfactoria.
		* @param tree Árbol de búsqueda.
		* @param nodeIndex Nodo inicial.
		* @param numNodesToVisit Cantidad de nodos a visitar.				
		*/
		void restoreParentRecursively(VlKDTree *tree, int nodeIndex, int *numNodesToVisit);

		/**
		* Calcula el histograma de palabras visuales.
		* @param clusterAssignments Lista de las asignaciones a los clusters de los vectores de características.
		* @param numDataX Dimensión del vector de asignaciones.		
		* @return Histograma de palabras visuales.
		*/
		float* computeHistogram(vl_size* clusterAssignments, int numDataX);

		/**
		* Normaliza el histograma de palabras visuales.
		* @param histogramsArray Hitograma de palabras visuales.
		* @param histogramCount Cantidad de histogramas contenidos en el histogramArray.
		* @return Histograma de palabras visuales normalizado.
		*/
		float* normalizeHistogram(float** histogramsArray, size_t histogramCount);

		/**
		* Ejecuta la clasificación de un histograma a través de una SVM.
		* @param histogram Histograma normalizado.
		* @return Devuelve 1 si es pistola, -1 en caso contrario.
		*/
		int computeSVMClassifier(float* histogram);

		/**
		* Libera la memoria acumulada por los punteros.
		* @param arr El arreglo de punteros a float.
		* @param arrCount Cantidad de elementos del arreglo.		
		*/
		void deletePtrArray(float** arr, size_t arrCount);
	};
}
