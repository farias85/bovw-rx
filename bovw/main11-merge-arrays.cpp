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

#include <iostream>

using namespace std;

namespace main11
{

	int main_arrays() //ESTE ES EL TIPO DE PHOW EN VENTANA DESLIZANTE
	{
		printf("hello \n");

		int const numData = 5;
		int const dimension = 2;

		int arr0[numData * dimension] = {
			1, 1,
			2, 2,
			3, 3,
			1, 1,
			2, 2
		};

		int arr1[numData * dimension] = {
			4, 4,
			5, 5,
			6, 6,
			4, 4,
			5, 5
		};

		int arr2[numData * dimension] = {
			7, 7,
			8, 8,
			9, 9,
			7, 7,
			8, 8
		};

		int const arr3Dimension = 3 * dimension;
		int arr3[numData * arr3Dimension] = { 0 };

		//for (int i = 0; i < numData; i++) {
		//	printf("center # %d:\n", i);
		//	for (int j = 0; j < dimension; j++) {
		//		printf("    coord [%d] = %d\n", j, arr2[dimension * i + j]);
		//	}
		//}

		for (int i = 0; i < numData; i++) {
			//printf("center # %d:\n", i);
			int* aux = new int[arr3Dimension];

			//int k = 0;
			//for (int j = 0; j < dimension; j++, k++) {
			//	//printf("    coord [%d] = %d\n", j, arr2[dimension * i + j]);
			//	aux[k] = arr1[dimension * i + j];
			//}
			//for (int j = 0; j < dimension; j++, k++) {
			//	//printf("    coord [%d] = %d\n", j, arr2[dimension * i + j]);
			//	aux[k] = arr2[dimension * i + j];
			//}

			copy(arr0 + (dimension * i), arr0 + dimension * (i + 1), aux);

			copy(arr1 + (dimension * i), arr1 + dimension * (i + 1), aux + dimension);

			copy(arr2 + (dimension * i), arr2 + dimension * (i + 1), aux + 2 * dimension);

			copy(aux, aux + arr3Dimension, arr3 + (arr3Dimension * i));
		}

		//for (int i = 0; i < numData; i++) {
		//	printf("center # %d:\n", i);
		//	for (int j = 0; j < arr3Dimension; j++) {
		//		printf("    coord [%d] = %d\n", j, arr3[arr3Dimension * i + j]);
		//	}
		//}

		int x = 0;
		int y = 2;
		cout << arr3[arr3Dimension * x + y] << endl;

		int const w = 3;
		int const h = 3;
		int arr4[h * w] = { 0 };
		int desplazamiento = arr3Dimension * y;

		for (int i = 0; i < h; i++) {

			int initY = desplazamiento + arr3Dimension * i;
			int go = initY + x;

			copy(arr3 + go, arr3 + go + w, arr4 + w * i);
		}

		for (int i = 0; i < h; i++) {
			printf("center # %d:\n", i);
			for (int j = 0; j < w; j++) {
				printf("    coord [%d] = %d\n", j, arr4[w * i + j]);
			}
		}

		//int* result = new int[size1 + size2];
		//copy(odd, odd + size1, result);
		//copy(even, even + size2, result + size1);

		system("pause");

		return 0;
	}


	void merge(int *input1, size_t sz1,
		int *input2, size_t sz2,
		int *output, size_t sz3) {
		int i = 0;
		int index1 = 0, index2 = 0;

		while (i < sz3 && index1 < sz1 && index2 < sz2)
		if (input1[index1] <= input2[index2])
			output[i++] = input1[index1++];
		else
			output[i++] = input2[index2++];

		if (index1 < sz1)
		for (; i < sz3 && index1 < sz1; ++i, ++index1)
			output[i] = input1[index1];
		else if (index2 < sz2)
		for (; i < sz3 && index2 < sz2; ++i, ++index2)
			output[i] = input2[index2];
	}

	int main()
	{
		printf("hello \n");

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

		int arr2[numData * dimension] = {
			7, 7,
			8, 8,
			9, 9
		};

		int const arr3Dimension = 3 * dimension;
		int arr3[numData * arr3Dimension] = { 0 };

		//for (int i = 0; i < numData; i++) {
		//	printf("center # %d:\n", i);
		//	for (int j = 0; j < dimension; j++) {
		//		printf("    coord [%d] = %d\n", j, arr2[dimension * i + j]);
		//	}
		//}

		for (int i = 0; i < numData; i++) {
			//printf("center # %d:\n", i);
			int* aux = new int[arr3Dimension];

			//int k = 0;
			//for (int j = 0; j < dimension; j++, k++) {
			//	//printf("    coord [%d] = %d\n", j, arr2[dimension * i + j]);
			//	aux[k] = arr1[dimension * i + j];
			//}
			//for (int j = 0; j < dimension; j++, k++) {
			//	//printf("    coord [%d] = %d\n", j, arr2[dimension * i + j]);
			//	aux[k] = arr2[dimension * i + j];
			//}

			copy(arr0 + (dimension * i), arr0 + dimension * (i + 1), aux);

			copy(arr1 + (dimension * i), arr1 + dimension * (i + 1), aux + dimension);

			copy(arr2 + (dimension * i), arr2 + dimension * (i + 1), aux + 2 * dimension);

			copy(aux, aux + arr3Dimension, arr3 + (arr3Dimension * i));
		}

		for (int i = 0; i < numData; i++) {
			printf("center # %d:\n", i);
			for (int j = 0; j < arr3Dimension; j++) {
				printf("    coord [%d] = %d\n", j, arr3[arr3Dimension * i + j]);
			}
		}

		//int* result = new int[size1 + size2];
		//copy(odd, odd + size1, result);
		//copy(even, even + size2, result + size1);

		system("pause");

		return 0;
	}

	void merge9()
	{
		int odd[] = { 1, 3, 5 };
		int even[] = { 2, 4, 6 };

		int* all[2] = { odd, even };

		// Here's a simple way of doing it... ;-)
		for (unsigned j = 0; j < 3; ++j)
		{
			for (unsigned i = 0; i < 2; i++)
				cout << all[i][j] << " ";
			cout << endl;
		}

		int* aux = *all;

		for (size_t i = 0; i < 6; i++)
		{
			cout << aux[i] << endl;
		}
	}

	void merge8()
	{
		int odd[] = { 1, 3, 5 };
		int even[] = { 2, 4, 6 };

		int size1 = 3;
		int size2 = 3;

		int* result = new int[size1 + size2];
		copy(odd, odd + size1, result);
		copy(even, even + size2, result + size1);

		for (size_t i = 0; i < 6; i++)
		{
			cout << result[i] << endl;
		}

	}

	void merge7()
	{
		int arr1[] = { 1, 2, 3 };
		int arr2[] = { 4, 5, 6 };

		std::basic_string<int> s1(arr1, 3);
		std::basic_string<int> s2(arr2, 3);

		std::basic_string<int> concat(s1 + s2);

		for (std::basic_string<int>::const_iterator i(concat.begin());
			i != concat.end();
			++i)
		{
			std::cout << *i << " ";
		}

		std::cout << std::endl;
	}

	void merge6NO()
	{
		int arraySize = 6;
		int aSize = 3;
		int array1[] = { 1, 2, 3 };
		int array2[] = { 4, 5, 6 };
		int* array3 = new int[arraySize];

		for (int i = 0; i < arraySize * 2; i++)
		{
			if (i < aSize)
			{
				*(array3 + i) = *(array1 + i);
			}
			else if (i >= arraySize)
			{
				*(array3 + i) = *(array2 + (i - arraySize));
			}
		}

		for (size_t i = 0; i < 6; i++)
		{
			cout << array3[i] << endl;
		}
	}

	void merge5()
	{
		int temp[] = { 1, 2, 3, 4 };
		int temp2[] = { 33, 55, 22 };
		int * arr1, *arr2, *arr3;
		int size1(4), size2(3);

		//size1 and size2 is how many elements you 
		//want to copy from the first and second array. In our case all.
		//arr1 = new int[size1]; // optional
		//arr2 = new int[size2];

		arr1 = temp;
		arr2 = temp2;

		arr3 = new int[size1 + size2];
		//or if you know the size: arr3 = new int[size1+size2];

		for (int i = 0; i < size1 + size2; i++){
			if (i < size1)
				arr3[i] = arr1[i];
			else
				arr3[i] = arr2[i - size1];
		}
		cout << endl;
		for (int i = 0; i < size1 + size2; i++) {
			cout << arr3[i] << ", ";
		}
	}

	void merge4()
	{
		int a[] = { 1, 2, 4, 6 };
		int b[] = { 3, 5, 7 };
		int c[100];     // fixme - production code would need a robust way
		//  to ensure c[] always big enough
		int nbr_a = sizeof(a) / sizeof(a[0]);
		int nbr_b = sizeof(b) / sizeof(b[0]);
		int i = 0, j = 0, k = 0;

		// Phase 1) 2 input arrays not exhausted
		while (i < nbr_a && j < nbr_b)
		{
			if (a[i] <= b[j])
				c[k++] = a[i++];
			else
				c[k++] = b[j++];
		}

		// Phase 2) 1 input array not exhausted
		while (i < nbr_a)
			c[k++] = a[i++];
		while (j < nbr_b)
			c[k++] = b[j++];

		for (size_t i = 0; i < nbr_a + nbr_b; i++)
		{
			cout << c[i] << endl;
		}
	}

	void merge3()
	{
		int tab1[] = { 1, 2, 4, 6 };
		int tab2[] = { 3, 5, 7 };

		int tabMerged[_countof(tab1) + _countof(tab2)];
		merge(tab1, _countof(tab1), tab2, _countof(tab2), tabMerged, _countof(tabMerged));

		for (size_t i = 0; i < 7; i++)
		{
			cout << tabMerged[i] << endl;
		}

	}

	void merge2NO()
	{
		double little_arrayA[] = { 1, 2, 3, 4, 5 };
		double little_arrayB[] = { 6, 7, 8, 9, 10 };

		int sizeA = 5;
		int sizeB = 5;

		double* BIG_array = new double[sizeA + sizeB];

		memcpy(BIG_array, little_arrayA, sizeA);
		memcpy(BIG_array + sizeA, little_arrayB, sizeB);

		for (size_t i = 0; i < 10; i++)
		{
			cout << BIG_array[i] << endl;
		}
	}

	void merge1()
	{
		int arSmallOwnerA[] = { 1, 2, 3, 4, 5, 6, 7, 8 };
		int arSmallOwnerB[] = { 1, 2, 3 };
		int* piConcatenation[_countof(arSmallOwnerA) + _countof(arSmallOwnerB) + 1] = { NULL };

		int i = 0;
		for (; i < _countof(arSmallOwnerA); i++)
		{
			piConcatenation[i] = &arSmallOwnerA[i];
		}

		for (int j = 0; i < _countof(arSmallOwnerA) + _countof(arSmallOwnerB); i++, j++)
		{
			piConcatenation[i] = &arSmallOwnerB[j];
		}

		for (int i = 0; i < _countof(piConcatenation) - 1; i++)
		{
			cout << *(piConcatenation[i]) << endl;
		}

		int* aux = *piConcatenation;

		for (size_t i = 0; i < 11; i++)
		{
			cout << aux[i] << endl;
		}
	}

}