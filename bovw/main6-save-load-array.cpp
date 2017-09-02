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

#include <algorithm>
#include <iterator>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

namespace main6
{

	void save1()
	{
		int const count = 10;
		double data[count] = {/*Init Array */ };

		for (int i = 0; i < count; ++i)
		{
			data[i] = i + 0.2f;
		}

		{
			// Write data too a file.
			ofstream outfile("C:/save1.bin", ios::out | ios::binary);
			copy(data, data + count, ostream_iterator<double>(outfile, " "));
		}

		//double data2[count];
		double *data2 = new double[count];

		{
			// Read data from a file
			ifstream infile("C:/save1.bin", ios::in | ios::binary);
			copy(istream_iterator<double>(infile), istream_iterator<double>(), data2);
		}

		for (int i = 0; i < count; ++i)
		{
			cout << data2[i] << endl;
		}
	}

	void save2()
	{
		ofstream ot("C:/save2.bin", ios::out | ios::binary);
		//ifstream in("C:/test.dat",ios::in|ios::binary);

		char** x = new char*[1000];
		for (int i = 0; i < 1000; i++)
			x[i] = new char[1000];

		for (int i = 0; i < 1000; i++)
		for (int j = 0; j < 1000; j++)
			x[i][j] = 10;

		//ot.write(x[0], 1000 * 1000);
		ot.write((char *)&x, sizeof x);
		//ot.write((char *) &x, 1000);

		ot.flush();
		ot.close();
	}

	void save3()
	{
		int nx = 1000, ny = 1000;

		long double *buff1 = new long double[nx * ny];
		long double *buff2 = new long double[nx * ny];

		long double **data = new long double *[nx];
		long double **data_read = new long double *[nx];

		for (int i = 0; i < nx; i++)
		{
			data[i] = buff1 + (i*ny);
			data_read[i] = buff2 + (i*ny);
		}

		data[4][4] = 10.0;
		printf("%f \n", data[4][4]);

		const char* dir = "C:/save3.bin";
		FILE *file;
		errno_t err;
		if ((err = fopen_s(&file, dir, "wb")) != 0) {
			fprintf(stderr, "Cannot open file %s!\n", dir);
		}
		else {
			fwrite(buff1, sizeof(*buff1), nx * ny, file);
			fclose(file);
		}

		if ((err = fopen_s(&file, "C:/save3.bin", "rb")) != 0)
			fprintf(stderr, "Cannot open file %s!\n", dir);
		else{
			fread(buff2, sizeof(*buff2), nx * ny, file);
			fclose(file);
		}

		printf("%f \n", data_read[4][4]);

		// delete pointer arrays
		delete[] data;
		delete[] data_read;

		// delete buffers
		delete[] buff1;
		delete[] buff2;
	}

	//Using a std::vector<> for an RAII Solution

	void save4()
	{
		int nx = 1000, ny = 1000;

		// buffers for allocation
		vector<long double> buff1(nx*ny);

		// holds pointers into original
		vector<long double*> data(nx);

		for (int i = 0; i < nx; i++)
		{
			data[i] = buff1.data() + (i*ny);
		}

		data[4][4] = 19.78;
		cout << data[4][4] << std::endl;

		ofstream ofp("C:/save4.bin", ios::out | ios::binary);
		ofp.write(reinterpret_cast<const char*>(buff1.data()), buff1.size() * sizeof(buff1[0]));
		ofp.close();
	}

	void load4()
	{
		int nx = 1000, ny = 1000;

		// buffers for allocation
		vector<long double> buff2(nx*ny);

		// holds pointers into original	
		vector<long double*> data_read(nx);

		for (int i = 0; i < nx; i++)
		{
			data_read[i] = buff2.data() + (i*ny);
		}

		ifstream ifp("C:/save4.bin", ios::in | ios::binary);
		ifp.read(reinterpret_cast<char*>(buff2.data()), buff2.size() * sizeof(buff2[0]));
		ifp.close();

		cout << data_read[4][4] << endl;
	}

	/**
	* Carga una matrix generada en matlab, la matrix debe ser de doubles
	* el separador debe ser un espacio en blanco
	* dlmwrite('myFile.txt', BoVW.vocabulary.words, 'delimiter', ' ');
	*/
	void loadTxtMatrixFromMatlab()
	{
		int nx = 1000, ny = 384;

		double *buff1 = new double[nx * ny];
		double *buff2 = new double[nx * ny];

		double **data = new double *[nx];
		double **data_read = new double *[nx];

		for (int i = 0; i < nx; i++)
		{
			data[i] = buff1 + (i*ny);
			data_read[i] = buff2 + (i*ny);
		}

		//Cargando el archivo generado en matlab
		char* dir = "C:/words.txt";
		fstream myfile(dir, ios_base::in);
		//i3 input("C:/words2.txt");

		float a;
		int countX = 0, countY = 0;
		while (myfile >> a)
		{
			//printf("Valor: %f \n", a);
			data[countX][countY] = a;

			countX++;
			if (countX = nx)
			{
				countX = 0;
				countY++;
			}
		}
		myfile.close();

		//for (int i = 0; i < nx; i++)
		//for (int j = 0; j < ny; j++)
		//	printf("cargado: [%d][%d] -> %f \n", i, j, data[i][j]);

		//salvando la matrix en un archivo en formato de bytes para su posterior carga
		const char* dirSave = "C:/save3.bin";
		FILE *file;
		errno_t err;
		if ((err = fopen_s(&file, dirSave, "wb")) != 0) {
			fprintf(stderr, "Cannot save file %s!\n", dirSave);
		}
		else {
			fwrite(buff1, sizeof(*buff1), nx * ny, file);
			fclose(file);
		}

		//cargando la matrix a su estado original
		if ((err = fopen_s(&file, dirSave, "rb")) != 0)
			fprintf(stderr, "Cannot open file %s!\n", dir);
		else{
			fread(buff2, sizeof(*buff2), nx * ny, file);
			fclose(file);
		}

		for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			printf("read: [%d][%d] -> %f \n", i, j, data_read[i][j]);
	}

	int main(int argc, const char * argv[])
	{
		printf("Este es el main %d \n", 6);

		//save1();
		//save2();
		//save3();
		//save4();
		//load4();
		
		return argc;
	}
}