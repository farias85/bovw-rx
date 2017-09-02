#include "HSITable.h"

namespace img
{
	HSITable::HSITable(void)
	{
		loadhCurve();
		loadiCurve();
		loadsMatrix();
	}


	HSITable::~HSITable(void)
	{
		
	}

	void load(char* dir, int nx, int ny, double** data_read, double* array_data)
	{	
		FILE *file;	
		errno_t err;
		if ((err = fopen_s(&file, dir, "rb")) != 0) {
			fprintf(stderr, "Cannot open file %s!\n", dir);
		}
		else {
			fread(array_data, sizeof(*array_data), nx * ny, file);
			fclose(file);
			printf("Cargado... %s! \n", dir);
		}
	}

	void HSITable::loadhCurve()
	{
		char* dir = bovw::Util::pathTo("hCurve.bin");

		const int nx = 1, ny = 256;

		double* array_data = new double[nx*ny];
		double** data = new double*[nx];
		for (int i = 0; i < nx; ++i)
			data[i] = array_data + ny*i;

		load(dir, nx, ny, data, array_data);

		this->hCurve = *data;
	}

	void HSITable::loadiCurve()
	{
		char* dir = bovw::Util::pathTo("iCurve.bin");

		const int nx = 1, ny = 2048;

		double* array_data = new double[nx*ny];
		double** data = new double*[nx];
		for (int i = 0; i < nx; ++i)
			data[i] = array_data + ny*i;

		load(dir, nx, ny, data, array_data);

		this->iCurve = *data;
	}

	void HSITable::loadsMatrix()
	{
		char* dir = bovw::Util::pathTo("sMatrix.bin");

		const int nx = 1, ny = 33 * 2048;

		double* array_data = new double[nx*ny];
		double** data = new double*[nx];
		for (int i = 0; i < nx; ++i)
			data[i] = array_data + ny*i;

		load(dir, nx, ny, data, array_data);

		this->sMatrix = *data;
	}
}