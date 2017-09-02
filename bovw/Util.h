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

#include <fstream>

namespace bovw
{
	class Util
	{
	public:
		Util();
		~Util();		

		/**
		* Carga una matrix de dos dimensiones en un arreglo secuencial en memoria. 
		* @param dir La direcci&oacute;n del archivo a cargar.
		* @param dataSize La cantidad de elementos del a cargar del archivo. si es una imagen dataSize = ancho * alto
		* @return Arreglo secuencial en memoria de los valores cargados
		*/
		static float* loadFloatArray(char* dir, int dataSize);

		/**
		* Carga una matrix de dos dimensiones en un arreglo secuencial en memoria.
		* @param dir La direcci&oacute;n del archivo a cargar.
		* @param dataSize La cantidad de elementos del a cargar del archivo. si es una imagen dataSize = ancho * alto
		* @return Arreglo secuencial en memoria de los valores cargados
		*/
		static int* loadIntegerArray(char* dir, int dataSize);

		static char* pathTo(const char *src);
	};
}
