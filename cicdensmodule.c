#include "Python.h"

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


static PyObject *cicdens_calcd_cic(PyObject *, PyObject *);


static PyMethodDef dMethods[] = {
    {"calcd_cic",	cicdens_calcd_cic,			METH_VARARGS, "calcd_cic"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initcicdens(void)
{
    PyObject *m;

    m = Py_InitModule("cicdens", dMethods);
    if (m == NULL)
        return;
  //  import_array();
    
}


static PyObject *cicdens_calcd_cic(PyObject *self, PyObject *args)
{
    
	
	
	int ngrid;
	double boxlen;
	long long int xpptr;
	long long int ypptr;
	long long int zpptr;
	long long int dptr;
	double *xpos, *ypos, *zpos;
	int numparts;
    	if (!PyArg_ParseTuple(args, "iidLLLL", &ngrid,&numparts, &boxlen, &xpptr, &ypptr, &zpptr, &dptr))
        	return NULL;
	xpos=(double *)xpptr;
	ypos=(double *)ypptr;
	zpos=(double *)zpptr;
	double *density=(double *)dptr;
	
	
	double cell_len=(boxlen)/((double)ngrid);
	int i;
	int j;
	double dngx, dngy, dngz;
	int cell[3];
	int neighbor[3];
	double posincell[3];
	for (i=0; i<numparts; i++) {
		
		cell[0]=(int)(xpos[i]/cell_len+0.5);
		cell[1]=(int)(ypos[i]/cell_len+0.5);
		cell[2]=(int)(zpos[i]/cell_len+0.5);
		
		posincell[0]=xpos[i]/(cell_len) - (double)cell[0];
		posincell[1]=ypos[i]/(cell_len) - (double)cell[1];
		posincell[2]=zpos[i]/(cell_len) - (double)cell[2];
		
		for (j=0; j<3; j++) {
			neighbor[j]=cell[j];
			if (posincell[j] >0)
				neighbor[j]+=1;
			else
				neighbor[j]-=1;
			
			if (neighbor[j] < 0)
				neighbor[j] = ngrid-1;
			if (neighbor[j] >= ngrid)
				neighbor[j]-=ngrid;
			if (cell[j] >= ngrid)
				cell[j]-=ngrid;
			if (neighbor[j] >= ngrid || neighbor[j] < 0) {
				printf("uh oh out of bounds\n");
			}
			if (cell[j] >=ngrid || cell[j] < 0) {
				printf("uh oh cell out of bounds\n");
			}
			
		}
		dngx=fabs(posincell[0]);
		dngy=fabs(posincell[1]);
		dngz=fabs(posincell[2]);
		
		
		
		density[cell[2] + ngrid * cell[1] + ngrid*ngrid* cell[0]]+=(1-dngx)*(1-dngy)*(1-dngz);
		density[neighbor[2] + ngrid * cell[1] + ngrid*ngrid *cell[0]]+=(1-dngx)*(1-dngy)*(dngz);
		density[neighbor[2]+ngrid*neighbor[1]+ngrid*ngrid *cell[0]]+=(1-dngx)*(dngy)*(dngz);
		density[neighbor[2]+ngrid*cell[1]+ngrid*ngrid*neighbor[0]]+=(dngx)*(1-dngy)*(dngz);
		density[neighbor[2]+ngrid*neighbor[1]+ngrid*ngrid*neighbor[0]]+=(dngx)*(dngy)*(dngz);
		density[cell[2]+ngrid*neighbor[1]+ngrid*ngrid*cell[0]]+=(1-dngx)*(dngy)*(1-dngz);
		density[cell[2]+ngrid*neighbor[1]+ngrid*ngrid*neighbor[0]]+=(dngx)*(dngy)*(1-dngz);
		density[cell[2]+ngrid*cell[1]+ngrid*ngrid*neighbor[0]]+=(dngx)*(1-dngy)*(1-dngz);
        
        
		
	}
	
	return Py_BuildValue("i", 1.0);
}


int main(int argc, char *argv[])
{
    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(argv[0]);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Add a static module */
    initcicdens();
    return 0;
}
