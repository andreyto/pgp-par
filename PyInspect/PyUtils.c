#include "CMemLeak.h"
#include "Python.h"
#include <stdarg.h>
#include "CMemLeak.h"
char PythonErrorString[2048];
extern PyObject* InspectError; // defined in PyInspect.c

FILE* LogFile = NULL;

// Simple interface for error-reporting to Python callers: 
// Print an error to PythonErrorString, then call ReportPythonError().
void ReportPythonError() 
{
    PyErr_SetString(InspectError, PythonErrorString);
}
